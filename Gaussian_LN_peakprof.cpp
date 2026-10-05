#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <charconv>
#include <cmath>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "fft.hpp"

// ============================================================================
// This file adds a lognormal-power-spectrum version of Gaussian_mono_peak.cpp
// (which uses a Dirac-delta/monochromatic power spectrum -- a single Fourier
// shell). Almost everything below -- the flat Grid/RGrid layout, the cached
// in-place FFTW3DPlan, the fast valuesToCSV/appendNum writer, the Laplacian +
// gradient + Lpeak section, and the sliding-window compaction + Cpeak section
// -- is ported unchanged from the already-optimized Gaussian_mono_peak.cpp;
// see that file's own header comment for the rationale behind those. The one
// genuinely new piece is buildGk() below, needed because a lognormal spectrum
// requires summing contributions from *every* Fourier shell (wavenumber =
// 1, 2, 3, ... up to the box's Nyquist corner), not just one:
//
// 8. The original (pre-optimization) Gaussian_LN_peak.cpp called a per-shell
//    dwk(wavenumber, engine) once for every wavenumber, and each call did
//    THREE full NL^3 sweeps (independent draws, reflection, normalization)
//    just to reach the few thousand points that are actually in that one
//    shell. With ~1.7*NL wavenumbers summed, that is O(NL^4) work overall --
//    the dominant cost of the whole program once NL gets large. buildGk()
//    below instead does a SINGLE NL^3 sweep that classifies every grid point
//    into its shell index once, bucket-sorts (counting sort) the point
//    indices by shell into one flat `order` array with each shell occupying
//    a contiguous range, and then loops over shells 1..maxW processing only
//    the points in that shell's slice. Total cost is O(NL^3) instead of
//    O(NL^4). Because dn = 1 makes shells disjoint (every grid point belongs
//    to at most one shell), each shell's slice of `order` is exactly the
//    "count" the original per-shell dwk() computed internally, and each
//    shell's contribution is written directly into its own (guaranteed
//    still-zero) slots of the persistent gk array -- no separate per-shell
//    array plus "+=" is needed, only a direct assignment.
//    RNG determinism is preserved exactly as in the original: shells are
//    still visited strictly in increasing wavenumber order (the shared
//    std::mt19937 engine advances across shells exactly as it did across the
//    original's sequential dwk() calls), and within a shell the independent
//    draws are still made in strictly increasing (i,j,k) row-major order
//    (the bucketing sweep that builds `order` runs single-threaded, in
//    increasing linear-index order, specifically so this holds). Only the
//    RNG-free parts (shell classification/histogram, reflection, and the
//    per-shell normalize+weight step) are parallelised.
//
// The FFTW plan is constructed *after* buildGk() returns (unlike
// Gaussian_mono_peak.cpp, which builds it first) so that buildGk()'s own
// transient bucketing arrays (freed automatically when it returns) are never
// alive at the same time as the plan's NL^3 scratch buffer.
// ============================================================================

Grid buildGk(int seed);
double WRTH(double z);
double powerspectrum(int wavenumber);
int shiftedindex(int n);                                  // shifted index
bool realpoint(int nx, int ny, int nz);                    // judge real point
bool complexpoint(int nx, int ny, int nz);                 // judge independent complex point

// random distribution
std::normal_distribution<> dist(0., 1.);

// imaginary unit
const std::complex<double> II(0, 1);

// real-space top-hat window
double WRTH(double z)
{
  if (z == 0)
  {
    return 1;
  }
  else
  {
    return 3 * (sin(z) - z * cos(z)) / pow(z, 3);
  }
}

// parameters
const int NL = 256;
const int nsigma = 16; //32;
const double s2 = 0.1; //0.2; // 0.01; 
const std::string s2value = "0,1"; // "0,2"; // "0,01";
const double As = 5e-3;
const double dn = 1; // Thickness of nsigma sphere shell
const std::string mapfileprefix = std::string("direct_256/LN_map_") + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_";
const std::string laplacianfileprefix = std::string("direct_256/LN_laplacian_") + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_";
const std::string Lpeakfileprefix = std::string("direct_256/LN_Lpeak_") + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_";
const std::string Cpeakfileprefix = std::string("direct_256/LN_Cpeak_") + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_";

// lognormal power spectrum, peaked around wavenumber = nsigma
double powerspectrum(int wavenumber)
{
  return exp(-pow(log(static_cast<double>(wavenumber)) - log(static_cast<double>(nsigma)), 2) / 2 / s2) / sqrt(2 * M_PI * s2);
}

// weight applied to shell `w`'s (unit-normalized) contribution before adding
// it into gk -- identical formula to the original file's main(): the i=1
// shell is special-cased (no "/i" factor), matching it exactly.
inline double gkWeight(int w)
{
  if (w == 1)
    return sqrt(powerspectrum(1) * dn);
  else
    return sqrt(powerspectrum(w) * dn / w);
}

// flat-index helper: matches the (i,j,k) -> i*NL*NL + j*NL + k layout used
// throughout (row-major, same iteration order as the original triple loop)
inline size_t IDX(int i, int j, int k)
{
  return (static_cast<size_t>(i) * NL + j) * NL + k;
}

inline int nextIndex(int n) { return (n == NL - 1) ? 0 : n + 1; }

// Shell index of grid point (i,j,k): the unique integer w such that
// w - dn/2 <= |shiftedindex(i,j,k)| < w + dn/2 -- i.e. exactly the condition
// the original innsigma(i,j,k,wavenumber) tested, but computed directly
// instead of checked against one candidate wavenumber at a time.
inline int computeShell(int i, int j, int k)
{
  double nxt = shiftedindex(i);
  double nyt = shiftedindex(j);
  double nzt = shiftedindex(k);
  double ntnorm = sqrt(nxt * nxt + nyt * nyt + nzt * nzt);
  return static_cast<int>(std::floor(ntnorm / dn + 0.5));
}

// Convert n values (accessed through getValue(idx)) to a comma-separated
// CSV string, using std::to_chars in parallel chunks. Produces byte-for-byte
// the same text as the original `stream << value` (general format, 6
// significant digits -- verified), but without paying iostream formatting
// overhead 16.7 million times over.
template <typename F>
std::string valuesToCSV(size_t n, F &&getValue)
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  if (nthreads < 1) nthreads = 1;

  std::vector<std::string> parts(nthreads);
  size_t chunk = (n + nthreads - 1) / static_cast<size_t>(nthreads);

#pragma omp parallel for schedule(static)
  for (int t = 0; t < nthreads; t++)
  {
    size_t begin = std::min(n, static_cast<size_t>(t) * chunk);
    size_t end = std::min(n, begin + chunk);
    std::string &s = parts[t];
    s.reserve((end - begin) * 14);
    char tmp[32];
    for (size_t idx = begin; idx < end; idx++)
    {
      auto res = std::to_chars(tmp, tmp + sizeof(tmp), getValue(idx), std::chars_format::general, 6);
      s.append(tmp, res.ptr - tmp);
      if (idx + 1 != n) s.push_back(',');
    }
  }

  size_t total = 0;
  for (auto &s : parts) total += s.size();
  std::string out;
  out.reserve(total);
  for (auto &s : parts) out += s;
  return out;
}

// Appends v to s using the same text format as `stream << v` (general
// format, 6 significant digits for doubles; plain decimal for ints) --
// verified byte-for-byte identical to ostream's default formatting.
inline void appendNum(std::string &s, double v)
{
  char tmp[32];
  auto res = std::to_chars(tmp, tmp + sizeof(tmp), v, std::chars_format::general, 6);
  s.append(tmp, res.ptr - tmp);
}
inline void appendNum(std::string &s, int v)
{
  char tmp[16];
  auto res = std::to_chars(tmp, tmp + sizeof(tmp), v);
  s.append(tmp, res.ptr - tmp);
}

// Computes the compaction field (real space) smoothed at radius rs, given
// the un-smoothed Fourier-space map gk. Identical to Gaussian_mono_peak.cpp's
// version -- it only depends on gk and the plan, not on how gk was built.
RGrid computeCompactionShell(const Grid &gk, const FFTW3DPlan &fftplan, int rs)
{
  Grid rzpk(gk.size());
#pragma omp parallel for collapse(2) schedule(static)
  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
      {
        size_t idx = IDX(i, j, k);
        int nxt = shiftedindex(i);
        int nyt = shiftedindex(j);
        int nzt = shiftedindex(k);
        double ntnorm = sqrt(nxt * nxt + nyt * nyt + nzt * nzt);
        double kr = 2 * M_PI * ntnorm * rs / NL;
        rzpk[idx] = gk[idx] * (-kr * kr / 3 * WRTH(kr) * sqrt(As));
      }

  Grid rzpx = fftplan.execute(rzpk);
  rzpk = Grid(); // free -- not needed once its transform is in hand

  RGrid compaction(gk.size());
#pragma omp parallel for collapse(2) schedule(static)
  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
      {
        size_t idx = IDX(i, j, k);
        double v = rzpx[idx].real();
        compaction[idx] = 2. / 3 * (1 - (1 + v) * (1 + v));
      }

  return compaction;
}

// Peak test across the radius direction combined with the spatial
// saddle-point test on C0's gradient. Identical to Gaussian_mono_peak.cpp's
// version.
void processCompactionPeaks(int r, const RGrid &C0, const RGrid &C1, const RGrid &C2, std::ofstream &Cpeakfile)
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  std::vector<std::string> buffers(nthreads);

#pragma omp parallel
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::string &buf = buffers[tid];
#pragma omp for collapse(2) schedule(static)
    for (int i = 0; i < NL; i++)
      for (int j = 0; j < NL; j++)
        for (int k = 0; k < NL; k++)
        {
          size_t idx = IDX(i, j, k);
          if (C0[idx] < C1[idx] && C1[idx] > C2[idx])
          {
            int ip = nextIndex(i), jp = nextIndex(j), kp = nextIndex(k);
            int ip2 = nextIndex(ip), jp2 = nextIndex(jp), kp2 = nextIndex(kp);
            double dx = C0[IDX(ip, j, k)] - C0[idx];
            double dxrot = C0[IDX(ip2, j, k)] - C0[IDX(ip, j, k)];
            double dy = C0[IDX(i, jp, k)] - C0[idx];
            double dyrot = C0[IDX(i, jp2, k)] - C0[IDX(i, jp, k)];
            double dz = C0[IDX(i, j, kp)] - C0[idx];
            double dzrot = C0[IDX(i, j, kp2)] - C0[IDX(i, j, kp)];
            if (dx * dxrot < 0 && dx > 0 &&
                dy * dyrot < 0 && dy > 0 &&
                dz * dzrot < 0 && dz > 0)
            {
              int rp = r + 1;
              double val = C1[IDX(ip, jp, kp)];
              appendNum(buf, rp); buf += ',';
              appendNum(buf, ip); buf += ',';
              appendNum(buf, jp); buf += ',';
              appendNum(buf, kp); buf += ',';
              appendNum(buf, val); buf += '\n';
            }
          }
        }
  }

  for (auto &buf : buffers) Cpeakfile << buf;
}


// ============================================================================
// Added main(): for every Laplacian peak with nu2 >= nuCut, evaluate the
// compaction profile C(r) (r = 1 .. Rbound, integer smoothing radius in grid
// units) at the peak position, and write per-peak
//   ip,jp,kp, mu2, k3, Cmax_c, rm_c, Cmax_n, rm_n, C(1),...,C(Rb), zeta(0),zeta(1),...,zeta(Rb)
// where mu2 = Laplacian value at the peak (= nu2*sigma2), k3 = sqrt(D2D2/D2)
// (same two columns as the Lpeak file), Cmax_c / rm_c = maximum over r of the
// compaction AT the peak grid point and its radius, Cmax_n / rm_n = same but
// taking the maximum over the 3x3x3 neighbourhood of the peak at each r, and
// C(r) = compaction at the peak grid point for each r (profile); zeta(r) = spherical
// shell average of zeta (= sqrt(As) g) at radius r about the peak (r=0: point value).
// Output: direct_256/LN_Cprof_<s2value>_<NL>_<nsigma>_<seed>.csv
// ============================================================================
// Spherical-shell-averaged zeta at smoothing radius rs (grid units): kernel sin(kr)/(kr), times sqrt(As).
// (v = r d(zeta)/dr of this field is exactly the quantity used in computeCompactionShell.)
RGrid computeZetaShell(const Grid &gk, const FFTW3DPlan &fftplan, int rs)
{
  Grid zk(gk.size());
#pragma omp parallel for collapse(2) schedule(static)
  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
      {
        size_t idx = IDX(i, j, k);
        int nxt = shiftedindex(i), nyt = shiftedindex(j), nzt = shiftedindex(k);
        double ntnorm = sqrt(nxt * nxt + nyt * nyt + nzt * nzt);
        double kr = 2 * M_PI * ntnorm * rs / NL;
        double w = (kr == 0) ? 1. : sin(kr) / kr;
        zk[idx] = gk[idx] * (w * sqrt(As));
      }
  Grid zx = fftplan.execute(zk);
  RGrid z(gk.size());
#pragma omp parallel for collapse(2) schedule(static)
  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
        z[IDX(i, j, k)] = zx[IDX(i, j, k)].real();
  return z;
}

const double nuCut = 3.; // keep peaks with nu2 = mu2/sigma2 >= nuCut

struct PeakRec { int i, j, k; double mu2, k3; };

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cerr << "Specify the noise file number correctly." << std::endl;
    return 1;
  }
  struct timeval Nv; struct timezone Nz;
  gettimeofday(&Nv, &Nz);
  double before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;

  int seed = atoi(argv[1]);
  const std::string profprefix = std::string("direct_256/LN_Cprof_") + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_";
  std::ofstream proffile(profprefix + std::to_string(seed) + ".csv");

  Grid gk = buildGk(seed);
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  FFTW3DPlan fftplan(NL, nthreads);

  // nu2 = mu2 / sigma2 with sigma2 = sqrt(exp(8 s2) kpeak^4) (theory; Laplacian of gk, no As factor)
  const double kpk = 2 * M_PI * nsigma / NL;
  const double sigma2th = sqrt(exp(2 * 4 * s2) * pow(kpk, 4));
  std::vector<PeakRec> peaks;
  {
    Grid D2gk(gk.size()), D2D2gk(gk.size());
#pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < NL; i++)
      for (int j = 0; j < NL; j++)
        for (int k = 0; k < NL; k++)
        {
          size_t idx = IDX(i, j, k);
          int nxt = shiftedindex(i), nyt = shiftedindex(j), nzt = shiftedindex(k);
          double ntnorm = sqrt(nxt * nxt + nyt * nyt + nzt * nzt);
          D2gk[idx] = gk[idx] * pow(2 * M_PI * ntnorm / NL, 2);
          D2D2gk[idx] = gk[idx] * pow(2 * M_PI * ntnorm / NL, 4);
        }
    Grid D2gx = fftplan.execute(D2gk);
    D2gk = Grid();
    Grid D2D2gx = fftplan.execute(D2D2gk);
    D2D2gk = Grid();

    // single-threaded scan (peak order deterministic); same peak test as Gaussian_LN_peak
    for (int i = 0; i < NL; i++)
      for (int j = 0; j < NL; j++)
        for (int k = 0; k < NL; k++)
        {
          size_t idx = IDX(i, j, k);
          int ip = nextIndex(i), jp = nextIndex(j), kp = nextIndex(k);
          int ip2 = nextIndex(ip), jp2 = nextIndex(jp), kp2 = nextIndex(kp);
          double dx = D2gx[IDX(ip, j, k)].real() - D2gx[idx].real();
          double dxrot = D2gx[IDX(ip2, j, k)].real() - D2gx[IDX(ip, j, k)].real();
          double dy = D2gx[IDX(i, jp, k)].real() - D2gx[idx].real();
          double dyrot = D2gx[IDX(i, jp2, k)].real() - D2gx[IDX(i, jp, k)].real();
          double dz = D2gx[IDX(i, j, kp)].real() - D2gx[idx].real();
          double dzrot = D2gx[IDX(i, j, kp2)].real() - D2gx[IDX(i, j, kp)].real();
          if (dx * dxrot < 0 && dx > 0 && dy * dyrot < 0 && dy > 0 && dz * dzrot < 0 && dz > 0)
          {
            size_t idxp = IDX(ip, jp, kp);
            double mu2 = D2gx[idxp].real();
            if (mu2 / sigma2th >= nuCut)
              peaks.push_back({ip, jp, kp, mu2, sqrt(D2D2gx[idxp].real() / mu2)});
          }
        }
  }
  std::cout << "peaks kept: " << peaks.size() << std::endl;

  double Rbound = 5. / (2 * M_PI * nsigma / NL);
  int Rb = (int)Rbound; // radii rs = 1..Rb
  size_t np = peaks.size();
  std::vector<double> prof(np * Rb), Cn(np * Rb), zprof(np * (Rb + 1));
  {
    RGrid Z0 = computeZetaShell(gk, fftplan, 0);
    for (size_t p = 0; p < np; p++) zprof[p * (Rb + 1)] = Z0[IDX(peaks[p].i, peaks[p].j, peaks[p].k)];
  }
  for (int rs = 1; rs <= Rb; rs++)
  {
    RGrid C = computeCompactionShell(gk, fftplan, rs);
    {
      RGrid Z = computeZetaShell(gk, fftplan, rs);
      for (size_t p = 0; p < np; p++) zprof[p * (Rb + 1) + rs] = Z[IDX(peaks[p].i, peaks[p].j, peaks[p].k)];
    }
#pragma omp parallel for schedule(static)
    for (size_t p = 0; p < np; p++)
    {
      const PeakRec &pk = peaks[p];
      prof[p * Rb + (rs - 1)] = C[IDX(pk.i, pk.j, pk.k)];
      double m = -1e300;
      for (int a = -1; a <= 1; a++)
        for (int b = -1; b <= 1; b++)
          for (int c = -1; c <= 1; c++)
          {
            int ii = (pk.i + a + NL) % NL, jj = (pk.j + b + NL) % NL, kk = (pk.k + c + NL) % NL;
            m = std::max(m, C[IDX(ii, jj, kk)]);
          }
      Cn[p * Rb + (rs - 1)] = m;
    }
  }

  std::string buf;
  for (size_t p = 0; p < np; p++)
  {
    const PeakRec &pk = peaks[p];
    int rc = 0, rn = 0;
    for (int r = 1; r < Rb; r++)
    {
      if (prof[p * Rb + r] > prof[p * Rb + rc]) rc = r;
      if (Cn[p * Rb + r] > Cn[p * Rb + rn]) rn = r;
    }
    appendNum(buf, pk.i); buf += ','; appendNum(buf, pk.j); buf += ','; appendNum(buf, pk.k); buf += ',';
    appendNum(buf, pk.mu2); buf += ','; appendNum(buf, pk.k3); buf += ',';
    appendNum(buf, prof[p * Rb + rc]); buf += ','; appendNum(buf, rc + 1); buf += ',';
    appendNum(buf, Cn[p * Rb + rn]); buf += ','; appendNum(buf, rn + 1);
    for (int r = 0; r < Rb; r++) { buf += ','; appendNum(buf, prof[p * Rb + r]); }
    for (int r = 0; r <= Rb; r++) { buf += ','; appendNum(buf, zprof[p * (Rb + 1) + r]); }
    buf += '\n';
  }
  proffile << buf;
  std::cout << "Exported to " << profprefix + std::to_string(seed) + ".csv" << std::endl;

  gettimeofday(&Nv, &Nz);
  std::cout << (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6 - before << " sec." << std::endl;
  return 0;
}


// -----------------------------------------------

// Builds gk = sum over wavenumber w of dwk_w * weight(w), where dwk_w is a
// unit-normalized random draw over shell w (the same lognormal-map
// construction as the original file's main() loop), but in O(NL^3) total
// instead of O(NL^4) -- see the header comment for the algorithm.
Grid buildGk(int seed)
{
  const size_t N = static_cast<size_t>(NL) * NL * NL;
  Grid gk(N, std::complex<double>(0, 0));

  // Upper bound on any reachable shell index: the largest |shiftedindex()|
  // component is NL/2, so the largest possible ntnorm is (NL/2)*sqrt(3).
  // Padded by a couple of shells purely as a defensive margin.
  const int maxW = static_cast<int>(std::ceil((NL / 2.0) * std::sqrt(3.0))) + 2;

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  if (nthreads < 1) nthreads = 1;

  // ---- Pass 1 (parallel, no RNG): histogram of shell membership. ----
  std::vector<std::vector<uint64_t>> localHist(nthreads, std::vector<uint64_t>(maxW + 1, 0));
#pragma omp parallel
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    std::vector<uint64_t> &hist = localHist[tid];
#pragma omp for schedule(static)
    for (long long idx = 0; idx < static_cast<long long>(N); idx++)
    {
      long long plane = static_cast<long long>(NL) * NL;
      int i = static_cast<int>(idx / plane);
      long long rem = idx % plane;
      int j = static_cast<int>(rem / NL);
      int k = static_cast<int>(rem % NL);
      int w = computeShell(i, j, k);
      if (w >= 1 && w <= maxW) hist[w]++;
    }
  }
  std::vector<uint64_t> shellCount(maxW + 1, 0);
  for (int t = 0; t < nthreads; t++)
    for (int w = 1; w <= maxW; w++)
      shellCount[w] += localHist[t][w];
  localHist.clear();
  localHist.shrink_to_fit();

  // Prefix sums: start[w] is where shell w's block begins in `order`.
  std::vector<uint64_t> start(maxW + 2, 0);
  for (int w = 1; w <= maxW; w++) start[w + 1] = start[w] + shellCount[w];
  uint64_t totalValid = start[maxW + 1];

  // ---- Pass 2 (serial, no RNG): bucket-sort point indices by shell. ----
  // Runs single-threaded, visiting idx in strictly increasing order, so
  // that within each shell's block the points end up in exactly the same
  // row-major order the original triple loop would have visited them in --
  // this is what pass 3 below relies on for bit-for-bit-identical draws.
  std::vector<uint32_t> order(totalValid);
  {
    std::vector<uint64_t> cur(start.begin(), start.begin() + maxW + 1); // cur[w], w=0..maxW (cur[0] unused)
    long long plane = static_cast<long long>(NL) * NL;
    for (long long idx = 0; idx < static_cast<long long>(N); idx++)
    {
      int i = static_cast<int>(idx / plane);
      long long rem = idx % plane;
      int j = static_cast<int>(rem / NL);
      int k = static_cast<int>(rem % NL);
      int w = computeShell(i, j, k);
      if (w >= 1 && w <= maxW)
      {
        order[cur[w]] = static_cast<uint32_t>(idx);
        cur[w]++;
      }
    }
  }

  // ---- Pass 3: process shells strictly in increasing order (the shared
  // RNG engine must advance across shells in exactly this order, matching
  // the original's sequential dwk(1,...), dwk(2,...), ... calls). ----
  std::mt19937 engine(std::hash<int>{}(seed));
  const long long plane = static_cast<long long>(NL) * NL;

  for (int w = 1; w <= maxW; w++)
  {
    uint64_t cnt = shellCount[w];
    if (cnt == 0) continue;
    uint64_t base = start[w];

    // (a) serial, RNG: independent draws in canonical row-major order --
    // identical draw sequence to the original per-wavenumber dwk().
    for (uint64_t p = 0; p < cnt; p++)
    {
      uint32_t idx = order[base + p];
      int i = static_cast<int>(idx / plane);
      long long rem = idx % plane;
      int j = static_cast<int>(rem / NL);
      int k = static_cast<int>(rem % NL);
      if (realpoint(i, j, k))
        gk[idx] = dist(engine);
      else if (complexpoint(i, j, k))
        gk[idx] = (dist(engine) + II * dist(engine)) / sqrt(2);
    }

    // (b) parallel, no RNG: dependent (conjugate-reflection) points. Each
    // one's mirror image lies in the SAME shell (reflection preserves
    // |n|) and was already written in (a), since realpoint/complexpoint
    // always pick the independent representative of each mirror pair.
#pragma omp parallel for schedule(static)
    for (uint64_t p = 0; p < cnt; p++)
    {
      uint32_t idx = order[base + p];
      int i = static_cast<int>(idx / plane);
      long long rem = idx % plane;
      int j = static_cast<int>(rem / NL);
      int k = static_cast<int>(rem % NL);
      if (!(realpoint(i, j, k) || complexpoint(i, j, k)))
      {
        int ip = (i == 0) ? 0 : NL - i;
        int jp = (j == 0) ? 0 : NL - j;
        int kp = (k == 0) ? 0 : NL - k;
        gk[idx] = conj(gk[IDX(ip, jp, kp)]);
      }
    }

    // (c) parallel, no RNG: normalize this shell to unit variance and fold
    // in the lognormal weight for wavenumber w. gk's entries for this
    // shell are finalized here and never touched by any other shell (dn=1
    // shells are disjoint).
    double scale = gkWeight(w) / sqrt(static_cast<double>(cnt));
#pragma omp parallel for schedule(static)
    for (uint64_t p = 0; p < cnt; p++)
    {
      uint32_t idx = order[base + p];
      gk[idx] *= scale;
    }
  }

  return gk;
}

int shiftedindex(int n)
{
  if (n <= NL / 2)
  {
    return n;
  }
  else
  {
    return n - NL;
  }
}

bool realpoint(int nx, int ny, int nz)
{
  return (nx == 0 || nx == NL / 2) && (ny == 0 || ny == NL / 2) && (nz == 0 || nz == NL / 2);
}

bool complexpoint(int nx, int ny, int nz)
{
  int nxt = shiftedindex(nx);
  int nyt = shiftedindex(ny);
  int nzt = shiftedindex(nz);

  return (1 <= nxt && nxt != NL / 2 && nyt != NL / 2 && nzt != NL / 2) ||
         (nxt == NL / 2 && nyt != NL / 2 && 1 <= nzt && nzt != NL / 2) || (nxt != NL / 2 && 1 <= nyt && nyt != NL / 2 && nzt == NL / 2) || (1 <= nxt && nxt != NL / 2 && nyt == NL / 2 && nzt != NL / 2) ||
         (nxt == 0 && nyt != NL / 2 && 1 <= nzt && nzt != NL / 2) ||
         (nxt == NL / 2 && nyt == NL / 2 && 1 <= nzt && nzt != NL / 2) || (nxt == NL / 2 && 1 <= nyt && nyt != NL / 2 && nzt == NL / 2) || (1 <= nxt && nxt != NL / 2 && nyt == NL / 2 && nzt == NL / 2) ||
         (nxt == 0 && 1 <= nyt && nyt != NL / 1 && nzt == 0) ||
         (nxt == NL / 2 && 1 <= nyt && nyt != NL / 2 && nzt == 0) || (1 <= nxt && nxt != NL / 2 && nyt == 0 && nzt == NL / 2) || (nxt == 0 && nyt == NL / 2 && 1 <= nzt && nzt != NL / 2);
}
