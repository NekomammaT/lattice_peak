#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <charconv>
#include <array>
#include <complex>
#include <cmath>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "fft.hpp"

// ============================================================================
// Monochromatic (Dirac-delta power spectrum, single Fourier shell k = nsigma)
// version of Gaussian_LN_peakprof.cpp: per-Laplacian-peak compaction profiles.
// The Fourier map gk, its shell construction dwk(), the Laplacian-peak test and
// the compaction/zeta shell computations are taken from Gaussian_mono_peak.cpp /
// Gaussian_LN_peakprof.cpp unchanged. Differences from the LN version:
//   * gk = dwk(nsigma, seed): one unit-normalised shell (sum |gk|^2 = 1).
//   * sigma2 (for nu2 = mu2/sigma2) is evaluated exactly from the realised
//     shell, sigma2 = sqrt(sum_k |gk|^2 k^4) (-> kpk^2 for a thin shell).
// Output: data/mono_Cprof_<NL>_<nsigma>_<seed>.csv, same columns as LN_Cprof.
// ============================================================================

Grid dwk(int wavenumber, int seed);
double WRTH(double z);
bool innsigma(int nx, int ny, int nz, double wavenumber);
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
const double As = 5e-3;
const double dn = 1; // Thickness of nsigma sphere shell
const std::string profprefix = std::string("data/mono_Cprof_") + std::to_string(NL) + "_" + std::to_string(nsigma) + "_";

// flat-index helper: matches the (i,j,k) -> i*NL*NL + j*NL + k layout used
// throughout (row-major, same iteration order as the original triple loop)
inline size_t IDX(int i, int j, int k)
{
  return (static_cast<size_t>(i) * NL + j) * NL + k;
}

inline int nextIndex(int n) { return (n == NL - 1) ? 0 : n + 1; }

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
// Output: data/mono_Cprof_<NL>_<nsigma>_<seed>.csv
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
  std::ofstream proffile(profprefix + std::to_string(seed) + ".csv");

  Grid gk = dwk(nsigma, seed);
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  FFTW3DPlan fftplan(NL, nthreads);

  // nu2 = mu2 / sigma2, sigma2 = sqrt(sum_k |gk|^2 k^4) from the realised shell (variance of Laplacian)
  double sigma2th;
  {
    double acc = 0;
#pragma omp parallel for reduction(+ : acc) collapse(2) schedule(static)
    for (int i = 0; i < NL; i++)
      for (int j = 0; j < NL; j++)
        for (int k = 0; k < NL; k++)
        {
          int nxt = shiftedindex(i), nyt = shiftedindex(j), nzt = shiftedindex(k);
          double k2 = pow(2 * M_PI / NL, 2) * (nxt * nxt + nyt * nyt + nzt * nzt);
          acc += std::norm(gk[IDX(i, j, k)]) * k2 * k2;
        }
    sigma2th = sqrt(acc);
  }
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

Grid dwk(int wavenumber, int seed)
{
  Grid dwk(static_cast<size_t>(NL) * NL * NL, std::complex<double>(0, 0));

  // Build the shell point list once (single NL^3 sweep, no RNG draws),
  // instead of re-scanning all NL^3 points on every one of the three
  // passes below. Iteration order matches the original triple loop
  // exactly (i outer, k inner), so the sequence of dist(engine) draws for
  // a given seed is unchanged -- this only removes redundant innsigma()
  // sweeps, never reorders random draws.
  std::vector<std::array<int, 3>> shellPoints;
#pragma omp parallel
  {
    std::vector<std::array<int, 3>> local;
#pragma omp for collapse(2) schedule(static) nowait
    for (int i = 0; i < NL; i++)
      for (int j = 0; j < NL; j++)
        for (int k = 0; k < NL; k++)
          if (innsigma(i, j, k, wavenumber))
            local.push_back({i, j, k});
#pragma omp critical
    shellPoints.insert(shellPoints.end(), local.begin(), local.end());
  }
  // Threads may finish their (i,j) chunks in any order, so re-sort into the
  // canonical row-major order to guarantee the RNG draw sequence below is
  // bit-for-bit identical to the fully serial original.
  std::sort(shellPoints.begin(), shellPoints.end());

  int count = 0;
  std::mt19937 engine(std::hash<int>{}(seed));

  // Pass 1: draw independent degrees of freedom (real points and
  // independent complex points), in exactly the original i,j,k order.
  // Kept single-threaded: this consumes the seeded RNG stream, so its
  // draw order must not change.
  for (auto &p : shellPoints)
  {
    int i = p[0], j = p[1], k = p[2];
    if (realpoint(i, j, k))
    {
      dwk[IDX(i, j, k)] = dist(engine);
      count++;
    }
    else if (complexpoint(i, j, k))
    {
      dwk[IDX(i, j, k)] = (dist(engine) + II * dist(engine)) / sqrt(2);
      count++;
    }
  }

  // Pass 2: reflection (no RNG draws -- safe to parallelise, each (i,j,k)
  // is written exactly once and only reads already-finalised source
  // points from pass 1).
  int ip, jp, kp;
#pragma omp parallel for private(ip, jp, kp) schedule(dynamic, 64)
  for (size_t idx = 0; idx < shellPoints.size(); idx++)
  {
    int i = shellPoints[idx][0], j = shellPoints[idx][1], k = shellPoints[idx][2];
    if (!(realpoint(i, j, k) || complexpoint(i, j, k)))
    {
      ip = (i == 0) ? 0 : NL - i;
      jp = (j == 0) ? 0 : NL - j;
      kp = (k == 0) ? 0 : NL - k;
      dwk[IDX(i, j, k)] = conj(dwk[IDX(ip, jp, kp)]);
#pragma omp atomic
      count++;
    }
  }

  if (count != 0)
  {
    double norm = sqrt(static_cast<double>(count));
#pragma omp parallel for schedule(static)
    for (size_t idx = 0; idx < shellPoints.size(); idx++)
    {
      int i = shellPoints[idx][0], j = shellPoints[idx][1], k = shellPoints[idx][2];
      dwk[IDX(i, j, k)] /= norm;
    }
  }

  return dwk;
}

double WRTH(double z);
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

// Fourier points
bool innsigma(int nx, int ny, int nz, double wavenumber)
{
  int nxt = shiftedindex(nx);
  int nyt = shiftedindex(ny);
  int nzt = shiftedindex(nz);

  double ntnorm = sqrt(nxt * nxt + nyt * nyt + nzt * nzt);

  return (wavenumber - dn / 2. <= ntnorm && ntnorm < wavenumber + dn / 2.);
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
