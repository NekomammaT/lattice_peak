#include <iostream>
#include <fstream>
#include <random>
#include <array>
#include <algorithm>
#include <charconv>
#include <cmath>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "fft.hpp"

// ============================================================================
// Performance notes (what changed vs. the original Gaussian_mono_peak.cpp)
// ----------------------------------------------------------------------
// 1. All std::vector<std::vector<std::vector<T>>> 3D arrays were replaced by
//    flat, row-major std::vector<T> (Grid / RGrid, index via IDX(i,j,k)).
//    The nested-vector version allocates NL*NL+NL small vectors per array
//    and needs 3 pointer dereferences per element access; the flat version
//    is one allocation and one multiply-add per access, which is both much
//    faster and makes safe OpenMP parallelisation straightforward.
// 2. fft.hpp no longer copies data in/out through a triple-nested loop of
//    vector-of-vector-of-vector indexing; it memcpy's the already-flat
//    buffer into an FFTW-aligned scratch buffer once. The FFTW plan itself
//    is created ONCE (FFTW_MEASURE) and reused for every one of the ~30
//    same-size transforms this program performs, and FFTW is linked against
//    its threaded backend (fftw3_threads) so each individual transform is
//    itself parallelised across cores.
// 3. The "compaction" section used to allocate one full NL^3 array per
//    smoothing radius (up to ~25 radii) for the field itself AND for its
//    x/y/z gradients AND for the "rotated" (neighbour-shifted) copies of
//    those gradients -- about 22 GB of RAM at NL=256. Peak-finding across
//    radius only ever needs three consecutive shells' values plus the
//    gradient of the middle one, so this is now computed with a 3-shell
//    sliding window, cutting peak memory by roughly an order of magnitude
//    and removing the resulting cache/allocator thrashing. The "rotated"
//    copies were removed entirely -- they are just a neighbour lookup into
//    an array we already have, so no need to store them at all.
// 4. All embarrassingly-parallel element-wise loops (spectral scaling,
//    finite-difference gradients, peak-condition tests) are parallelised
//    with OpenMP over the (i, j) plane, keeping the innermost k loop
//    sequential and contiguous in memory for cache-friendliness. Writes to
//    the (rare) peak output files are protected with `omp critical`.
// 5. The random-number generation in dwk() is INTENTIONALLY left serial
//    (single-threaded, same draw order as the original): parallelising a
//    shared std::mt19937 stream would change which normal deviates land on
//    which grid point for a given seed, silently changing every downstream
//    result. That loop only touches a thin shell of the k-grid (a few
//    thousand points, not NL^3) so it is not a performance bottleneck
//    anyway; it has simply been restructured to visit the shell points via
//    a pre-built list instead of re-scanning all NL^3 points three times.
// 6. CSV writing of the two full-box fields (map, laplacian -- 16.7M
//    numbers each at NL=256) used ostream's `operator<<` one value at a
//    time, which is dominated by iostream formatting overhead. It is
//    replaced by std::to_chars (verified byte-for-byte identical to the
//    original `operator<<` output, i.e. still 6 significant digits,
//    general format) run in parallel per chunk, then concatenated once.
// 7. Additional memory trims for large boxes (e.g. NL=1024, where a single
//    NL^3 complex array is ~17 GB):
//    - Define NO_FULL_FIELD_OUTPUT to skip the unbiased-map field entirely
//      (it was only ever used to write mapfile) and skip building the
//      laplacianfile CSV (D2gx/D2D2gx themselves are still computed --
//      Lpeakfile still needs them -- only the giant CSV dump is skipped).
//      valuesToCSV's own temporary string buffer is itself a full extra
//      NL^3-sized allocation, so skipping it is a real memory saving, not
//      just an I/O one.
//    - D2gk/D2D2gk (the spectral fields) are now released (`.clear()` +
//      swap-with-empty) the moment their transform has been taken, instead
//      of sitting alive until the end of main() unused. Likewise D2gx,
//      D2D2gx are scoped to a block that ends right after Laplacian-peak
//      detection, since the compaction section afterwards never touches
//      them.
//    - The DxD2gx/DyD2gx/DzD2gx gradient arrays (and, in the compaction
//      peak test, Dx/Dy/Dz) are gone entirely: a finite difference and its
//      "rotated" (neighbour) value are just two subtractions of values
//      already sitting in D2gx / the compaction shell, so they are now
//      computed on the fly with a couple of extra neighbour look-ups
//      instead of being pre-computed into three more full-size arrays.
//    - In the compaction sliding window, the shell about to fall out of
//      the window is freed *before* the next shell is generated rather
//      than after, so the transient peak (while generating a new shell)
//      only ever holds 2 old shells + 1 new one, not 3 + 1.
//    See the reply that introduced these for the resulting memory budget.
//
// Build with the accompanying Makefile: -fopenmp -lfftw3_threads -lfftw3.
// Control the thread count with the OMP_NUM_THREADS environment variable
// if you run several seeds concurrently on the same machine.
// ============================================================================

Grid dwk(int wavenumber, int seed);
double WRTH(double z);
int shiftedindex(int n); // shifted index
bool innsigma(int nx, int ny, int nz, double wavenumber); // judge if point is in nsigma sphere shell
bool realpoint(int nx, int ny, int nz);                   // judge real point
bool complexpoint(int nx, int ny, int nz);                // judge independent complex point

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
const int NL = 512; //256; // Box size NL
const int nsigma = 32; //16;
const double As = 1e-2; //3.625e-3;
const double dn = 1; // Thickness of nsigma sphere shell
const std::string mapfileprefix = std::string("data/mono_map_") + std::to_string(NL) + std::string("_") + std::to_string(nsigma) + std::string("_");
const std::string laplacianfileprefix = std::string("data/mono_laplacian_") + std::to_string(NL) + std::string("_") + std::to_string(nsigma) + std::string("_");
const std::string Lpeakfileprefix = std::string("data/mono_Lpeak_") + std::to_string(NL) + std::string("_") + std::to_string(nsigma) + std::string("_");
const std::string Cpeakfileprefix = std::string("data/mono_Cpeak_") + std::to_string(NL) + std::string("_") + std::to_string(nsigma) + std::string("_");

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
// the un-smoothed Fourier-space map gk. This is exactly the per-radius body
// of the original compaction loop, factored out so it can be called from a
// 3-shell sliding window instead of building the full history up front.
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

// Peak test across the radius direction (r-1, r, r+1 in the *original*
// r-loop's indexing, here passed in as C0/C1/C2) combined with the spatial
// saddle-point test on C0's gradient. This is exactly the body of the
// original `for (r ...)` peak-detection loop, just reading its inputs from
// the sliding window instead of fully materialised per-radius arrays.
void processCompactionPeaks(int r, const RGrid &C0, const RGrid &C1, const RGrid &C2, std::ofstream &Cpeakfile)
{
  // Each thread appends its hits, in the order it encounters them, to its
  // own buffer; buffers are then written out in thread-index order. Since
  // OpenMP's static schedule hands out contiguous, increasing blocks of the
  // collapsed (i,j) space to increasing thread numbers, this reproduces
  // exactly the same line order the original fully-serial i/j/k loop would
  // have produced, while still doing the scan itself in parallel.
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
            // Dx/Dy/Dz (the finite-difference gradient of C0) and their
            // "rotated" (one-neighbour-over) values are just subtractions
            // of C0 at up to 3 points along each axis (i, i+1, i+2 with
            // wraparound) -- computed here instead of in three
            // full-size pre-computed arrays, saving 3 * NL^3 doubles.
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

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cerr << "Specify the noise file number correctly." << std::endl;
    return 1;
  }

  // ---------- start timer ----------
  struct timeval Nv;
  struct timezone Nz;
  double before, after;

  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  // --------------------------------------

  int seed = atoi(argv[1]);
#ifndef NO_FULL_FIELD_OUTPUT
  std::ofstream mapfile(mapfileprefix + std::to_string(seed) + ".csv");
  std::ofstream laplacianfile(laplacianfileprefix + std::to_string(seed) + ".csv");
#endif
  std::ofstream Lpeakfile(Lpeakfileprefix + std::to_string(seed) + ".csv");
  std::ofstream Cpeakfile(Cpeakfileprefix + std::to_string(seed) + ".csv");

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  FFTW3DPlan fftplan(NL, nthreads);

  // ----------- unbiased map -----------
  // gx (the real-space unbiased map) is only ever used to write mapfile, so
  // with NO_FULL_FIELD_OUTPUT defined we skip computing it entirely -- one
  // fewer NL^3 complex array (~17 GB at NL=1024) and one fewer transform.
  Grid gk = dwk(nsigma, seed);
#ifndef NO_FULL_FIELD_OUTPUT
  {
    Grid gx = fftplan.execute(gk);
    mapfile << valuesToCSV(gx.size(), [&](size_t idx) { return gx[idx].real(); }) << std::endl;
  }
  std::cout << "Exported to " << mapfileprefix + std::to_string(seed) + ".csv" << std::endl;
#endif

  // ----------- laplacian + gradient + Laplacian-peak detection -----------
  // D2gx/D2D2gx (and, before this rewrite, DxD2gx/DyD2gx/DzD2gx) are scoped
  // to this block: nothing outside it -- the compaction section below --
  // ever touches them, so their memory (2-3 more NL^3 arrays) is released
  // as soon as the block ends, instead of sitting alive, unused, until the
  // program exits.
  {
    Grid D2gk(gk.size()), D2D2gk(gk.size());
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
          D2gk[idx] = gk[idx] * pow(2 * M_PI * ntnorm / NL, 2);
          D2D2gk[idx] = gk[idx] * pow(2 * M_PI * ntnorm / NL, 4);
        }

    // Take the two transforms one at a time and drop each D2*gk the moment
    // its transform is in hand, rather than keeping both spectral copies
    // around until after both transforms are done -- one fewer NL^3
    // complex array alive at the peak instant.
    Grid D2gx = fftplan.execute(D2gk);
    D2gk = Grid(); // free
    Grid D2D2gx = fftplan.execute(D2D2gk);
    D2D2gk = Grid(); // free

#ifndef NO_FULL_FIELD_OUTPUT
    laplacianfile << valuesToCSV(D2gx.size(), [&](size_t idx) { return D2gx[idx].real(); }) << std::endl;
    std::cout << "Exported to " << laplacianfileprefix + std::to_string(seed) + ".csv" << std::endl;
#endif

    // Laplacian peak detection. The finite-difference gradient of D2gx and
    // its "rotated" (one-neighbour-over) value are just subtractions of
    // D2gx at up to 3 points along each axis (i, i+1, i+2 with wraparound)
    // -- computed here on the fly instead of in three full-size
    // pre-computed DxD2gx/DyD2gx/DzD2gx arrays, saving 3 * NL^3 doubles
    // (~26 GB at NL=1024). As with the compaction peaks below, each thread
    // buffers its hits in encounter order and buffers are concatenated in
    // thread order afterwards, so the file ends up in the same row order
    // the serial original produced.
    int nthreadsL = 1;
#ifdef _OPENMP
    nthreadsL = omp_get_max_threads();
#endif
    std::vector<std::string> lbuffers(nthreadsL);
#pragma omp parallel
    {
      int tid = 0;
#ifdef _OPENMP
      tid = omp_get_thread_num();
#endif
      std::string &buf = lbuffers[tid];
#pragma omp for collapse(2) schedule(static)
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
            if (dx * dxrot < 0 && dx > 0 &&
                dy * dyrot < 0 && dy > 0 &&
                dz * dzrot < 0 && dz > 0)
            {
              size_t idxp = IDX(ip, jp, kp);
              appendNum(buf, ip); buf += ',';
              appendNum(buf, jp); buf += ',';
              appendNum(buf, kp); buf += ',';
              appendNum(buf, D2gx[idxp].real()); buf += ',';
              appendNum(buf, sqrt(D2D2gx[idxp].real() / D2gx[idxp].real())); buf += '\n';
            }
          }
    }
    for (auto &buf : lbuffers) Lpeakfile << buf;
    // D2gx, D2D2gx (and lbuffers) go out of scope here and are freed before
    // the compaction section below starts.
  }

  std::cout << "Exported to " << Lpeakfileprefix + std::to_string(seed) + ".csv" << std::endl;

  // ----------- compaction (3-shell sliding window over radius) ----------
  double Rbound = 5. / (2 * M_PI * nsigma / NL); //10. / (2 * M_PI * nsigma / NL);
  RGrid C0, C1, C2;
  int haveCount = 0;

  for (int rs = 1; rs <= Rbound; rs++)
  {
    if (haveCount == 3)
    {
      // C0 is about to fall out of the window on this iteration's slide
      // (it was already consumed by the previous peak check and won't be
      // read again) -- free it *before* generating the next shell rather
      // than after, so the transient while generating a new shell only
      // ever holds 2 old shells + 1 new one, not 3 + 1.
      C0 = RGrid();
    }
    RGrid Cnew = computeCompactionShell(gk, fftplan, rs);
    int shellIdx = rs - 1; // 0-based, matches original compaction[rs-1]

    if (haveCount < 3)
    {
      if (haveCount == 0) C0 = std::move(Cnew);
      else if (haveCount == 1) C1 = std::move(Cnew);
      else C2 = std::move(Cnew);
      haveCount++;
    }
    else
    {
      C0 = std::move(C1);
      C1 = std::move(C2);
      C2 = std::move(Cnew);
    }

    if (haveCount == 3)
    {
      int r = shellIdx - 2; // matches the original loop variable r (0-based)
      processCompactionPeaks(r, C0, C1, C2, Cpeakfile);
    }
  }

  std::cout << "Exported to " << Cpeakfileprefix + std::to_string(seed) + ".csv" << std::endl;

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------

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
