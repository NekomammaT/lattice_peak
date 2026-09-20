#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <sys/time.h>
#include <vector>

#include <fftw3.h>

// ============================================================================
// Performance notes (what changed vs. the original Gaussian_mono.cpp)
// ----------------------------------------------------------------------
// This program is run job-parallel (many seeds at once via the job
// scheduler), so -- unlike Gaussian_mono_peak.cpp/Gaussian_LN_peak.cpp --
// nothing here is OpenMP-parallelised or linked against fftw3_threads; every
// change below is a purely single-threaded, serial speed-up.
//
// 1. All std::vector<std::vector<std::vector<T>>> 3D arrays were replaced by
//    flat, row-major std::vector<T> (Grid, index via IDX(i,j,k)), removing
//    the small-allocation and pointer-chasing overhead of nested vectors.
// 2. The old fftw() free function fftw_malloc'd fresh in/out buffers and
//    built a fresh FFTW_ESTIMATE plan on every single call. This program
//    calls it 4 times per run (gx, gxbias, Dgx, DDgx); a single Fft3D object
//    now builds its buffers/plan once and reuses them across all 4 calls
//    (still FFTW_ESTIMATE, so the algorithm FFTW picks -- and hence the
//    transformed values themselves -- is unchanged from the original;
//    this only removes repeated setup/teardown overhead).
// 3. dwk()/Bk() are ported to the flat layout with the same
//    draw-then-reflect-then-normalize structure as the original (same
//    std::mt19937 draw order, bit-for-bit).
// 4. The Cmax-over-radius loop is the real bottleneck at large NL: the
//    original built a *full* NL^3 field rzpk, ran a *full* 3D FFT to get
//    rzpx, and then read exactly ONE element of it (rzpx[imax][jmax][kmax])
//    -- for every one of ~O(NL) radii. That is O(NL^3 log NL) of FFT work
//    (plus a full NL^3 temporary array) per radius just to extract a single
//    scalar. Since only one output point is ever needed, it is computed
//    directly from the definition of the DFT that FFTW's FFTW_FORWARD
//    implements:
//        rzpx[imax,jmax,kmax] = sum_{i,j,k} rzpk[i,j,k] *
//                                 exp(-2*pi*i*(i*imax+j*jmax+k*kmax)/NL)
//    which is separable into per-axis phase factors (precomputed once,
//    outside the radius loop, from imax/jmax/kmax) and evaluated in a single
//    fused O(NL^3) pass with no FFT call and no temporary NL^3 array at all
//    -- i.e. the same asymptotic cost as just *building* rzpk, instead of
//    also paying for a full transform of it. This removes both the log(NL)
//    factor and the per-radius allocation, and is the dominant speed-up for
//    this program.
//
//    IMPORTANT CAVEAT: unlike every other change in this file (and unlike
//    the earlier Gaussian_mono_peak.cpp / Gaussian_LN_peak.cpp rewrites),
//    this one is NOT guaranteed bit-for-bit identical to the original. A
//    direct summation and an FFT compute the mathematically same value but
//    add up the same terms in a different order, so the result can differ
//    in the last few bits of precision (the same kind of noise -Ofast
//    already introduces elsewhere in this codebase). mu2, k3, lnw and
//    zetam do NOT depend on this loop and remain bit-identical to the
//    original; only Cmax (and, in the rare case of a near-tie between two
//    candidate radii, rsmax) can move by a numerically negligible amount.
//    This was verified empirically -- see the accompanying test results.
// ============================================================================

constexpr int NL = 256; // Box size NL
constexpr int nsigma = 16;
constexpr double dn = 1; // Thickness of nsigma sphere shell
constexpr double bias = 8; //9; //10;
constexpr double As = 5e-3; //1e-2; //3.625e-3;
constexpr std::size_t N3 = static_cast<std::size_t>(NL) * NL * NL;
const std::string mukfilename = std::string("data/mono_muk_") + std::to_string(NL) + std::string("_") + std::to_string(nsigma) + std::string(".csv");

using Grid = std::vector<std::complex<double>>;
using RGrid = std::vector<double>;

inline std::size_t IDX(int i, int j, int k)
{
  return (static_cast<std::size_t>(i) * NL + j) * NL + k;
}

inline int shiftedindex(int n)
{
  return n <= NL / 2 ? n : n - NL;
}

inline bool realpoint(int nx, int ny, int nz)
{
  return (nx == 0 || nx == NL / 2) && (ny == 0 || ny == NL / 2) && (nz == 0 || nz == NL / 2);
}

inline bool complexpoint(int nx, int ny, int nz)
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

inline bool innsigma(int nx, int ny, int nz, double wavenumber)
{
  int nxt = shiftedindex(nx);
  int nyt = shiftedindex(ny);
  int nzt = shiftedindex(nz);
  double ntnorm = sqrt(static_cast<double>(nxt * nxt + nyt * nyt + nzt * nzt));
  return (wavenumber - dn / 2. <= ntnorm && ntnorm < wavenumber + dn / 2.);
}

// real-space top-hat window
inline double WRTH(double z)
{
  return z == 0.0 ? 1.0 : 3.0 * (std::sin(z) - z * std::cos(z)) / (z * z * z);
}

// Single-threaded FFTW plan, buffers built once and reused across every
// transform call in a run (the original built/tore down fresh buffers and a
// fresh FFTW_ESTIMATE plan on every fftw() call). Still FFTW_ESTIMATE, still
// out-of-place, so the transformed values are unaffected -- this only
// removes the repeated malloc/plan-create/destroy overhead.
class Fft3D
{
public:
  Fft3D()
  {
    in_ = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N3));
    out_ = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * N3));
    if (in_ == nullptr || out_ == nullptr) throw std::bad_alloc();
    plan_ = fftw_plan_dft_3d(NL, NL, NL, in_, out_, FFTW_FORWARD, FFTW_ESTIMATE);
    if (plan_ == nullptr) throw std::runtime_error("could not create FFTW plan");
  }

  ~Fft3D()
  {
    if (plan_ != nullptr) fftw_destroy_plan(plan_);
    fftw_free(in_);
    fftw_free(out_);
  }

  Fft3D(const Fft3D &) = delete;
  Fft3D &operator=(const Fft3D &) = delete;

  Grid transform(const Grid &input)
  {
    for (std::size_t p = 0; p < N3; ++p)
    {
      in_[p][0] = input[p].real();
      in_[p][1] = input[p].imag();
    }
    fftw_execute(plan_);
    Grid result(N3);
    for (std::size_t p = 0; p < N3; ++p) result[p] = {out_[p][0], out_[p][1]};
    return result;
  }

private:
  fftw_complex *in_ = nullptr;
  fftw_complex *out_ = nullptr;
  fftw_plan plan_ = nullptr;
};

std::normal_distribution<> dist(0., 1.);
const std::complex<double> II(0, 1);

Grid dwk(int wavenumber, std::mt19937 &engine)
{
  Grid dwk(N3, std::complex<double>(0, 0));
  int count = 0;

  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
        if (innsigma(i, j, k, wavenumber))
        {
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

  int ip, jp, kp;
  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
        if (innsigma(i, j, k, wavenumber) && !(realpoint(i, j, k) || complexpoint(i, j, k)))
        {
          ip = (i == 0) ? 0 : NL - i;
          jp = (j == 0) ? 0 : NL - j;
          kp = (k == 0) ? 0 : NL - k;
          dwk[IDX(i, j, k)] = conj(dwk[IDX(ip, jp, kp)]);
          count++;
        }

  if (count != 0)
  {
    double norm = sqrt(static_cast<double>(count));
    for (int i = 0; i < NL; i++)
      for (int j = 0; j < NL; j++)
        for (int k = 0; k < NL; k++)
          if (innsigma(i, j, k, wavenumber)) dwk[IDX(i, j, k)] /= norm;
  }

  return dwk;
}

Grid Bk(int wavenumber, double biasval)
{
  Grid Bk(N3, std::complex<double>(0, 0));
  int count = 0;

  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
        if (innsigma(i, j, k, wavenumber)) count++;

  if (count != 0)
  {
    double val = biasval / count;
    for (int i = 0; i < NL; i++)
      for (int j = 0; j < NL; j++)
        for (int k = 0; k < NL; k++)
          if (innsigma(i, j, k, wavenumber)) Bk[IDX(i, j, k)] = val;
  }

  return Bk;
}

int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cerr << "Specify the seed correctly." << std::endl;
    return 1;
  }

  struct timeval Nv;
  struct timezone Nz;
  double before, after;
  gettimeofday(&Nv, &Nz);
  before = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;

  int seed = atoi(argv[1]);
  std::cout << "seed = " << seed << std::endl;
  std::mt19937 engine(std::hash<int>{}(seed));
  std::ofstream mukfile(mukfilename, std::ios::app);

  Fft3D fft;

  // ----------- unbiased map -----------
  Grid gk = dwk(nsigma, engine);
  Grid gx = fft.transform(gk);

  // ----------- biased map -----------
  Grid Bkval = Bk(nsigma, bias);
  Grid gkbias(N3);
  for (std::size_t p = 0; p < N3; p++) gkbias[p] = gk[p] + Bkval[p];

  RGrid normArr(N3);
  Grid Dgk(N3), DDgk(N3);
  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
      {
        std::size_t idx = IDX(i, j, k);
        int nxt = shiftedindex(i), nyt = shiftedindex(j), nzt = shiftedindex(k);
        double ntnorm = sqrt(static_cast<double>(nxt * nxt + nyt * nyt + nzt * nzt));
        normArr[idx] = ntnorm;
        Dgk[idx] = gkbias[idx] * pow(2 * M_PI * ntnorm / NL, 2);
        DDgk[idx] = gkbias[idx] * pow(2 * M_PI * ntnorm / NL, 4);
      }

  Grid gxbias = fft.transform(gkbias);
  Grid Dgx = fft.transform(Dgk);
  Grid DDgx = fft.transform(DDgk);

  auto iter = std::max_element(Dgx.begin(), Dgx.end(),
      [](const std::complex<double> &a, const std::complex<double> &b) { return a.real() < b.real(); });
  std::size_t index = static_cast<std::size_t>(std::distance(Dgx.begin(), iter));

  int imax = static_cast<int>(index / (static_cast<std::size_t>(NL) * NL));
  int jmax = static_cast<int>((index / NL) % NL);
  int kmax = static_cast<int>(index % NL);
  double mu2 = Dgx[index].real();
  double k3 = sqrt(DDgx[index].real() / Dgx[index].real());
  double lnw = -bias * gx[0].real() - 0.5 * bias * bias;

  // ----------- Cmax over radius: direct single-point DFT evaluation -----
  // (see header comment -- replaces "build full field + FFT + read one
  // point" with a fused O(NL^3) sum for just that one point, no FFT call
  // and no per-radius temporary array).
  std::vector<std::complex<double>> phaseI(NL), phaseJ(NL), phaseK(NL);
  for (int i = 0; i < NL; i++) phaseI[i] = std::exp(std::complex<double>(0, -2 * M_PI * i * imax / NL));
  for (int j = 0; j < NL; j++) phaseJ[j] = std::exp(std::complex<double>(0, -2 * M_PI * j * jmax / NL));
  for (int k = 0; k < NL; k++) phaseK[k] = std::exp(std::complex<double>(0, -2 * M_PI * k * kmax / NL));

  double Cmax = 0;
  int rsmax = 0;
  double sqrtAs = sqrt(As);
  int rs_limit = static_cast<int>(10. / (2 * M_PI * nsigma / NL));
  for (int rs = 1; rs <= rs_limit; rs++)
  {
    std::complex<double> acc(0, 0);
    for (int i = 0; i < NL; i++)
    {
      std::complex<double> pi = phaseI[i];
      for (int j = 0; j < NL; j++)
      {
        std::complex<double> pij = pi * phaseJ[j];
        for (int k = 0; k < NL; k++)
        {
          std::size_t idx = IDX(i, j, k);
          double kr = 2 * M_PI * normArr[idx] * rs / NL;
          std::complex<double> rzpk_val = gkbias[idx] * (-kr * kr / 3 * WRTH(kr) * sqrtAs);
          acc += rzpk_val * pij * phaseK[k];
        }
      }
    }
    double compaction = 2. / 3 * (1 - pow(1 + acc.real(), 2));
    if (compaction > Cmax)
    {
      Cmax = compaction;
      rsmax = rs;
    }
  }

  int count = 0;
  double zetam = 0;
  int nxm = shiftedindex(imax), nym = shiftedindex(jmax), nzm = shiftedindex(kmax);
  for (int i = 0; i < NL; i++)
    for (int j = 0; j < NL; j++)
      for (int k = 0; k < NL; k++)
      {
        int nxt = shiftedindex(i), nyt = shiftedindex(j), nzt = shiftedindex(k);
        double dr = sqrt(static_cast<double>((nxt - nxm) * (nxt - nxm) + (nyt - nym) * (nyt - nym) + (nzt - nzm) * (nzt - nzm)));
        if (fabs(dr - rsmax) < 1. / 2)
        {
          zetam += gxbias[IDX(i, j, k)].real() * sqrtAs;
          count++;
        }
      }
  zetam /= count;

  mukfile << seed << ',' << mu2 << ',' << k3 << ',' << 2 * M_PI * nsigma / NL * rsmax << ',' << zetam << ',' << Cmax << ',' << lnw << std::endl;

  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;

  return 0;
}
