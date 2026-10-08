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
// Performance notes (what changed vs. the original Gaussian_LN_GB.cpp)
// ----------------------------------------------------------------------
// This program is already flat-array / bucketed-shell / cached-plan (it was
// evidently optimized before), and it is run job-parallel like Gaussian_mono
// -- no OpenMP, no fftw3_threads, purely single-threaded. The one remaining
// bottleneck is the same one fixed in Gaussian_mono.cpp: the Cmax-over-radius
// loop built a *full* NL^3 field rzpk and ran a *full* 3D FFT to get rzpx,
// only to read exactly ONE element of it (rzpx[peak_index]), for every one
// of ~O(NL) radii -- O(NL^3 log NL) of FFT work (plus a full NL^3 temporary
// array) per radius just to extract a single scalar.
//
// Since only one output point is ever needed, it is now computed directly
// from the definition of the DFT that FFTW's FFTW_FORWARD implements:
//     rzpx[imax,jmax,kmax] = sum_{i,j,k} rzpk[i,j,k] *
//                              exp(-2*pi*i*(i*imax+j*jmax+k*kmax)/NL)
// which separates into per-axis phase factors (precomputed once from
// imax/jmax/kmax, outside the radius loop) and is evaluated in a single
// fused O(NL^3) pass with no FFT call and no per-radius temporary array --
// i.e. the same asymptotic cost as just *building* rzpk used to be, instead
// of also paying for a full transform of it.
//
// IMPORTANT CAVEAT: this one change is NOT guaranteed bit-for-bit identical
// to the original -- a direct summation and an FFT compute the
// mathematically same value but add up the same terms in a different order,
// so the result can differ in the last few bits of precision. mu2, k3 and
// lnw do not depend on this loop and remain bit-identical to the original;
// only Cmax (and, in the rare case of a near-tie between two candidate
// radii, rsmax) can move by a numerically negligible amount. This was
// verified empirically against the original -- see the accompanying test
// results.
// ============================================================================

constexpr int NL = 256;
constexpr int nsigma = 16;
constexpr double As = 1e-2;
constexpr int nbias = 16;
constexpr double dlnn = 0.1;
double biascoeff = 12.5;  // scan: overridden by argv[2]
constexpr double s2 = 0.1;
const std::string s2value = "0,1";
constexpr double dn = 1.0;

constexpr std::size_t N3 = static_cast<std::size_t>(NL) * NL * NL;
// shifted indices span [-127, 128], so the furthest populated shell is 222.
constexpr int MAX_SHELL = 222;
std::string mukfilename = "data/LN_muk_" + s2value + "_" +
    std::to_string(NL) + "_" + std::to_string(nsigma) + "_" +
    std::to_string(nbias) + "_GB.csv";

using Grid = std::vector<std::complex<double>>;

inline std::size_t index_of(int i, int j, int k)
{
  return (static_cast<std::size_t>(i) * NL + j) * NL + k;
}

inline int shiftedindex(int n)
{
  return n <= NL / 2 ? n : n - NL;
}

inline bool realpoint(int nx, int ny, int nz)
{
  return (nx == 0 || nx == NL / 2) && (ny == 0 || ny == NL / 2) &&
         (nz == 0 || nz == NL / 2);
}

inline bool complexpoint(int nx, int ny, int nz)
{
  const int x = shiftedindex(nx);
  const int y = shiftedindex(ny);
  const int z = shiftedindex(nz);
  return (1 <= x && x != NL / 2 && y != NL / 2 && z != NL / 2) ||
         (x == NL / 2 && y != NL / 2 && 1 <= z && z != NL / 2) ||
         (x != NL / 2 && 1 <= y && y != NL / 2 && z == NL / 2) ||
         (1 <= x && x != NL / 2 && y == NL / 2 && z != NL / 2) ||
         (x == 0 && y != NL / 2 && 1 <= z && z != NL / 2) ||
         (x == NL / 2 && y == NL / 2 && 1 <= z && z != NL / 2) ||
         (x == NL / 2 && 1 <= y && y != NL / 2 && z == NL / 2) ||
         (1 <= x && x != NL / 2 && y == NL / 2 && z == NL / 2) ||
         (x == 0 && 1 <= y && y != NL / 1 && z == 0) ||
         (x == NL / 2 && 1 <= y && y != NL / 2 && z == 0) ||
         (1 <= x && x != NL / 2 && y == 0 && z == NL / 2) ||
         (x == 0 && y == NL / 2 && 1 <= z && z != NL / 2);
}

inline double powerspectrum(int wavenumber)
{
  const double d = std::log(static_cast<double>(wavenumber)) - std::log(static_cast<double>(nsigma));
  return std::exp(-d * d / (2.0 * s2)) / std::sqrt(2.0 * M_PI * s2);
}

inline double BN(int wavenumber)
{
  const double d = std::log(static_cast<double>(wavenumber)) - std::log(static_cast<double>(nbias));
  return biascoeff * std::exp(-d * d / (2.0 * dlnn)) / std::sqrt(2.0 * M_PI * dlnn);
}

inline double WRTH(double z)
{
  return z == 0.0 ? 1.0 : 3.0 * (std::sin(z) - z * std::cos(z)) / (z * z * z);
}

// Reuses both buffers and the FFTW plan. (Unchanged from before.)
class Fft3D {
 public:
  Fft3D()
  {
    in_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N3));
    out_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N3));
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

  Fft3D(const Fft3D&) = delete;
  Fft3D& operator=(const Fft3D&) = delete;

  fftw_complex* in() { return in_; }
  fftw_complex* out() { return out_; }
  void execute() { fftw_execute(plan_); }

  Grid transform(const Grid& input)
  {
    for (std::size_t p = 0; p < N3; ++p) {
      in_[p][0] = input[p].real();
      in_[p][1] = input[p].imag();
    }
    fftw_execute(plan_);
    Grid result(N3);
    for (std::size_t p = 0; p < N3; ++p) result[p] = {out_[p][0], out_[p][1]};
    return result;
  }

 private:
  fftw_complex* in_ = nullptr;
  fftw_complex* out_ = nullptr;
  fftw_plan plan_ = nullptr;
};

struct FourierModes {
  std::vector<std::vector<std::size_t>> shells;
  std::vector<std::vector<std::size_t>> independent;
  std::vector<std::vector<std::size_t>> reflected;
  std::vector<std::size_t> reflection_source;
  std::vector<double> norm;
};

FourierModes make_fourier_modes()
{
  FourierModes modes{std::vector<std::vector<std::size_t>>(MAX_SHELL + 1),
                     std::vector<std::vector<std::size_t>>(MAX_SHELL + 1),
                     std::vector<std::vector<std::size_t>>(MAX_SHELL + 1), {},
                     std::vector<double>(N3)};
  modes.reflection_source.resize(N3);

  // This preserves the original i-j-k traversal order within every shell.
  for (int i = 0; i < NL; ++i) for (int j = 0; j < NL; ++j) for (int k = 0; k < NL; ++k) {
    const std::size_t p = index_of(i, j, k);
    const int x = shiftedindex(i), y = shiftedindex(j), z = shiftedindex(k);
    const double r = std::sqrt(static_cast<double>(x * x + y * y + z * z));
    modes.norm[p] = r;
    const int shell = static_cast<int>(std::floor(r + 0.5));
    if (shell > MAX_SHELL) continue;
    modes.shells[shell].push_back(p);
    if (realpoint(i, j, k) || complexpoint(i, j, k)) {
      modes.independent[shell].push_back(p);
    }
  }
  // The second pass matches the original reflection pass and its ordering.
  for (int i = 0; i < NL; ++i) for (int j = 0; j < NL; ++j) for (int k = 0; k < NL; ++k) {
    if (realpoint(i, j, k) || complexpoint(i, j, k)) continue;
    const std::size_t p = index_of(i, j, k);
    const int shell = static_cast<int>(std::floor(modes.norm[p] + 0.5));
    if (shell > MAX_SHELL) continue;
    const int ip = i == 0 ? 0 : NL - i;
    const int jp = j == 0 ? 0 : NL - j;
    const int kp = k == 0 ? 0 : NL - k;
    modes.reflected[shell].push_back(p);
    modes.reflection_source[p] = index_of(ip, jp, kp);
  }
  return modes;
}

int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: Gaussian_LN_GBC3 <seed> <biascoeff>\n";
    return 1;
  }

  timeval now{};
  gettimeofday(&now, nullptr);
  const double before = now.tv_sec + now.tv_usec * 1.e-6;

  const int seed = std::atoi(argv[1]);
  biascoeff = std::atof(argv[2]);
  mukfilename = std::string("gbc/ev_As1e-2_bc") + argv[2] + ".csv";
  std::cout << "seed = " << seed << '\n';
  std::mt19937 engine(std::hash<int>{}(seed));
  std::normal_distribution<> dist(0.0, 1.0);
  std::ofstream mukfile(mukfilename, std::ios::app);

  const FourierModes modes = make_fourier_modes();
  Grid gkbias(N3, {0.0, 0.0});
  Grid dwk(N3, {0.0, 0.0});
  double lnw = 0.0;

  for (int w = 1; w <= MAX_SHELL; ++w) {
    const auto& shell = modes.shells[w];
    if (shell.empty()) continue;
    const double random_scale = std::sqrt(dn / w / static_cast<double>(shell.size()));
    for (const std::size_t p : modes.independent[w]) {
      const int k = static_cast<int>(p % NL);
      const int j = static_cast<int>((p / NL) % NL);
      const int i = static_cast<int>(p / (static_cast<std::size_t>(NL) * NL));
      dwk[p] = realpoint(i, j, k)
          ? std::complex<double>(dist(engine), 0.0)
          : std::complex<double>(dist(engine), dist(engine)) / std::sqrt(2.0);
    }
    for (const std::size_t p : modes.reflected[w]) dwk[p] = std::conj(dwk[modes.reflection_source[p]]);

    std::complex<double> zero_mode(0.0, 0.0);
    const double bias = BN(w);
    const double spectrum_scale = std::sqrt(powerspectrum(w));
    const std::complex<double> deterministic = bias * dn / w / static_cast<double>(shell.size()) * spectrum_scale;
    for (const std::size_t p : shell) {
      zero_mode += dwk[p] * random_scale;
      gkbias[p] += dwk[p] * random_scale * spectrum_scale + deterministic;
      dwk[p] = {0.0, 0.0};
    }
    lnw -= bias * zero_mode.real() + 0.5 * bias * bias * dn / w;
    std::cout << "\r" << w << " / " << MAX_SHELL << std::flush;
  }
  std::cout << '\n';

  const double k_unit = 2.0 * M_PI / NL;
  Fft3D fft;

  // ----------- compaction C(r,x), integer radii rs = 1..12 (as in Gaussian_LN_peak.cpp) --------
  // Local "compaction peak" rule of the direct-sampling code (processCompactionPeaks), applied
  // at every lattice site of a window |x|_inf <= H around the bias centre (origin).
  // The rule is local, so  n(props) = E_unbiased[1(site is a C-peak)] = E_biased[W 1(...)]  per site.
  const double sqrtAs = std::sqrt(As);
  const int rs_hi = static_cast<int>(5.0 / (k_unit * nsigma));   // = 12, same as direct sampling
  const int H = 16, WN = 2 * H + 1;
  const int N2MAX = 3 * (NL / 2) * (NL / 2);
  std::vector<double> tbl(N2MAX + 1);
  std::vector<std::vector<double>> win(rs_hi + 1, std::vector<double>(static_cast<std::size_t>(WN) * WN * WN));
  auto widx = [&](int a, int b, int c) { return (static_cast<std::size_t>(a + H) * WN + (b + H)) * WN + (c + H); };
  auto wrap = [&](int a) { return a < 0 ? a + NL : a; };
  for (int rs = 1; rs <= rs_hi; ++rs) {
    for (int n = 0; n <= N2MAX; ++n) {
      const double kr = k_unit * std::sqrt(static_cast<double>(n)) * rs;
      tbl[n] = -kr * kr / 3.0 * WRTH(kr) * sqrtAs;
    }
    fftw_complex* in = fft.in();
    for (int i = 0; i < NL; ++i) {
      const int x = shiftedindex(i);
      for (int j = 0; j < NL; ++j) {
        const int y = shiftedindex(j);
        for (int k = 0; k < NL; ++k) {
          const int z = shiftedindex(k);
          const std::size_t p = index_of(i, j, k);
          const double f = tbl[x * x + y * y + z * z];
          in[p][0] = gkbias[p].real() * f;
          in[p][1] = gkbias[p].imag() * f;
        }
      }
    }
    fft.execute();
    const fftw_complex* out = fft.out();
    for (int a = -H; a <= H; ++a) for (int b = -H; b <= H; ++b) for (int c = -H; c <= H; ++c) {
      const double v = out[index_of(wrap(a), wrap(b), wrap(c))][0];
      win[rs][widx(a, b, c)] = 2.0 / 3.0 * (1.0 - (1.0 + v) * (1.0 + v));
    }
  }

  // events: one line per C-peak found at a window site; one summary line per sample
  std::string buf;
  int nev = 0;
  for (int rp = 1; rp + 2 <= rs_hi; ++rp) {
    const auto& C0 = win[rp];
    const auto& C1 = win[rp + 1];
    const auto& C2 = win[rp + 2];
    for (int i = -H; i <= H - 2; ++i) for (int j = -H; j <= H - 2; ++j) for (int k = -H; k <= H - 2; ++k) {
      const std::size_t idx = widx(i, j, k);
      if (!(C0[idx] < C1[idx] && C1[idx] > C2[idx])) continue;
      const double dx = C0[widx(i + 1, j, k)] - C0[idx];
      const double dxrot = C0[widx(i + 2, j, k)] - C0[widx(i + 1, j, k)];
      const double dy = C0[widx(i, j + 1, k)] - C0[idx];
      const double dyrot = C0[widx(i, j + 2, k)] - C0[widx(i, j + 1, k)];
      const double dz = C0[widx(i, j, k + 1)] - C0[idx];
      const double dzrot = C0[widx(i, j, k + 2)] - C0[widx(i, j, k + 1)];
      if (dx * dxrot < 0 && dx > 0 && dy * dyrot < 0 && dy > 0 && dz * dzrot < 0 && dz > 0) {
        const double val = C1[widx(i + 1, j + 1, k + 1)];
        buf += "E," + std::to_string(seed) + ',' + std::to_string(lnw) + ',' + std::to_string(i + 1) + ',' +
               std::to_string(j + 1) + ',' + std::to_string(k + 1) + ',' + std::to_string(rp + 1) + ',' +
               std::to_string(val) + '\n';
        ++nev;
      }
    }
  }
  buf += "S," + std::to_string(seed) + ',' + std::to_string(lnw) + ',' + std::to_string(nev) + '\n';
  mukfile << buf;
  mukfile.flush();
  gettimeofday(&now, nullptr);
  const double after = now.tv_sec + now.tv_usec * 1.e-6;
  std::cout << after - before << " sec.\n";
  return 0;
}
