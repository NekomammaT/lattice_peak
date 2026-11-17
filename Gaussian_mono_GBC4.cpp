#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <sys/time.h>
#include <vector>

#include <fftw3.h>

// ============================================================================
// Gaussian_mono_GBC4.cpp  (monochromatic spectrum, As = 5e-3): importance-sampling program,
// the monochromatic counterpart of Gaussian_LN_GBC4.cpp (same output format and window H = 8).
//
// Field: one shell |n| in [nsigma-1/2, nsigma+1/2) (exactly the field of Gaussian_mono.cpp /
// Gaussian_mono_peak.cpp: same draw order, same 1/sqrt(#modes) normalisation, so the unbiased
// part of seed s is identical to the one of the direct-sampling programs).
// Bias: the shell zero mode is shifted by `bc` (Gaussian_mono.cpp: bias = 8), spread uniformly
// over the shell modes (bc/#modes each); likelihood ratio  lnW = -bc*z0 - bc^2/2,
// z0 = zero mode (= field at the origin) of the *unbiased* noise.  Bias centre = origin.
//
// In the window |x|_inf <= H (H = 8) around the origin it records
//   P : every Laplacian peak (direct-sampling peak rule on D2g, mu2 > 0):
//       position, mu2, k3 = sqrt(D2D2g/D2g), and the maximum over rs = 1..12 of the
//       compaction centred on the peak: Cmax_c (at the peak site), rm_c, and
//       Cmax_n (max over the 3x3x3 neighbourhood at each rs), rm_n
//   E : every compaction peak (direct-sampling local rule) at radius rp+1 with value val
//   S : one summary line per sample.
// File: data/mono_muk_<NL>_<nsigma>_bc<bc>.csv, appended (one write per sample);
// usage ./Gaussian_mono_GBC4 seed bc
//   P,seed,lnW,x,y,z,mu2,k3,Cmax_c,rm_c,Cmax_n,rm_n
//   E,seed,lnW,x,y,z,rp+1,val
//   S,seed,lnW,nP,nE
// Both site-count (site0) and k3^3-weighted estimators can be formed from these lines
// (analysis exactly as for LN_muk_*; sigma2 = k^2 of the shell, computed from the data
// of the direct program or kpk^2 up to the shell-width effect).
// Single-threaded (job-parallel over seeds), like Gaussian_LN_GBC4.
// ============================================================================

constexpr int NL = 256;
constexpr int nsigma = 16;
constexpr double As = 5e-3;
constexpr double dn = 1.0;     // shell thickness
double biascoeff = 8.0;        // overridden by argv[2]

constexpr std::size_t N3 = static_cast<std::size_t>(NL) * NL * NL;
const std::string mukfileprefix = "data/mono_muk_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_bc";

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

// same shell definition as innsigma() of Gaussian_mono*.cpp
inline bool innsigma(int nx, int ny, int nz, double wavenumber)
{
  const int x = shiftedindex(nx), y = shiftedindex(ny), z = shiftedindex(nz);
  const double nt = std::sqrt(static_cast<double>(x * x + y * y + z * z));
  return (wavenumber - dn / 2. <= nt && nt < wavenumber + dn / 2.);
}

inline double WRTH(double z)
{
  return z == 0.0 ? 1.0 : 3.0 * (std::sin(z) - z * std::cos(z)) / (z * z * z);
}

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
  void clear_in() { std::memset(in_, 0, sizeof(fftw_complex) * N3); }
  void execute() { fftw_execute(plan_); }
 private:
  fftw_complex* in_ = nullptr;
  fftw_complex* out_ = nullptr;
  fftw_plan plan_ = nullptr;
};

int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: Gaussian_mono_GBC4 <seed> <biascoeff>\n";
    return 1;
  }

  timeval now{};
  gettimeofday(&now, nullptr);
  const double before = now.tv_sec + now.tv_usec * 1.e-6;

  const int seed = std::atoi(argv[1]);
  biascoeff = std::atof(argv[2]);
  const std::string mukfilename = mukfileprefix + argv[2] + ".csv";
  std::cout << "seed = " << seed << '\n';
  std::mt19937 engine(std::hash<int>{}(seed));
  std::normal_distribution<> dist(0.0, 1.0);
  std::ofstream mukfile(mukfilename, std::ios::app);

  // ----------- shell modes (i-j-k traversal order, as in Gaussian_mono.cpp dwk()) -----------
  std::vector<std::size_t> shell, independent, reflected, reflection_source;
  std::vector<int> sx, sy, sz;   // shifted integer wave vector of every shell mode
  for (int i = 0; i < NL; ++i) for (int j = 0; j < NL; ++j) for (int k = 0; k < NL; ++k) {
    if (!innsigma(i, j, k, nsigma)) continue;
    const std::size_t p = index_of(i, j, k);
    shell.push_back(p);
    sx.push_back(shiftedindex(i)); sy.push_back(shiftedindex(j)); sz.push_back(shiftedindex(k));
    if (realpoint(i, j, k) || complexpoint(i, j, k)) independent.push_back(p);
    else {
      reflected.push_back(p);
      const int ip = i == 0 ? 0 : NL - i, jp = j == 0 ? 0 : NL - j, kp = k == 0 ? 0 : NL - k;
      reflection_source.push_back(index_of(ip, jp, kp));
    }
  }
  const double nmodes = static_cast<double>(shell.size());

  Grid gk(N3, {0.0, 0.0});
  for (const std::size_t p : independent) {
    const int k = static_cast<int>(p % NL);
    const int j = static_cast<int>((p / NL) % NL);
    const int i = static_cast<int>(p / (static_cast<std::size_t>(NL) * NL));
    gk[p] = realpoint(i, j, k)
        ? std::complex<double>(dist(engine), 0.0)
        : std::complex<double>(dist(engine), dist(engine)) / std::sqrt(2.0);
  }
  for (std::size_t q = 0; q < reflected.size(); ++q) gk[reflected[q]] = std::conj(gk[reflection_source[q]]);
  std::complex<double> zero_mode(0.0, 0.0);
  const double norm = std::sqrt(nmodes);
  for (const std::size_t p : shell) { gk[p] /= norm; zero_mode += gk[p]; }
  // gx[0] of Gaussian_mono.cpp = zero_mode;  bias shifts it by biascoeff
  const double lnw = -biascoeff * zero_mode.real() - 0.5 * biascoeff * biascoeff;
  Grid gkbias = std::move(gk);
  for (const std::size_t p : shell) gkbias[p] += biascoeff / nmodes;

  const double k_unit = 2.0 * M_PI / NL;
  Fft3D fft;

  const int H = 8, WN = 2 * H + 1;   // window |x|_inf <= H around the bias centre (origin)
  const int WN3 = WN * WN * WN;
  auto widx = [&](int a, int b, int c) { return (static_cast<std::size_t>(a + H) * WN + (b + H)) * WN + (c + H); };
  auto wrap = [&](int a) { return a < 0 ? a + NL : a; };

  // ----------- Laplacian D2g and D4g = k^4 g on the window (for Laplacian peaks) -----------
  std::vector<double> winD2(WN3), winD4(WN3);
  for (int pw = 1; pw <= 2; ++pw) {
    fft.clear_in();
    fftw_complex* in = fft.in();
    for (std::size_t q = 0; q < shell.size(); ++q) {
      const std::size_t p = shell[q];
      const double k2 = k_unit * k_unit * static_cast<double>(sx[q] * sx[q] + sy[q] * sy[q] + sz[q] * sz[q]);
      const double f = (pw == 1) ? k2 : k2 * k2;
      in[p][0] = gkbias[p].real() * f;
      in[p][1] = gkbias[p].imag() * f;
    }
    fft.execute();
    const fftw_complex* out = fft.out();
    std::vector<double>& dst = (pw == 1) ? winD2 : winD4;
    for (int a = -H; a <= H; ++a) for (int b = -H; b <= H; ++b) for (int c = -H; c <= H; ++c)
      dst[widx(a, b, c)] = out[index_of(wrap(a), wrap(b), wrap(c))][0];
  }
  struct LPk { int x, y, z; double mu2, k3, Cc, Cn; int rc, rn; };
  std::vector<LPk> lpk;
  for (int i = -H; i <= H - 2; ++i) for (int j = -H; j <= H - 2; ++j) for (int k = -H; k <= H - 2; ++k) {
    const std::size_t idx = widx(i, j, k);
    const double dx = winD2[widx(i + 1, j, k)] - winD2[idx];
    const double dxrot = winD2[widx(i + 2, j, k)] - winD2[widx(i + 1, j, k)];
    const double dy = winD2[widx(i, j + 1, k)] - winD2[idx];
    const double dyrot = winD2[widx(i, j + 2, k)] - winD2[widx(i, j + 1, k)];
    const double dz = winD2[widx(i, j, k + 1)] - winD2[idx];
    const double dzrot = winD2[widx(i, j, k + 2)] - winD2[widx(i, j, k + 1)];
    if (dx * dxrot < 0 && dx > 0 && dy * dyrot < 0 && dy > 0 && dz * dzrot < 0 && dz > 0) {
      const std::size_t ip = widx(i + 1, j + 1, k + 1);
      const double mu2 = winD2[ip];
      if (mu2 <= 0) continue;
      lpk.push_back({i + 1, j + 1, k + 1, mu2, std::sqrt(winD4[ip] / mu2), -1e300, -1e300, 0, 0});
    }
  }

  // ----------- compaction C(r,x), integer radii rs = 1..12 (as in Gaussian_mono_peak.cpp) --------
  // Local "compaction peak" rule of the direct-sampling code, applied at every lattice site of
  // the window; n(props) = E_unbiased[1(site is a C-peak)] = E_biased[W 1(...)] per site.
  const double sqrtAs = std::sqrt(As);
  const int rs_hi = static_cast<int>(5.0 / (k_unit * nsigma));   // = 12, same as direct sampling
  std::vector<std::vector<double>> win(rs_hi + 1, std::vector<double>(static_cast<std::size_t>(WN3)));
  for (int rs = 1; rs <= rs_hi; ++rs) {
    fft.clear_in();
    fftw_complex* in = fft.in();
    for (std::size_t q = 0; q < shell.size(); ++q) {
      const std::size_t p = shell[q];
      const double kr = k_unit * std::sqrt(static_cast<double>(sx[q] * sx[q] + sy[q] * sy[q] + sz[q] * sz[q])) * rs;
      const double f = -kr * kr / 3.0 * WRTH(kr) * sqrtAs;
      in[p][0] = gkbias[p].real() * f;
      in[p][1] = gkbias[p].imag() * f;
    }
    fft.execute();
    const fftw_complex* out = fft.out();
    for (int a = -H; a <= H; ++a) for (int b = -H; b <= H; ++b) for (int c = -H; c <= H; ++c) {
      const double v = out[index_of(wrap(a), wrap(b), wrap(c))][0];
      win[rs][widx(a, b, c)] = 2.0 / 3.0 * (1.0 - (1.0 + v) * (1.0 + v));
    }
  }

  // Laplacian peaks: Cmax over rs = 1..rs_hi of C at the peak site and of its 3x3x3 neighbourhood max
  std::string buf;
  for (LPk& q : lpk) {
    for (int rs = 1; rs <= rs_hi; ++rs) {
      const auto& C = win[rs];
      const double c0 = C[widx(q.x, q.y, q.z)];
      if (c0 > q.Cc) { q.Cc = c0; q.rc = rs; }
      double m = -1e300;
      for (int a = -1; a <= 1; ++a) for (int b = -1; b <= 1; ++b) for (int c = -1; c <= 1; ++c)
        m = std::max(m, C[widx(q.x + a, q.y + b, q.z + c)]);
      if (m > q.Cn) { q.Cn = m; q.rn = rs; }
    }
    buf += "P," + std::to_string(seed) + ',' + std::to_string(lnw) + ',' + std::to_string(q.x) + ',' +
           std::to_string(q.y) + ',' + std::to_string(q.z) + ',' + std::to_string(q.mu2) + ',' +
           std::to_string(q.k3) + ',' + std::to_string(q.Cc) + ',' + std::to_string(q.rc) + ',' +
           std::to_string(q.Cn) + ',' + std::to_string(q.rn) + '\n';
  }
  // events: one line per C-peak found at a window site; one summary line per sample
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
  buf += "S," + std::to_string(seed) + ',' + std::to_string(lnw) + ',' + std::to_string(lpk.size()) + ',' + std::to_string(nev) + '\n';
  mukfile << buf;
  mukfile.flush();
  gettimeofday(&now, nullptr);
  const double after = now.tv_sec + now.tv_usec * 1.e-6;
  std::cout << after - before << " sec.\n";
  return 0;
}
