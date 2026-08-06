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

constexpr int NL = 256;
constexpr int nsigma = 16;
constexpr double As = 3.5e-3;
constexpr int nbias = 16;
constexpr double dlnn = 0.1;
constexpr double biascoeff = 12.5;
constexpr double s2 = 0.01;
const std::string s2value = "0,01";
constexpr double dn = 1.0;

constexpr std::size_t N3 = static_cast<std::size_t>(NL) * NL * NL;
// shifted indices span [-127, 128], so the furthest populated shell is 222.
constexpr int MAX_SHELL = 222;
const std::string mukfilename = "data/LN_muk_" + s2value + "_" +
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

// Reuses both buffers and the FFTW plan.  The original code recreated all three
// for every transform, which is particularly expensive for the radius scan.
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
  if (argc != 2) {
    std::cerr << "Specify the seed correctly.\n";
    return 1;
  }

  timeval now{};
  gettimeofday(&now, nullptr);
  const double before = now.tv_sec + now.tv_usec * 1.e-6;

  const int seed = std::atoi(argv[1]);
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
  Grid Dgk(N3);
  for (std::size_t p = 0; p < N3; ++p) {
    const double k2 = k_unit * k_unit * modes.norm[p] * modes.norm[p];
    Dgk[p] = gkbias[p] * k2;
  }

  Fft3D fft;
  const Grid gxbias = fft.transform(gkbias);
  Grid Dgx = fft.transform(Dgk);
  const auto peak = std::max_element(Dgx.begin(), Dgx.end(),
      [](const auto& a, const auto& b) { return a.real() < b.real(); });
  const std::size_t peak_index = static_cast<std::size_t>(std::distance(Dgx.begin(), peak));
  const int imax = static_cast<int>(peak_index / (NL * NL));
  const int jmax = static_cast<int>((peak_index / NL) % NL);
  const int kmax = static_cast<int>(peak_index % NL);
  const double mu2 = peak->real();
  for (std::size_t p = 0; p < N3; ++p) {
    const double k2 = k_unit * k_unit * modes.norm[p] * modes.norm[p];
    Dgk[p] *= k2;
  }
  const double k3 = std::sqrt(fft.transform(Dgk)[peak_index].real() / mu2);
  Grid().swap(Dgk);
  Grid().swap(Dgx);

  double Cmax = 0.0;
  int rsmax = 0;
  Grid rzpk(N3);
  const int rs_limit = static_cast<int>(10.0 / (k_unit * nsigma));
  for (int rs = 1; rs <= rs_limit; ++rs) {
    for (std::size_t p = 0; p < N3; ++p) {
      const double kr = k_unit * modes.norm[p] * rs;
      rzpk[p] = gkbias[p] * (-kr * kr / 3.0 * WRTH(kr) * std::sqrt(As));
    }
    const Grid rzpx = fft.transform(rzpk);
    const double compaction = 2.0 / 3.0 * (1.0 - std::pow(1.0 + rzpx[peak_index].real(), 2));
    if (compaction > Cmax) {
      Cmax = compaction;
      rsmax = rs;
    }
  }

  const int nxm = shiftedindex(imax);
  const int nym = shiftedindex(jmax);
  const int nzm = shiftedindex(kmax);
  int count = 0;
  double zetam = 0.0;
  for (int i = 0; i < NL; ++i) for (int j = 0; j < NL; ++j) for (int k = 0; k < NL; ++k) {
    const int dx = shiftedindex(i) - nxm;
    const int dy = shiftedindex(j) - nym;
    const int dz = shiftedindex(k) - nzm;
    const double dr = std::sqrt(static_cast<double>(dx * dx + dy * dy + dz * dz));
    if (std::fabs(dr - rsmax) < 0.5) {
      zetam += gxbias[index_of(i, j, k)].real() * std::sqrt(As);
      ++count;
    }
  }
  zetam /= count;

  mukfile << seed << ',' << mu2 << ',' << k3 << ',' << k_unit * nsigma * rsmax
          << ',' << zetam << ',' << Cmax << ',' << lnw << '\n';
  gettimeofday(&now, nullptr);
  const double after = now.tv_sec + now.tv_usec * 1.e-6;
  std::cout << after - before << " sec.\n";
  return 0;
}
