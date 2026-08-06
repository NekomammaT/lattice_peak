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
constexpr double biascoeff = 10;
constexpr double s2 = 0.1;
const std::string s2value = "0,1";
constexpr double dn = 1.0;
constexpr std::size_t N3 = static_cast<std::size_t>(NL) * NL * NL;
// shifted indices span [-127, 128], so the furthest populated shell is 222.
constexpr int MAX_SHELL = 222;

const std::string mapfileprefix = "data/LN_map_" + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_";
const std::string biasedfileprefix = "data/LN_biased_" + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_" + std::to_string(nbias) + "_";
const std::string laplacianfileprefix = "data/LN_laplacian_" + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_" + std::to_string(nbias) + "_";
const std::string powerfileprefix = "data/LN_power_" + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_";
const std::string Cfileprefix = "data/LN_compaction_" + s2value + "_" + std::to_string(NL) + "_" + std::to_string(nsigma) + "_" + std::to_string(nbias) + "_";

using Grid = std::vector<std::complex<double>>;

inline std::size_t index_of(int i, int j, int k) { return (static_cast<std::size_t>(i) * NL + j) * NL + k; }
inline int shiftedindex(int n) { return n <= NL / 2 ? n : n - NL; }

inline bool realpoint(int x, int y, int z)
{
  return (x == 0 || x == NL / 2) && (y == 0 || y == NL / 2) && (z == 0 || z == NL / 2);
}

inline bool complexpoint(int nx, int ny, int nz)
{
  const int x = shiftedindex(nx), y = shiftedindex(ny), z = shiftedindex(nz);
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

inline double powerspectrum(int w)
{
  const double d = std::log(static_cast<double>(w)) - std::log(static_cast<double>(nsigma));
  return std::exp(-d * d / (2.0 * s2)) / std::sqrt(2.0 * M_PI * s2);
}

inline double BN(int w)
{
  const double d = std::log(static_cast<double>(w)) - std::log(static_cast<double>(nbias));
  return biascoeff * std::exp(-d * d / (2.0 * dlnn)) / std::sqrt(2.0 * M_PI * dlnn);
}

inline double WRTH(double z) { return z == 0.0 ? 1.0 : 3.0 * (std::sin(z) - z * std::cos(z)) / (z * z * z); }

class Fft3D {
 public:
  Fft3D()
  {
    in_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N3));
    out_ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * N3));
    if (!in_ || !out_) throw std::bad_alloc();
    plan_ = fftw_plan_dft_3d(NL, NL, NL, in_, out_, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!plan_) throw std::runtime_error("could not create FFTW plan");
  }
  ~Fft3D() { if (plan_) fftw_destroy_plan(plan_); fftw_free(in_); fftw_free(out_); }
  Fft3D(const Fft3D&) = delete;
  Fft3D& operator=(const Fft3D&) = delete;

  Grid transform(const Grid& input)
  {
    for (std::size_t p = 0; p < N3; ++p) { in_[p][0] = input[p].real(); in_[p][1] = input[p].imag(); }
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
  std::vector<std::vector<std::size_t>> shells, independent, reflected;
  std::vector<std::size_t> reflection_source;
  std::vector<double> norm;
};

FourierModes make_fourier_modes()
{
  FourierModes modes{std::vector<std::vector<std::size_t>>(MAX_SHELL + 1),
                     std::vector<std::vector<std::size_t>>(MAX_SHELL + 1),
                     std::vector<std::vector<std::size_t>>(MAX_SHELL + 1),
                     std::vector<std::size_t>(N3), std::vector<double>(N3)};
  // Keep the original i-j-k order, including the distinct reflection pass.
  for (int i = 0; i < NL; ++i) for (int j = 0; j < NL; ++j) for (int k = 0; k < NL; ++k) {
    const std::size_t p = index_of(i, j, k);
    const int x = shiftedindex(i), y = shiftedindex(j), z = shiftedindex(k);
    modes.norm[p] = std::sqrt(static_cast<double>(x * x + y * y + z * z));
    const int w = static_cast<int>(std::floor(modes.norm[p] + 0.5));
    modes.shells[w].push_back(p);
    if (realpoint(i, j, k) || complexpoint(i, j, k)) modes.independent[w].push_back(p);
  }
  for (int i = 0; i < NL; ++i) for (int j = 0; j < NL; ++j) for (int k = 0; k < NL; ++k) {
    if (realpoint(i, j, k) || complexpoint(i, j, k)) continue;
    const std::size_t p = index_of(i, j, k);
    const int w = static_cast<int>(std::floor(modes.norm[p] + 0.5));
    modes.reflected[w].push_back(p);
    modes.reflection_source[p] = index_of(i == 0 ? 0 : NL - i, j == 0 ? 0 : NL - j, k == 0 ? 0 : NL - k);
  }
  return modes;
}

Grid make_random_field(const FourierModes& modes, std::mt19937& engine,
                       std::normal_distribution<>& dist, bool biased)
{
  Grid field(N3, {0.0, 0.0});
  Grid shell_values(N3, {0.0, 0.0});
  for (int w = 1; w <= MAX_SHELL; ++w) {
    const auto& shell = modes.shells[w];
    const double random_scale = std::sqrt(dn / w / static_cast<double>(shell.size()));
    for (std::size_t p : modes.independent[w]) {
      const int k = static_cast<int>(p % NL), j = static_cast<int>((p / NL) % NL), i = static_cast<int>(p / (static_cast<std::size_t>(NL) * NL));
      shell_values[p] = realpoint(i, j, k) ? std::complex<double>(dist(engine), 0.0)
                                             : std::complex<double>(dist(engine), dist(engine)) / std::sqrt(2.0);
    }
    for (std::size_t p : modes.reflected[w]) shell_values[p] = std::conj(shell_values[modes.reflection_source[p]]);
    const double spectrum_scale = std::sqrt(powerspectrum(w));
    const std::complex<double> deterministic = biased ? std::complex<double>(BN(w) * dn / w / shell.size() * spectrum_scale, 0.0) : std::complex<double>(0.0, 0.0);
    for (std::size_t p : shell) {
      field[p] += shell_values[p] * random_scale * spectrum_scale + deterministic;
      shell_values[p] = {0.0, 0.0};
    }
    std::cout << "\r" << w << " / " << MAX_SHELL << std::flush;
  }
  std::cout << '\n';
  return field;
}

void write_real_grid(std::ofstream& file, const Grid& grid)
{
  for (std::size_t p = 0; p < N3; ++p) { file << grid[p].real(); if (p + 1 != N3) file << ','; }
  file << '\n';
}

int main(int argc, char* argv[])
{
  if (argc != 2) { std::cerr << "Specify the noise file number correctly.\n"; return 1; }
  timeval now{};
  gettimeofday(&now, nullptr);
  const double before = now.tv_sec + now.tv_usec * 1.e-6;
  const uint32_t seed = static_cast<uint32_t>(std::atoi(argv[1]));
  std::cout << "seed = " << seed << '\n';
  std::mt19937 engine(std::hash<int>{}(seed));
  std::normal_distribution<> dist(0.0, 1.0);
  const FourierModes modes = make_fourier_modes();
  Fft3D fft;

  /*
  // These grids are no longer needed once the unbiased outputs are written.
  // Keeping this work in a scope releases roughly 768 MiB before the biased run.
  {
    std::ofstream mapfile(mapfileprefix + std::to_string(seed) + ".csv");
    const Grid gk = make_random_field(modes, engine, dist, false);
    const Grid gx = fft.transform(gk);
    write_real_grid(mapfile, gx);
    std::cout << "Exported to " << mapfileprefix + std::to_string(seed) + ".csv\n";

    std::ofstream powerfile(powerfileprefix + std::to_string(seed) + ".csv");
    const Grid gkp = fft.transform(gx);
    const int powerlistlength = static_cast<int>(std::round(std::sqrt(3.0) * (NL - 1))) + 1;
    std::vector<double> calPg(powerlistlength, 0.0);
    const double inverse_n3 = 1.0 / N3;
    for (std::size_t p = 0; p < N3; ++p) calPg[static_cast<int>(std::round(modes.norm[p]))] += std::norm(gkp[p] * inverse_n3);
    for (int w = 0; w < powerlistlength; ++w) { powerfile << calPg[w] * w; if (w + 1 != powerlistlength) powerfile << ','; }
    powerfile << '\n';
    std::cout << "Exported to " << powerfileprefix + std::to_string(seed) + ".csv\n";
  }
    */

  std::ofstream biasedfile(biasedfileprefix + std::to_string(seed) + ".csv");
  std::ofstream laplacianfile(laplacianfileprefix + std::to_string(seed) + ".csv");
  std::ofstream compactionfile(Cfileprefix + std::to_string(seed) + ".csv");
  const Grid gkbias = make_random_field(modes, engine, dist, true);
  const double k_unit = 2.0 * M_PI / NL;
  Grid Dgk(N3);
  for (std::size_t p = 0; p < N3; ++p) Dgk[p] = gkbias[p] * (k_unit * k_unit * modes.norm[p] * modes.norm[p]);
  const Grid gxbias = fft.transform(gkbias);
  Grid Dgx = fft.transform(Dgk);
  write_real_grid(biasedfile, gxbias);
  write_real_grid(laplacianfile, Dgx);
  std::cout << "Exported to " << biasedfileprefix + std::to_string(seed) + ".csv\n";
  std::cout << "Exported to " << laplacianfileprefix + std::to_string(seed) + ".csv\n";

  const auto peak = std::max_element(Dgx.begin(), Dgx.end(), [](const auto& a, const auto& b) { return a.real() < b.real(); });
  const std::size_t peak_index = static_cast<std::size_t>(std::distance(Dgx.begin(), peak));
  const double mu2 = peak->real();
  const int imax = static_cast<int>(peak_index / (NL * NL)), jmax = static_cast<int>((peak_index / NL) % NL), kmax = static_cast<int>(peak_index % NL);
  for (std::size_t p = 0; p < N3; ++p) Dgk[p] *= k_unit * k_unit * modes.norm[p] * modes.norm[p];
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
    const double compaction = 2.0 / 3.0 * (1.0 - std::pow(1.0 + fft.transform(rzpk)[peak_index].real(), 2));
    if (compaction > Cmax) { Cmax = compaction; rsmax = rs; }
    compactionfile << rs << ',' << compaction << '\n';
  }
  const int nxm = shiftedindex(imax), nym = shiftedindex(jmax), nzm = shiftedindex(kmax);
  int count = 0;
  double zetam = 0.0;
  for (int i = 0; i < NL; ++i) for (int j = 0; j < NL; ++j) for (int k = 0; k < NL; ++k) {
    const int dx = shiftedindex(i) - nxm, dy = shiftedindex(j) - nym, dz = shiftedindex(k) - nzm;
    if (std::fabs(std::sqrt(static_cast<double>(dx * dx + dy * dy + dz * dz)) - rsmax) < 0.5) { zetam += gxbias[index_of(i, j, k)].real() * std::sqrt(As); ++count; }
  }
  zetam /= count;
  compactionfile << mu2 << ',' << k3 << ',' << k_unit * nsigma * rsmax << ',' << zetam << '\n';
  std::cout << "Exported to " << Cfileprefix + std::to_string(seed) + ".csv\n";
  gettimeofday(&now, nullptr);
  std::cout << now.tv_sec + now.tv_usec * 1.e-6 - before << " sec.\n";
  return 0;
}
