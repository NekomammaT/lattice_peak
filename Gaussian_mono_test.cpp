#include <iostream>
#include <fstream>
#include <random>
#include <sys/time.h>
#include "fft.hpp"

std::vector<std::vector<std::vector<std::complex<double>>>> dwk(int wavenumber, double bias, int seed);
int shiftedindex(int n); // shifted index
bool innsigma(int nx, int ny, int nz, double wavenumber); // judge if point is in nsigma sphere shell
bool realpoint(int nx, int ny, int nz);                   // judge real point
bool complexpoint(int nx, int ny, int nz);                // judge independent complex point

// useful macro
#define LOOP                     \
  for (int i = 0; i < NL; i++)   \
    for (int j = 0; j < NL; j++) \
      for (int k = 0; k < NL; k++)

// random distribution
std::normal_distribution<> dist(0., 1.);

// imaginary unit
const std::complex<double> II(0, 1);

// parameters
const int NL = pow(2, 8); // Box size NL
const int nsigma = pow(2, 4);
const double dn = 1; // Thickness of nsigma sphere shell
const double bias = 10;
const std::string mapfileprefix = "data/mono_map_";
const std::string biasedfileprefix = "data/mono_biased_";
const std::string laplacianfileprefix = "data/mono_laplacian_";
const std::string powerfileprefix = "data/mono_power_";

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
  std::ofstream mapfile(mapfileprefix + std::to_string(seed) + ".csv");
  std::ofstream biasedfile(biasedfileprefix + std::to_string(seed) + ".csv");
  std::ofstream laplacianfile(laplacianfileprefix + std::to_string(seed) + ".csv");
  std::ofstream powerfile(powerfileprefix + std::to_string(seed) + ".csv");

  // ----------- unbiased map -----------
  std::vector<std::vector<std::vector<std::complex<double>>>> gk = dwk(nsigma, 0., seed);
  std::vector<std::vector<std::vector<std::complex<double>>>> gx = fftw(gk);
  
  LOOP
  {
    mapfile << gx[i][j][k].real() << ' '; 

    if (i != NL-1 || j != NL-1 || k != NL-1) {
      mapfile << ','; 
    }
  }
  mapfile << std::endl;

  std::cout << "Exported to " << mapfileprefix + std::to_string(seed) + ".csv" << std::endl;

  // ----------- power spectrum ----------
  std::vector<std::vector<std::vector<std::complex<double>>>> gkp = fftw(gx);
  int powerlistlength = std::round(sqrt(3)*(NL-1)) + 1;
  std::vector<std::vector<double>> powerdata(powerlistlength, std::vector<double>{});
  double sigma1sq = 0;
  double sigma2sq = 0;
  double sigma4sq = 0;

  LOOP
  {
    int nxt = shiftedindex(i);
    int nyt = shiftedindex(j);
    int nzt = shiftedindex(k);
    double ntnorm = sqrt(nxt*nxt+nyt*nyt+nzt*nzt);

    powerdata[std::round(ntnorm)].push_back(std::norm(gkp[i][j][k]));
    sigma1sq += pow(ntnorm,2)*std::norm(gkp[i][j][k]) / pow(NL, 6);
    sigma2sq += pow(ntnorm,4)*std::norm(gkp[i][j][k]) / pow(NL, 6);
    sigma4sq += pow(ntnorm,8)*std::norm(gkp[i][j][k]) / pow(NL, 6);
  }
  std::cout << "sigma1sq = " << sigma1sq << std::endl;
  std::cout << "sigma2sq = " << sigma2sq << std::endl;
  std::cout << "sigma4sq = " << sigma4sq << std::endl;

  std::vector<double> calPg(powerlistlength, 0);
  for (int i = 0; i < calPg.size(); i++)
  {
    if (powerdata[i].size() != 0)
    {
      for (int j = 0; j < powerdata[i].size(); j++)
      {
        calPg[i] += powerdata[i][j];
      }
      calPg[i] *= i / pow(NL, 6);
    }

    powerfile << calPg[i] << ' ';
    if (i != calPg.size()-1) powerfile << ',';
  }
  powerfile << std::endl;

  std::cout << "Exported to " << powerfileprefix + std::to_string(seed) + ".csv" << std::endl;
  
  // ----------- biased map -----------
  std::vector<std::vector<std::vector<std::complex<double>>>> gkbias = dwk(nsigma, bias, seed);
  std::vector<std::vector<std::vector<std::complex<double>>>> Dgk = gkbias;
  std::vector<std::vector<std::vector<std::complex<double>>>> DDgk = gkbias;
  LOOP
  {
    int nxt = shiftedindex(i);
    int nyt = shiftedindex(j);
    int nzt = shiftedindex(k);
    double ntnorm = sqrt(nxt*nxt+nyt*nyt+nzt*nzt);

    Dgk[i][j][k] *= pow(ntnorm,2);
    DDgk[i][j][k] *= pow(ntnorm,4);
  }
  std::vector<std::vector<std::vector<std::complex<double>>>> gxbias = fftw(gkbias);
  std::vector<std::vector<std::vector<std::complex<double>>>> Dgx = fftw(Dgk);
  std::vector<std::vector<std::vector<std::complex<double>>>> DDgx = fftw(DDgk);

  std::vector<double> Dgx1d(NL*NL*NL, 0);
  LOOP
  {
    Dgx1d[i*NL*NL+j*NL+k] = Dgx[i][j][k].real();

    biasedfile << gxbias[i][j][k].real() << ' ';
    laplacianfile << Dgx[i][j][k].real() << ' ';
    if (i != NL-1 || j != NL-1 || k != NL-1) {
      biasedfile << ',';
      laplacianfile << ',';
    }
  }
  biasedfile << std::endl;
  laplacianfile << std::endl;
  std::cout << "Exported to " << biasedfileprefix + std::to_string(seed) + ".csv" << std::endl;
  std::cout << "Exported to " << laplacianfileprefix + std::to_string(seed) + ".csv" << std::endl;

  auto iter = std::max_element(Dgx1d.begin(), Dgx1d.end());
  size_t index = std::distance(Dgx1d.begin(), iter);
  std::cout << "Maximal value is " << *iter << " at " << index << std::endl;

  int imax = index / (NL * NL);
  int jmax = (index - imax * NL * NL) / NL;
  int kmax = index - imax * NL * NL - jmax * NL;
  double mu2 = Dgx[imax][jmax][kmax].real() * sigma1sq/sigma2sq;
  double k3 = sqrt(DDgx[imax][jmax][kmax].real() / Dgx[imax][jmax][kmax].real() * sqrt(sigma2sq/sigma4sq));
  std::cout << "mu2 is " << mu2 << " at (" << imax << ", " << jmax << ", " << kmax << ")" << std::endl;
  std::cout << "k3 is " << k3 << std::endl;

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------

  return 0;
}

std::vector<std::vector<std::vector<std::complex<double>>>> dwk(int wavenumber, double bias, int seed)
{
  std::vector<std::vector<std::vector<std::complex<double>>>> dwk(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));

  int count = 0;
  std::mt19937 engine(std::hash<int>{}(seed));

  LOOP
  {
    if (innsigma(i, j, k, wavenumber))
    {
      if (realpoint(i, j, k))
      {
        dwk[i][j][k] = dist(engine);
        count++;
      }
      else if (complexpoint(i, j, k))
      {
        dwk[i][j][k] = (dist(engine) + II * dist(engine)) / sqrt(2);
        count++;
      }
    }
  }

  // reflection
  int ip, jp, kp; // reflected index
  LOOP
  {
    if (innsigma(i, j, k, wavenumber))
    {
      if (!(realpoint(i, j, k) || complexpoint(i, j, k)))
      {
        if (i == 0)
        {
          ip = 0;
        }
        else
        {
          ip = NL - i;
        }

        if (j == 0)
        {
          jp = 0;
        }
        else
        {
          jp = NL - j;
        }

        if (k == 0)
        {
          kp = 0;
        }
        else
        {
          kp = NL - k;
        }

        dwk[i][j][k] = conj(dwk[ip][jp][kp]);
        count++;
      }
    }
  }

  if (count != 0)
  {
    LOOP{
      if (innsigma(i,j,k,wavenumber)) {
        dwk[i][j][k] /= sqrt(count);
        dwk[i][j][k] += bias/count;
      }
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

  return (wavenumber-dn/2. <= ntnorm && ntnorm < wavenumber+dn/2.);
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
