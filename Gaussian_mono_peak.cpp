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
// const std::string biasedfileprefix = "data/mono_biased_";
const std::string laplacianfileprefix = "data/mono_laplacian_";
// const std::string powerfileprefix = "data/mono_power_";
const std::string peakfileprefix = "data/mono_peak_";

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
  // std::ofstream biasedfile(biasedfileprefix + std::to_string(seed) + ".dat");
  // std::ofstream laplacianfile(laplacianfileprefix + std::to_string(seed) + ".dat");
  // std::ofstream powerfile(powerfileprefix + std::to_string(seed) + ".dat");
  std::ofstream peakfile(peakfileprefix + std::to_string(seed) + ".csv");

  // ----------- unbiased map -----------
  std::vector<std::vector<std::vector<std::complex<double>>>> gk = dwk(nsigma, 0., seed);
  std::vector<std::vector<std::vector<std::complex<double>>>> gx = fftw(gk);
  
  LOOP
  {
    mapfile << gx[i][j][k].real();
    if (i != NL-1 || j != NL-1 || k != NL-1) mapfile << ','; 
  }
  mapfile << std::endl;

  std::cout << "Exported to " << mapfileprefix + std::to_string(seed) + ".csv" << std::endl;

  // ----------- gradient ------------
  std::vector<std::vector<std::vector<double>>> Dxgx(NL, std::vector<std::vector<double>>(NL, std::vector<double>(NL, 0)));
  std::vector<std::vector<std::vector<double>>> Dygx = Dxgx;
  std::vector<std::vector<std::vector<double>>> Dzgx = Dxgx;
  LOOP
    {
      if (i == NL-1) {
	Dxgx[i][j][k] = gx[0][j][k].real() - gx[i][j][k].real();
      } else {
	Dxgx[i][j][k] = gx[i+1][j][k].real() - gx[i][j][k].real();
      }

      if (j == NL-1) {
	Dygx[i][j][k] = gx[i][0][k].real() - gx[i][j][k].real();
      } else {
	Dygx[i][j][k] = gx[i][j+1][k].real() - gx[i][j][k].real();
      }

      if (k == NL-1) {
	Dzgx[i][j][k] = gx[i][j][0].real() - gx[i][j][k].real();
      } else {
	Dzgx[i][j][k] = gx[i][j][k+1].real() - gx[i][j][k].real();
      }
    }

  std::vector<std::vector<std::vector<double>>> Dxgxrot = Dxgx;
  std::vector<std::vector<std::vector<double>>> Dygxrot = Dygx;
  std::vector<std::vector<std::vector<double>>> Dzgxrot = Dzgx;

  LOOP
    {
      if (i == NL-1) {
	Dxgxrot[i][j][k] = Dxgx[0][j][k];
      } else {
	Dxgxrot[i][j][k] = Dxgx[i+1][j][k];
      }

      if (j == NL-1) {
	Dygxrot[i][j][k] = Dygx[i][0][k];
      } else {
	Dygxrot[i][j][k] = Dygx[i][j+1][k];
      }

      if (k == NL-1) {
	Dzgxrot[i][j][k] = Dzgx[i][j][0];
      } else {
	Dzgxrot[i][j][k] = Dzgx[i][j][k+1];
      }
    }

  int ip, jp, kp;
  LOOP
    {
      if (Dxgx[i][j][k] * Dxgxrot[i][j][k] < 0 && Dxgx[i][j][k] > 0 && 
	  Dygx[i][j][k] * Dygxrot[i][j][k] < 0 && Dygx[i][j][k] > 0 &&
	  Dzgx[i][j][k] * Dzgxrot[i][j][k] < 0 && Dzgx[i][j][k] > 0) {
	if (i==NL-1) {
	  ip = 0;
	} else {
	  ip = i+1;
	}

	if (j==NL-1) {
	  jp = 0;
	} else {
	  jp = j+1;
	}

	if (k==NL-1) {
	  kp = 0;
	} else {
	  kp = k+1;
	}

	peakfile << ip << ',' << jp << ',' << kp << std::endl;
      }
    }
  

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------

  return 0;
}



// -----------------------------------------------

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

