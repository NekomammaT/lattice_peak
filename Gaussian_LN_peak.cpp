#include <iostream>
#include <fstream>
#include <random>
#include <sys/time.h>
#include "fft.hpp"
#include "vec_op.hpp"

std::vector<std::vector<std::vector<std::complex<double>>>> dwk(int wavenumber, double bias, int seed);
double powerspectrum(double wavenumber);
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
const int NL = 256; // Box size NL
const int nsigma = 16;
const double s2 = 0.1;
const double dn = 1; // Thickness of nsigma sphere shell
const double bias = 0.3;
const std::string mapfileprefix = std::string("data/LN") //+ std::to_string((int)s2) 
  + std::string("0,1")
  + std::string("_map_");
// const std::string biasedfileprefix = "data/mono_biased_";
const std::string laplacianfileprefix = std::string("data/LN") //+ std::to_string((int)s2) 
  + std::string("0,1")
  + std::string("_laplacian_");
// const std::string powerfileprefix = "data/mono_power_";
const std::string peakfileprefix = std::string("data/LN") //+ std::to_string((int)s2) 
  + std::string("0,1")
  + std::string("_peak_");

// power spectrum
double powerspectrum(int wavenumber)
{
  return exp(-pow(log(wavenumber)-log(nsigma),2)/2/s2) / sqrt(2*M_PI*s2);
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
  std::ofstream mapfile(mapfileprefix + std::to_string(seed) + ".csv");
  // std::ofstream biasedfile(biasedfileprefix + std::to_string(seed) + ".dat");
  std::ofstream laplacianfile(laplacianfileprefix + std::to_string(seed) + ".csv");
  // std::ofstream powerfile(powerfileprefix + std::to_string(seed) + ".dat");
  std::ofstream peakfile(peakfileprefix + std::to_string(seed) + ".csv");

  // ----------- unbiased map -----------
  std::vector<std::vector<std::vector<std::complex<double>>>> gk = dwk(1, 0., seed)*sqrt(powerspectrum(1)*dn);
  
  for (int i = 2; i < sqrt(3)*NL; i++)
  {
    gk += dwk(i, 0., seed)*sqrt(powerspectrum(i)*dn/i);
    std::cout << "\r" << i << " / " << std::floor(sqrt(3)*NL) << std::flush;
  }
  std::cout << std::endl;
  
  std::vector<std::vector<std::vector<std::complex<double>>>> gx = fftw(gk);
  double sigma1sq = exp(2*s2)*pow(2*M_PI*nsigma/NL,2);
  double sigma2sq = exp(2*2*2*s2)*pow(2*M_PI*nsigma/NL,4);
  double sigma4sq = exp(2*4*4*s2)*pow(2*M_PI*nsigma/NL,8);
  
  LOOP
  {
    mapfile << gx[i][j][k].real();
    if (i != NL-1 || j != NL-1 || k != NL-1) mapfile << ','; 
  }
  mapfile << std::endl;

  std::cout << "Exported to " << mapfileprefix + std::to_string(seed) + ".csv" << std::endl;

  // ----------- laplacian -----------
  std::vector<std::vector<std::vector<std::complex<double>>>> D2gk = gk;
  std::vector<std::vector<std::vector<std::complex<double>>>> D2D2gk = gk;
  LOOP
    {
      int nxt = shiftedindex(i);
      int nyt = shiftedindex(j);
      int nzt = shiftedindex(k);
      double ntnorm = sqrt(nxt*nxt+nyt*nyt+nzt*nzt);
      
      D2gk[i][j][k] *= pow(2*M_PI*ntnorm/NL,2);
      D2D2gk[i][j][k] *= pow(2*M_PI*ntnorm/NL,4);
    }
  std::vector<std::vector<std::vector<std::complex<double>>>> D2gx = fftw(D2gk);
  std::vector<std::vector<std::vector<std::complex<double>>>> D2D2gx = fftw(D2D2gk);

  LOOP
  {
    laplacianfile << D2gx[i][j][k].real(); //* sigma1sq/sigma2sq;
    if (i != NL-1 || j != NL-1 || k != NL-1) laplacianfile << ','; 
  }
  laplacianfile << std::endl;

  std::cout << "Exported to " << laplacianfileprefix + std::to_string(seed) + ".csv" << std::endl;

  // ----------- gradient ------------
  std::vector<std::vector<std::vector<double>>> DxD2gx(NL, std::vector<std::vector<double>>(NL, std::vector<double>(NL, 0)));
  std::vector<std::vector<std::vector<double>>> DyD2gx = DxD2gx;
  std::vector<std::vector<std::vector<double>>> DzD2gx = DxD2gx;
  LOOP
    {
      if (i == NL-1) {
	DxD2gx[i][j][k] = D2gx[0][j][k].real() - D2gx[i][j][k].real();
      } else {
	DxD2gx[i][j][k] = D2gx[i+1][j][k].real() - D2gx[i][j][k].real();
      }

      if (j == NL-1) {
	DyD2gx[i][j][k] = D2gx[i][0][k].real() - D2gx[i][j][k].real();
      } else {
	DyD2gx[i][j][k] = D2gx[i][j+1][k].real() - D2gx[i][j][k].real();
      }

      if (k == NL-1) {
	DzD2gx[i][j][k] = D2gx[i][j][0].real() - D2gx[i][j][k].real();
      } else {
	DzD2gx[i][j][k] = D2gx[i][j][k+1].real() - D2gx[i][j][k].real();
      }
    }

  std::vector<std::vector<std::vector<double>>> DxD2gxrot = DxD2gx;
  std::vector<std::vector<std::vector<double>>> DyD2gxrot = DyD2gx;
  std::vector<std::vector<std::vector<double>>> DzD2gxrot = DzD2gx;

  LOOP
    {
      if (i == NL-1) {
	DxD2gxrot[i][j][k] = DxD2gx[0][j][k];
      } else {
	DxD2gxrot[i][j][k] = DxD2gx[i+1][j][k];
      }

      if (j == NL-1) {
	DyD2gxrot[i][j][k] = DyD2gx[i][0][k];
      } else {
	DyD2gxrot[i][j][k] = DyD2gx[i][j+1][k];
      }

      if (k == NL-1) {
	DzD2gxrot[i][j][k] = DzD2gx[i][j][0];
      } else {
	DzD2gxrot[i][j][k] = DzD2gx[i][j][k+1];
      }
    }

  int ip, jp, kp;
  LOOP
    {
      if (DxD2gx[i][j][k] * DxD2gxrot[i][j][k] < 0 && DxD2gx[i][j][k] > 0 && 
	  DyD2gx[i][j][k] * DyD2gxrot[i][j][k] < 0 && DyD2gx[i][j][k] > 0 &&
	  DzD2gx[i][j][k] * DzD2gxrot[i][j][k] < 0 && DzD2gx[i][j][k] > 0) {
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

	peakfile << ip << ',' << jp << ',' << kp << ','
		 << D2gx[ip][jp][kp].real() << ',' // * sigma1sq/sigma2sq << ','
		 << sqrt(D2D2gx[ip][jp][kp].real()/D2gx[ip][jp][kp].real())
			 //* sqrt(sigma2sq/sigma4sq))
		 << std::endl;
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

