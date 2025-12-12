#include <iostream>
#include <fstream>
#include <random>
#include <sys/time.h>
#include "fft.hpp"
#include "vec_op.hpp"

std::vector<std::vector<std::vector<std::complex<double>>>> dwk(int wavenumber, std::mt19937& engine);
std::vector<std::vector<std::vector<std::complex<double>>>> Bk(int wavenumber, double bias);
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
const int nbias = 16;
const double s2 = 0.01;
const std::string s2value = //std::to_string((int)s2); 
"0,01";
const double dn = 1; // Thickness of nsigma sphere shell
const double biascoeff = 10; //0.48;
const std::string mukfilename = std::string("data/LN_muk_") + s2value + std::string("_") + std::to_string(NL) + std::string("_") + std::to_string(nsigma) + std::string("_") + std::to_string(nbias) + std::string(".csv");

// power spectrum
double powerspectrum(int wavenumber)
{
  return exp(-pow(log(wavenumber)-log(nsigma),2)/2/s2) / sqrt(2*M_PI*s2);
}


int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    std::cerr << "Specify the seed correctly." << std::endl;
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
  std::cout << "seed = " << seed << std::endl;
  std::mt19937 engine(std::hash<int>{}(seed));
  std::ofstream mukfile(mukfilename, std::ios::app);

  double bias = biascoeff / sqrt(powerspectrum(nbias)*dn/nbias);

  // ----------- unbiased/biased map -----------
  std::vector<std::vector<std::vector<std::complex<double>>>> Wk(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));
  std::vector<std::vector<std::vector<std::complex<double>>>> gkbias = dwk(1, engine)*sqrt(powerspectrum(1)*dn);
  
  for (int i = 2; i < sqrt(3)*NL; i++)
  {
    std::vector<std::vector<std::vector<std::complex<double>>>> dwkt = dwk(i, engine)*(sqrt(powerspectrum(i)*dn/i));
    gkbias += dwkt;

    if (i == nbias) {
      Wk += dwkt;
    }

    std::cout << "\r" << i << " / " << std::floor(sqrt(3)*NL) << std::flush;
  }
  std::cout << std::endl;

  gkbias += Bk(nbias, bias)*sqrt(powerspectrum(nbias)*dn/nbias);
  
  std::vector<std::vector<std::vector<std::complex<double>>>> Dgk = gkbias;
  std::vector<std::vector<std::vector<std::complex<double>>>> DDgk = gkbias;
  LOOP
  {
    int nxt = shiftedindex(i);
    int nyt = shiftedindex(j);
    int nzt = shiftedindex(k);
    double ntnorm = sqrt(nxt*nxt+nyt*nyt+nzt*nzt);

    Dgk[i][j][k] *= pow(2*M_PI*ntnorm/NL,2);
    DDgk[i][j][k] *= pow(2*M_PI*ntnorm/NL,4);
  }
  std::vector<std::vector<std::vector<std::complex<double>>>> Dgx = fftw(Dgk);
  std::vector<std::vector<std::vector<std::complex<double>>>> DDgx = fftw(DDgk);

  std::vector<double> Dgx1d(NL*NL*NL, 0);
  LOOP
  {
    Dgx1d[i*NL*NL+j*NL+k] = Dgx[i][j][k].real();
  }
  auto iter = std::max_element(Dgx1d.begin(), Dgx1d.end());
  size_t index = std::distance(Dgx1d.begin(), iter);

  int imax = index / (NL * NL);
  int jmax = (index - imax * NL * NL) / NL;
  int kmax = index - imax * NL * NL - jmax * NL;
  double mu2 = Dgx[imax][jmax][kmax].real(); 
  double k3 = sqrt(DDgx[imax][jmax][kmax].real() / Dgx[imax][jmax][kmax].real()); 

  std::vector<std::vector<std::vector<std::complex<double>>>> Wx = fftw(Wk);
  double lnw = -bias*Wx[0][0][0].real() - 0.5*bias*bias;

  //mukfile << seed << ',' << mu2 << ',' << k3 << ',' << lnw << std::endl;
  std::cout << mu2 << ',' << imax << ',' << jmax << ',' << kmax << ',' << Wx[0][0][0].real() << ',' << lnw << std::endl;

  // ---------- stop timer ----------
  gettimeofday(&Nv, &Nz);
  after = (double)Nv.tv_sec + (double)Nv.tv_usec * 1.e-6;
  std::cout << after - before << " sec." << std::endl;
  // -------------------------------------

  return 0;
}


// -----------------------------------------------

std::vector<std::vector<std::vector<std::complex<double>>>> dwk(int wavenumber, std::mt19937& engine)
{
  std::vector<std::vector<std::vector<std::complex<double>>>> dwk(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));

  int count = 0;

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
      if (innsigma(i,j,k,wavenumber)) dwk[i][j][k] /= sqrt(count);
    }
  }

  return dwk;
}

std::vector<std::vector<std::vector<std::complex<double>>>> Bk(int wavenumber, double bias)
{
  std::vector<std::vector<std::vector<std::complex<double>>>> Bk(NL, std::vector<std::vector<std::complex<double>>>(NL, std::vector<std::complex<double>>(NL, 0)));

  int count = 0;

  LOOP
  {
    if (innsigma(i, j, k, wavenumber)) count++;
  }

  if (count != 0)
  {
    LOOP{
      if (innsigma(i,j,k,wavenumber)) Bk[i][j][k] = bias/count;
    }
  }

  return Bk;
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
