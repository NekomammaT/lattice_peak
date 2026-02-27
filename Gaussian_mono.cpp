#include <iostream>
#include <fstream>
#include <random>
#include <sys/time.h>
#include "fft.hpp"
#include "vec_op.hpp"

std::vector<std::vector<std::vector<std::complex<double>>>> dwk(int wavenumber, std::mt19937& engine);
std::vector<std::vector<std::vector<std::complex<double>>>> Bk(int wavenumber, double bias);
double WRTH(double z);
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
const double dn = 1; // Thickness of nsigma sphere shell
const double bias = 8; //10;
const double As = 5e-3;
const std::string mukfilename = std::string("data/mono_muk_") + std::to_string(NL) + std::string("_") + std::to_string(nsigma) + std::string("_2.csv");

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

  // ----------- unbiased map -----------
  std::vector<std::vector<std::vector<std::complex<double>>>> gk = dwk(nsigma, engine);
  std::vector<std::vector<std::vector<std::complex<double>>>> gx = fftw(gk);
  
  // ----------- biased map -----------
  std::vector<std::vector<std::vector<std::complex<double>>>> gkbias = gk + Bk(nsigma, bias);
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
  std::vector<std::vector<std::vector<std::complex<double>>>> gxbias = fftw(gkbias);
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
  double lnw = -bias*gx[0][0][0].real() - 0.5*bias*bias;

  double Cmax = 0;
  int rsmax;
  double sigma1 = 2*M_PI*nsigma/NL;
  double sigma2 = pow(2*M_PI*nsigma/NL,2);
  for (int rs = 1; rs <= NL/2; rs++) {
    std::vector<std::vector<std::vector<std::complex<double>>>> rzpk = gkbias;
    LOOP
    {
      int nxt = shiftedindex(i);
      int nyt = shiftedindex(j);
      int nzt = shiftedindex(k);
      double ntnorm = sqrt(nxt*nxt+nyt*nyt+nzt*nzt);
      double kr = 2*M_PI*ntnorm*rs/NL;

      rzpk[i][j][k] *= -kr*kr/3*WRTH(kr)*sqrt(As);
    }
    std::vector<std::vector<std::vector<std::complex<double>>>> rzpx = fftw(rzpk);
    double compaction = 2./3*(1-pow(1+rzpx[imax][jmax][kmax].real(),2));
    if (compaction > Cmax) {
      Cmax = compaction;
      rsmax = rs;
    }
  }
  int count = 0;
  double zetam = 0;
  LOOP
  {
    int nxt = shiftedindex(i);
    int nyt = shiftedindex(j);
    int nzt = shiftedindex(k);
    int nxm = shiftedindex(imax);
    int nym = shiftedindex(jmax);
    int nzm = shiftedindex(kmax);
    double dr = sqrt((nxt-nxm)*(nxt-nxm)+(nyt-nym)*(nyt-nym)+(nzt-nzm)*(nzt-nzm));
    if (fabs(dr-rsmax) < 1./2) {
      zetam += gxbias[i][j][k].real() * sqrt(As);
      count++;
    }
  }
  zetam /= count;

  mukfile << seed << ',' << mu2 << ',' << k3 << ',' << k3*rsmax << ',' << zetam << ',' << Cmax << ',' << lnw << std::endl;
  std::cout << seed << ',' << mu2 << ',' << k3 << ',' << k3*rsmax << ',' << zetam << ',' << Cmax << ',' << lnw << std::endl;

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
