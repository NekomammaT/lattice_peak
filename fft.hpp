#ifndef INCLUDED_fft_hpp_
#define INCLUDED_fft_hpp_

#include <complex>
#include <vector>
#include <fftw3.h>

std::vector<std::vector<std::vector<std::complex<double>>>> fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk); // FFTW


// FFTW
std::vector<std::vector<std::vector<std::complex<double>>>> fftw(const std::vector<std::vector<std::vector<std::complex<double>>>>& bk) {

  const int BoxSize = bk.size();
  fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * BoxSize * BoxSize * BoxSize);
  fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * BoxSize * BoxSize * BoxSize);

  // bk (3D std::vector) -> in (1D fftw_complex) に変換
  int idx = 0;
  for (int i = 0; i < BoxSize; ++i) {
    for (int j = 0; j < BoxSize; ++j) {
      for (int k = 0; k < BoxSize; ++k) {
        in[idx][0] = bk[i][j][k].real();
        in[idx][1] = bk[i][j][k].imag();
        idx++;
      }
    }
  }

  // FFTWプランを作成
  fftw_plan plan = fftw_plan_dft_3d(BoxSize, BoxSize, BoxSize, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  // FFTを実行
  fftw_execute(plan);

  // out (1D fftw_complex) -> bx (3D std::vector) に変換
  std::vector<std::vector<std::vector<std::complex<double>>>> bx(BoxSize, std::vector<std::vector<std::complex<double>>>(BoxSize, std::vector<std::complex<double>>(BoxSize)));

  idx = 0;
  for (int i = 0; i < BoxSize; ++i) {
    for (int j = 0; j < BoxSize; ++j) {
      for (int k = 0; k < BoxSize; ++k) {
        bx[i][j][k] = std::complex<double>(out[idx][0], out[idx][1]);
        idx++;
      }
    }
  }

  // メモリ解放
  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);

  return bx;
}


#endif
