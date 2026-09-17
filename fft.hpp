#ifndef INCLUDED_fft_hpp_
#define INCLUDED_fft_hpp_

#include <complex>
#include <vector>
#include <cstring>
#include <fftw3.h>

// Flat (row-major) NL*NL*NL complex grid: index(i,j,k) = (i*NL + j)*NL + k.
// Replacing std::vector<std::vector<std::vector<T>>> with a single flat
// vector removes ~NL*NL+NL small heap allocations and two levels of
// pointer-chasing on every access, which is what made the original
// element-wise loops and the FFT copy-loops so slow and cache-hostile.
using Grid  = std::vector<std::complex<double>>;
using RGrid = std::vector<double>;

// A reusable, thread-parallel FFTW plan for forward 3D complex DFTs of a
// fixed box size. The program performs ~30 transforms of the *same* size
// (NL), so instead of fftw_malloc'ing fresh buffers and building a fresh
// FFTW_ESTIMATE plan every single call (as the original fftw() did), we
// allocate the aligned scratch buffers and build one FFTW_MEASURE plan
// once, then reuse it for every call via fftw_execute_dft's "new-array
// execute" interface (memcpy'ing into/out of the aligned scratch, so the
// plan's alignment assumptions are always satisfied regardless of how the
// caller's std::vector happens to be aligned).
class FFTW3DPlan
{
public:
  explicit FFTW3DPlan(int boxSize, int nthreads) : N(boxSize)
  {
#ifndef NO_FFTW_THREADS
    // fftw_init_threads()/fftw_plan_with_nthreads() live in libfftw3_threads
    // (or libfftw3_omp), not in plain libfftw3. If that library isn't
    // available the Makefile defines NO_FFTW_THREADS and skips these calls
    // entirely so the link still succeeds; each transform then just runs
    // single-threaded inside FFTW, everything else in this program is
    // unaffected.
    static bool threadsInited = false;
    if (!threadsInited)
    {
      fftw_init_threads();
      threadsInited = true;
    }
    fftw_plan_with_nthreads(nthreads > 0 ? nthreads : 1);
#else
    (void)nthreads;
#endif

    const size_t n = static_cast<size_t>(N) * N * N;
    in_  = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * n));
    out_ = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * n));
    // FFTW_MEASURE destroys the contents of in_/out_ while timing candidate
    // algorithms; that is fine because in_/out_ here are just scratch, and
    // the extra one-time planning cost is amortized over ~30 executions.
    plan_ = fftw_plan_dft_3d(N, N, N, in_, out_, FFTW_FORWARD, FFTW_MEASURE);
  }

  ~FFTW3DPlan()
  {
    fftw_destroy_plan(plan_);
    fftw_free(in_);
    fftw_free(out_);
  }

  FFTW3DPlan(const FFTW3DPlan &) = delete;
  FFTW3DPlan &operator=(const FFTW3DPlan &) = delete;

  // Runs the cached plan on bk (flat, row-major N*N*N) and returns the result.
  Grid execute(const Grid &bk) const
  {
    const size_t n = static_cast<size_t>(N) * N * N;
    // std::complex<double> is guaranteed layout-compatible with double[2]
    // (and hence with fftw_complex); the void* casts just keep -Wclass-memaccess
    // quiet about memcpy-ing a "non-trivial" class type, which is a false
    // positive for this specific, standard-sanctioned layout.
    std::memcpy(in_, static_cast<const void *>(bk.data()), sizeof(fftw_complex) * n);
    fftw_execute_dft(plan_, in_, out_);
    Grid bx(n);
    std::memcpy(static_cast<void *>(bx.data()), out_, sizeof(fftw_complex) * n);
    return bx;
  }

private:
  int N;
  fftw_complex *in_;
  fftw_complex *out_;
  fftw_plan plan_;
};

#endif
