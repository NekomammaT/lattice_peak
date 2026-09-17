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
// fixed box size. The program performs many transforms of the *same* size
// (NL), so instead of fftw_malloc'ing fresh buffers and building a fresh
// FFTW_ESTIMATE plan every single call (as the original fftw() did), we
// allocate ONE aligned scratch buffer and build one FFTW_MEASURE plan once,
// then reuse it for every call.
//
// The transform is done IN PLACE (a single scratch buffer used as both
// FFTW's "in" and "out"), not out-of-place with two buffers. For a
// complex-to-complex transform this is mathematically identical -- FFTW
// fully supports in-place c2c DFTs -- and it halves this class's own fixed
// memory footprint (one NL^3 complex buffer instead of two). At NL=1024
// that's the difference between reserving ~17 GB and ~34 GB for the plan's
// scratch space alone, on top of whatever Grids the caller is holding, so
// it matters a lot for large boxes even though it makes no difference for
// small ones. In-place transforms can occasionally pick a slightly
// different internal algorithm than out-of-place ones; the result is still
// the correct DFT, just possibly not bit-for-bit identical in the last
// rounding digit, which is well within the noise -Ofast already introduces.
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
    buf_ = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * n));
    // FFTW_MEASURE destroys the contents of buf_ while timing candidate
    // algorithms; that is fine because buf_ here is just scratch, and the
    // extra one-time planning cost is amortized over many executions.
    plan_ = fftw_plan_dft_3d(N, N, N, buf_, buf_, FFTW_FORWARD, FFTW_MEASURE);
  }

  ~FFTW3DPlan()
  {
    fftw_destroy_plan(plan_);
    fftw_free(buf_);
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
    std::memcpy(buf_, static_cast<const void *>(bk.data()), sizeof(fftw_complex) * n);
    fftw_execute(plan_); // in-place: buf_ serves as both input and output
    Grid bx(n);
    std::memcpy(static_cast<void *>(bx.data()), buf_, sizeof(fftw_complex) * n);
    return bx;
  }

private:
  int N;
  fftw_complex *buf_;
  fftw_plan plan_;
};

#endif
