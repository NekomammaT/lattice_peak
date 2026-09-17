MODEL = Gaussian_mono_peak

#CXX := FCCpx

# --- Compiler selection -----------------------------------------------
# On macOS, the `g++`/`gcc` on PATH are normally just Apple Clang under an
# alias (Apple stopped shipping real GCC years ago). Apple Clang does not
# support -fopenmp out of the box, so plain `CXX := g++` silently builds
# with Clang and then fails (or ignores) -fopenmp.
#
# This picks a real GNU g++ if one is installed via Homebrew
# (`brew install gcc`), trying the versioned names Homebrew actually
# installs (g++-15 down to g++-11), and only falls back to whatever `g++`
# resolves to (Apple Clang, most likely) if none of those exist.
GNU_GXX := $(firstword $(wildcard \
	/opt/homebrew/bin/g++-15 /opt/homebrew/bin/g++-14 /opt/homebrew/bin/g++-13 /opt/homebrew/bin/g++-12 /opt/homebrew/bin/g++-11 \
	/usr/local/bin/g++-15 /usr/local/bin/g++-14 /usr/local/bin/g++-13 /usr/local/bin/g++-12 /usr/local/bin/g++-11))

ifdef GNU_GXX
CXX := $(GNU_GXX)
else
CXX := g++
endif

IS_APPLE_CLANG := $(shell $(CXX) --version 2>/dev/null | grep -qi "Apple clang" && echo yes)

ifeq ($(IS_APPLE_CLANG),yes)
# No OpenMP support without extra setup (libomp + -Xpreprocessor). The
# program still builds and runs correctly this way -- it just skips the
# OpenMP-parallelised loops (the #pragma omp lines are simply ignored) and
# falls back to running them single-threaded. It still gets the flat-array
# rewrite, the FFTW threading below, and the algorithmic changes.
CXXFLAGS := -std=c++20 -Ofast
$(info === Building with Apple Clang ($(CXX)): no -fopenmp support, loops will run single-threaded. ===)
$(info === For full multi-core speed: brew install gcc   (this Makefile will then pick it up automatically) ===)
else
CXXFLAGS := -std=c++20 -Ofast -march=native -fopenmp
endif

# --- FFTW threading backend --------------------------------------------
# Not every FFTW install ships the threaded library (Homebrew's normally
# does, but this checks rather than assuming). Falls back to plain -lfftw3
# if -lfftw3_threads isn't linkable, in which case each individual FFT call
# runs single-threaded but everything else (memory layout, sliding window,
# cached plan, fast CSV writer, and any OpenMP loops above) still applies.
HAVE_FFTW_THREADS := $(shell echo 'int main(){return 0;}' > /tmp/_fftwcheck.c && \
	$(CXX) /tmp/_fftwcheck.c -lfftw3_threads -lfftw3 -lpthread -o /tmp/_fftwcheck 2>/dev/null && echo yes; \
	rm -f /tmp/_fftwcheck.c /tmp/_fftwcheck)

ifeq ($(HAVE_FFTW_THREADS),yes)
LDLIBS := -lfftw3_threads -lfftw3 -lpthread -lm
else
LDLIBS := -lfftw3 -lm
CXXFLAGS += -DNO_FFTW_THREADS
$(info === libfftw3_threads not found: building without FFTW-level threading. ===)
endif

# EXTRA_DEFS lets you switch on extra memory-saving behaviour without
# touching this file, e.g. for a large box:
#   make clean && make EXTRA_DEFS=-DNO_FULL_FIELD_OUTPUT
# -DNO_FULL_FIELD_OUTPUT skips computing/writing mapfile and skips writing
# the laplacianfile CSV (Lpeakfile/Cpeakfile, the peak catalogues, are
# unaffected and still written). See the comment block at the top of
# Gaussian_mono_peak.cpp for exactly what this saves.
CXXFLAGS += $(EXTRA_DEFS)

all: $(MODEL)
$(MODEL): $(MODEL).o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(LDLIBS)

$(MODEL).o: fft.hpp

clean:
	$(RM) *.o
	$(RM) $(MODEL)

# Tip: if you ever run several seeds of this program at once (the way the
# Gaussian_LN_GB*.sh scripts do for a different binary), cap the thread
# count per process so they don't all fight over every core, e.g.:
#   OMP_NUM_THREADS=4 ./Gaussian_mono_peak $i &
# Gaussian_mono_peak.sh itself runs seeds one at a time, so by default this
# program will use all available cores for each seed in turn.
