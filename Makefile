MODEL = Gaussian_mono

CXX := g++
#CXX := FCCpx

CXXFLAGS := -std=c++11 -Ofast -lfftw3 -lm
#CXXFLAGS := -Nclang -Ofast -fopenmp -lfftw3_threads -lfftw3 -lm

all: $(MODEL)
$(MODEL): $(MODEL).o
	$(CXX) $(CXXFLAGS) -o $(MODEL) $(MODEL).o

$(MODEL).o: fft.hpp

clean:
	$(RM) *.o
	$(RM) $(MODEL)