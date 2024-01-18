CXX = g++
NVCC = nvcc
CXXFLAGS = -std=c++11
NVCCFLAGS = -std=c++11

all: seq kernel

seq: seq.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lstdc++

kernel: kernel.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -f seq kernel
