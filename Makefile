CXX = g++
NVCC = nvcc
CXXFLAGS = -std=c++11
NVCCFLAGS = -std=c++11

all: seq.exe kernel.exe

seq.exe: seq.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lstdc++

kernel.exe: kernel.cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

clean:
	rm -f seq.exe kernel.exe
