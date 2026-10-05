CC = /usr/bin/clang
CXX = /usr/bin/clang++
FC = /opt/homebrew/bin/gfortran
FFLAGS = -fPIC -O3 -std=legacy -w -mno-outline-atomics
FFLAGS_DYN = -shared -fPIC
CFLAGS = -fPIC -O3 -arch arm64 -std=c99
CXXFLAGS = -std=c++11 -DSCTL_PROFILE=-1 -fPIC -O3 -arch arm64
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp
