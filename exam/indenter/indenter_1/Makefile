
INCLUDE_DIRS := $(PYTHON_INCLUDE) /usr/local/include /~/Downloads/hdf5-1.10.1/hdf5/include/
gsl = /usr/local/Cellar/gsl/1.16
fftw = /home/yang/Downloads/fftw-3.3.6-pl2/
hdf5 = /home/yang/Downloads/hdf5-1.8.18/hdf5/

default: 
	mpic++ -Wall -O2 -c kd_alloc.cc
	mpic++ -Wall -O2 -c parameter_file.cc
	mpic++ -Wall -O2 -c log.cc
	mpic++ -Wall -O2 -c initialize.cc
	mpic++ -Wall -O2 -c -std=c++11 main_eigen.cc -I$(fftw)/include -I$(hdf5)/include
	mpic++ -Wall  kd_alloc.o parameter_file.o log.o initialize.o main_eigen.o -L$(fftw)/lib -L$(hdf5)/lib -lfftw3_mpi -lfftw3 -lhdf5

