

#include <fftw3-mpi.h>
#include <stdio.h>

void log_greens_function(double *** G, double ** kxy, ptrdiff_t local_n0, ptrdiff_t N1);
void log_elastic_tensors(double **** lam, double **** epsT);
