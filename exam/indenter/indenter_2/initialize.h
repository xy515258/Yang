
#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <fftw3-mpi.h>
#include <math.h>
#include "sdf.h"

void initialize(double ** eta, double ** eta_old, ptrdiff_t local_n0, ptrdiff_t N1);
void initialize_phi_0(double * phi, ptrdiff_t local_n0, ptrdiff_t N1);
void initialize_phi_1(double * phi, ptrdiff_t local_n0, ptrdiff_t N1);
void initialize_lsf_stripe(double * lsf, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1);
void initialize_lsf_oblique(double * lsf, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1);
void initialize_lsf_circle(double * lsf, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1);
void initialize_lsf_zigzag(double * lsf, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1);
void diffuse_lsf(double * lsf, ptrdiff_t local_n0, ptrdiff_t N1);
void copy_lsf(double * lsf, double * phi, ptrdiff_t local_n0, ptrdiff_t N1);
void initialize_w(double * w, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1);


#endif
