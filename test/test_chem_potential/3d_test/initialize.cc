
#include "initialize.h"

void initialize(double ** eta, double ** eta_old, ptrdiff_t local_n0, ptrdiff_t N1)
{
    const int N1r = 2*(N1/2+1);
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        eta[0][ndx] = 0.00;
        eta[1][ndx] = 0.00;
        eta[2][ndx] = 0.00;

        eta_old[0][ndx] = eta[0][ndx];
        eta_old[1][ndx] = eta[1][ndx];
        eta_old[2][ndx] = eta[2][ndx];
    }
}

void initialize_phi_0(double * phi, ptrdiff_t local_n0, ptrdiff_t N1)
{
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        phi[ndx] = 0.05;
    }
}

void initialize_phi_1(double * phi, ptrdiff_t local_n0, ptrdiff_t N1)
{
    const int N1r = 2*(N1/2+1);
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        phi[ndx] = 1;
    }
}

void initialize_lsf_stripe(double * lsf, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1)
{
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1 + j;
        int x = local_0_start + i;

        lsf[ndx] = -1;
        if (x > 25 && x < 75) lsf[ndx] = 1;
    }
}

void initialize_lsf_circle(double * lsf, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1)
{
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1 + j;
        int x = local_0_start + i;

        double rx = x - local_n0/2.0;
        double ry = j - N1 / 2.0;
        double rr = sqrt(rx*rx + ry*ry);

        lsf[ndx] = 1;
        if (rr<50) lsf[ndx] = -1;
    }
}

void initialize_lsf_zigzag(double * lsf, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1)
{
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1 + j;
        int x = local_0_start + i;

        double y1 = x;
        double y2 = -x+N1;

        if (j<y1 && j>y2) lsf[ndx] = 1;
        else if (j>y1 && j<y2) lsf[ndx] = 1;
        else lsf[ndx] = -1;

    }
}

void diffuse_lsf(double * lsf, ptrdiff_t local_n0, ptrdiff_t N1)
{
    double width = 5.0;
    double half_width = 0.5*width;

    SDF sdf;
    double band = 20;
    double tolerance = 0.01;
    int dims[2] = {local_n0, N1};
    sdf.construct(lsf, dims, 2, band, tolerance);

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1 + j;
        lsf[ndx] = tanh(lsf[ndx]/half_width);
        lsf[ndx] = 0.5*(lsf[ndx]+1);
    }
}

void copy_lsf(double * lsf, double * phi, ptrdiff_t local_n0, ptrdiff_t N1)
{
    int N1r = 2*(N1/2 + 1);

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx_lsf = i*N1 + j;
        int ndx_phi = i*N1r + j;
      if(lsf[ndx_lsf] <= 0.05)
       phi[ndx_phi] = 0.05;
      else 
       phi[ndx_phi] = lsf[ndx_lsf];
      
    }
}

void initialize_w(double * w, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1)
{   
    int N1r = 2*(N1/2 + 1);
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        int x = local_0_start + i;
        int y = j;

        w[ndx] = 0.48*sin(2*3.1415926*y/20);
    }
}

