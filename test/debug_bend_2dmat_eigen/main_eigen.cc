
#include <stdio.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cstring>

#include <fftw3-mpi.h>

#include "parameter_file.h"
#include "kd_alloc.h"
#include "h5_file.h"
#include "log.h"
#include "initialize.h"

const int Re = 0;
const int Im = 1;

struct input_parameters {

    int Nx, Ny;
    int nsteps;
    int out_freq;

    double dx, dt;
    double epsx;
    double epsy;
    double beta;
    double gamma;
    double gamma_eps;
    double alpha;
    double alpha_eps;
    double kappa_el0;
    double kappa_el1;
    double L_eps00;
    double change_etap_thresh;
    double change_eps00_thresh;
    double mu_el0;
    double nu_el0;
    double mu_el1;
    double nu_el1;

    double M0_chem_a;
    double M0_chem_b;
    double M0_chem_c;

    double M1_chem_a;
    double M1_chem_b;
    double M1_chem_c;

    double M0_2H_a;
    double M0_2H_b;
    double M0_Tp_a;
    double M0_Tp_b;

    double M1_2H_a;
    double M1_2H_b;
    double M1_Tp_a;
    double M1_Tp_b;

    double M0_norm;
    double M1_norm;
};

fftw_plan planF_eta[3];
fftw_plan planB_lap[3];
fftw_plan planF_s0n2[3][3];
fftw_plan planF_sig00[2][2];
fftw_plan planB_ft_eps00[2][2];

fftw_plan planB_ux;
fftw_plan planB_uy;

fftw_plan plan_strain_xx;
fftw_plan plan_strain_yy;
fftw_plan plan_strain_xy;

fftw_plan planF_bend_1;
fftw_plan planF_bend_2;



void calc_greens_function(double *** G, double ** kxy, ptrdiff_t local_n0, ptrdiff_t local_0_start, ptrdiff_t N1, struct input_parameters ip)
{
    double pi = 3.14159265359;

    int Nx = ip.Nx;
    int Ny = ip.Ny;
    double dx = ip.dx;
    double mu = ip.mu_el0;
    double nu = ip.nu_el0;
    double radius = 4*pi/(Nx*dx)/(Ny*dx);

    double * kx = new double [local_n0];
    double * ky = new double [N1/2+1];

    double Lx = Nx*dx;
    double Ly = Ny*dx;

    for (int i=0; i<local_n0; i++)
    {
        int local_x = local_0_start + i;
        if ( local_x < Nx/2+1 ) kx[i] = local_x * 2*pi/Lx;
        else kx[i] = (local_x - Nx) * 2*pi/Lx;
    }

    for (int j=0; j<N1/2+1; j++)
        ky[j] = j * 2*pi/Ly;

    for (int i=0; i<local_n0; i++) 
    for (int j=0; j<N1/2+1; j++)
    {
        int ndx = i*(N1/2+1) + j;

        double k2 = kx[i]*kx[i] + ky[j]*ky[j];
        double norm = sqrt(k2);
        double n0, n1;

        if (k2 < radius) {
            n0 = 0;
            n1 = 0;
        } else {
            n0 = kx[i]/norm;
            n1 = ky[j]/norm;
        }

        G[0][0][ndx] = 1.0/mu - (1+nu)*n0*n0/(2*mu);
        G[1][1][ndx] = 1.0/mu - (1+nu)*n1*n1/(2*mu);

        G[0][1][ndx] = -(1+nu)*n0*n1/(2.0*mu);
        G[1][0][ndx] = -(1+nu)*n1*n0/(2.0*mu);

        if (k2 >= radius) {
            G[0][0][ndx] /= k2;
            G[0][1][ndx] /= k2;
            G[1][0][ndx] /= k2;
            G[1][1][ndx] /= k2;
        } else {
            G[0][0][ndx] = 0;
            G[0][1][ndx] = 0;
            G[1][0][ndx] = 0;
            G[1][1][ndx] = 0;
        }

        kxy[0][ndx] = kx[i];
        kxy[1][ndx] = ky[j];
    }

    delete [] kx;
    delete [] ky;
}

void calc_transformation_strains(double **** epsT, struct input_parameters ip)
{
    const int M0 = 0;
    const int M1 = 1;
    const int X = 0;
    const int Y = 1;
    const double Pi = 3.14159265359;

    double e0[2][2][2];
    double thetav[2][3];

    e0[M0][X][X] = (ip.M0_Tp_a - ip.M0_2H_a)/(ip.M0_2H_a*ip.M0_norm*ip.M0_norm);
    e0[M0][Y][Y] = (ip.M0_Tp_b - ip.M0_2H_b)/(ip.M0_2H_b*ip.M0_norm*ip.M0_norm);
    e0[M0][X][Y] = 0.0;
    e0[M0][Y][X] = 0.0;

    e0[M1][X][X] = (ip.M1_Tp_a - ip.M1_2H_a)/(ip.M1_2H_a*ip.M1_norm*ip.M1_norm);
    e0[M1][Y][Y] = (ip.M1_Tp_b - ip.M1_2H_b)/(ip.M1_2H_b*ip.M1_norm*ip.M1_norm);
    e0[M1][X][Y] = 0.0;
    e0[M1][Y][X] = 0.0;

    thetav[M0][0] = 0.0;
    thetav[M0][1] = 2.0*Pi/3.0;
    thetav[M0][2] = -2.0*Pi/3.0;;

    thetav[M1][0] = thetav[M0][0];
    thetav[M1][1] = thetav[M0][1];
    thetav[M1][2] = thetav[M0][2];

    for (int mat=M0; mat<=M1; mat++)
    for (int pp=0; pp<3; pp++)
    {
        double r[2][2];

        r[0][0] =  cos(thetav[mat][pp]);
        r[1][1] =  cos(thetav[mat][pp]);
        r[0][1] = -sin(thetav[mat][pp]);
        r[1][0] =  sin(thetav[mat][pp]);

        for (int mm=0; mm<2; mm++)
        for (int nn=0; nn<2; nn++)
        {
            epsT[mat][pp][mm][nn] = 0;

            for (int ii=0; ii<2; ii++)
            for (int jj=0; jj<2; jj++)
                epsT[mat][pp][mm][nn] += r[mm][ii]*r[nn][jj]*e0[mat][ii][jj];
        }
    }
}

void calc_elastic_tensors(double **** lam, double **** eps0, double **** sig0, double *** sigeps, double mu1, double nu1, ptrdiff_t local_n0, ptrdiff_t N1, double ***** delta_S, double mu2, double nu2, double * phi)
{
    const int N1r = 2*(N1/2+1);

    double E1 = 2*mu1*(1+nu1);
    lam[0][0][0][0] = E1/(1-nu1*nu1);
    lam[0][0][1][1] = E1*nu1/(1-nu1*nu1);
    lam[0][0][1][0] = 0.0;
    lam[0][0][0][1] = 0.0;
    lam[1][1][0][0] = E1*nu1/(1-nu1*nu1);
    lam[1][1][1][1] = E1/(1-nu1*nu1);
    lam[1][1][0][1] = 0;
    lam[1][1][1][0] = 0;
    lam[0][1][0][0] = 0;
    lam[0][1][1][1] = 0;
    lam[0][1][0][1] = mu1;
    lam[0][1][1][0] = mu1;
    lam[1][0][0][0] = 0;
    lam[1][0][1][1] = 0;
    lam[1][0][0][1] = mu1;
    lam[1][0][1][0] = mu1;

    double E2 = 2*mu2*(1+nu2);
    double v1s = 1-nu1*nu1;
    double v2s = 1-nu2*nu2;

    for (int i=0; i<local_n0; i++)
        for (int j=0; j<N1; j++)
        {
            int ndx = i*N1r + j;
     
      delta_S[0][0][0][0][ndx] = ((E1/v1s-E2/v2s)/((E1/v1s-E2/v2s)*(E1/v1s-E2/v2s)-(E1*nu1/v1s-E2*nu2/v2s)*(E1*nu1/v1s-E2*nu2/v2s)))/phi[ndx];
      delta_S[0][0][1][1][ndx] = (-1)*((E1*nu1/v1s-E2*nu2/v2s)/((E1/v1s-E2/v2s)*(E1/v1s-E2/v2s)-(E1*nu1/v1s-E2*nu2/v2s)*(E1*nu1/v1s-E2*nu2/v2s)))/phi[ndx];
      delta_S[0][0][1][0][ndx] = 0.0;
      delta_S[0][0][0][1][ndx] = 0.0;
      delta_S[1][1][0][0][ndx] = (-1)*((E1*nu1/v1s-E2*nu2/v2s)/((E1/v1s-E2/v2s)*(E1/v1s-E2/v2s)-(E1*nu1/v1s-E2*nu2/v2s)*(E1*nu1/v1s-E2*nu2/v2s)))/phi[ndx];
      delta_S[1][1][1][1][ndx] = ((E1/v1s-E2/v2s)/((E1/v1s-E2/v2s)*(E1/v1s-E2/v2s)-(E1*nu1/v1s-E2*nu2/v2s)*(E1*nu1/v1s-E2*nu2/v2s)))/phi[ndx];
      delta_S[1][1][0][1][ndx] = 0;
      delta_S[1][1][1][0][ndx] = 0;
      delta_S[0][1][0][0][ndx] = 0;
      delta_S[0][1][1][1][ndx] = 0;
      delta_S[0][1][0][1][ndx] = 1/(E1/(2*(1+nu1))-E2/(2*(1+nu2)))/phi[ndx];
      delta_S[0][1][1][0][ndx] = 1/(E1/(2*(1+nu1))-E2/(2*(1+nu2)))/phi[ndx];
      delta_S[1][0][0][0][ndx] = 0;
      delta_S[1][0][1][1][ndx] = 0;
      delta_S[1][0][0][1][ndx] = 1/(E1/(2*(1+nu1))-E2/(2*(1+nu2)))/phi[ndx];
      delta_S[1][0][1][0][ndx] = 1/(E1/(2*(1+nu1))-E2/(2*(1+nu2)))/phi[ndx];
}

    for (int pp=0; pp<3; pp++)
    {
        for (int i=0; i<local_n0; i++)
        for (int j=0; j<N1; j++)
        {
            int ndx = i*N1r + j;

            for (int ii=0; ii<2; ii++)
            for (int jj=0; jj<2; jj++)
            {
                sig0[pp][ii][jj][ndx] = 0.0;

                for (int kk=0; kk<2; kk++)
                for (int ll=0; ll<2; ll++)
                    sig0[pp][ii][jj][ndx] += lam[ii][jj][kk][ll] * eps0[pp][kk][ll][ndx];
            }
        }
    }

    for (int p=0; p<3; p++)
    for (int q=0; q<3; q++)
    {

        for (int i=0; i<local_n0; i++)
        for (int j=0; j<N1; j++)
        {
            int ndx = i*N1r + j;
            sigeps[p][q][ndx] = 0;

            for (int ii=0; ii<2; ii++)
            for (int jj=0; jj<2; jj++)
                sigeps[p][q][ndx] += sig0[p][ii][jj][ndx]*eps0[q][ii][jj][ndx];
        }
    }
}

void normalize(double * array, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0)
{
    const int N1r = 2*(N1/2+1);
    const double area = (double) (N0*N1);
    for (ptrdiff_t i=0; i<local_n0; i++)
    for (ptrdiff_t j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        array[ndx] /= area;
    }
}


void introduce_noise(double ** eta, ptrdiff_t local_n0, ptrdiff_t N1)
{
    const int N1r = 2*(N1/2+1);
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        double r;
        r = 2*((float)rand())/((float)RAND_MAX) - 1;
        eta[0][ndx] += 0.0003*r;
        r = 2*((float)rand())/((float)RAND_MAX) - 1;
        eta[1][ndx] += 0.0003*r;
        r = 2*((float)rand())/((float)RAND_MAX) - 1;
        eta[2][ndx] += 0.0003*r;
    }
}

void introduce_noise_eps00(double *** eps00, ptrdiff_t local_n0, ptrdiff_t N1)
{
    const int N1r = 2*(N1/2+1);
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        double r;
        r = 2*((float)rand())/((float)RAND_MAX) - 1;
        eps00[0][0][ndx] += 0.000003*r;
        r = 2*((float)rand())/((float)RAND_MAX) - 1;
        eps00[0][1][ndx] += 0.000003*r;
        r = 2*((float)rand())/((float)RAND_MAX) - 1;
        eps00[1][0][ndx] += 0.000003*r;
        r = 2*((float)rand())/((float)RAND_MAX) - 1;
        eps00[1][1][ndx] += 0.000003*r;
    }
}

void calc_chemical_potential(double *** sig00, double **** eps0, double ** chem, double ** eta, double *** eps00, double ** epsbar, double **** sig0, double *** eps, double ** lap, double * phi, double ** dw, ptrdiff_t local_n0, ptrdiff_t N1, struct input_parameters ip)
{

    double EelAppl[3] = {0, 0, 0};
    const int N1r = 2*(N1/2+1);

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;

        double f_bulk[3]    = {0,0,0};
        double f_squeeze[3] = {0,0,0};
        double f_homo[3]    = {0,0,0};
        double f_hetero[3]  = {0,0,0};

        double eta_sum = eta[0][ndx]*eta[0][ndx] + eta[1][ndx]*eta[1][ndx] + eta[2][ndx]*eta[2][ndx];


        // bulk free energy
        for (int p=0; p<3; p++)
        {
            double eta_sq = eta[p][ndx]*eta[p][ndx];
            double a = (1-phi[ndx])*ip.M0_chem_a + phi[ndx]*ip.M1_chem_a;
            double b = (1-phi[ndx])*ip.M0_chem_b + phi[ndx]*ip.M1_chem_b;
            double c = (1-phi[ndx])*ip.M0_chem_c + phi[ndx]*ip.M1_chem_c;

            f_bulk[p] = eta[p][ndx]*(a - b*eta_sq + c*eta_sum*eta_sum);
        }


        // stress-free strain (squeeze) part of free energy (double sum)
        for (int p=0; p<3; p++)
        for (int kk=0; kk<2; kk++)
        for (int ll=0; ll<2; ll++)
            f_squeeze[p] += 2*eps0[p][kk][ll][ndx]*sig00[kk][ll][ndx]*eta[p][ndx];

   /*     for (int p=0; p<3; p++)
        for (int q=0; q<3; q++)
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
        f_squeeze[p] += 2*sig0[p][ii][jj][ndx]*eps0[q][ii][jj][ndx]*eta[p][ndx]*eta[q][ndx]*eta[q][ndx]; */

        // homogenous, macroscropic strain part of free energy
        for (int p=0; p<3; p++)
        {
            EelAppl[p] = -2*( sig0[p][0][0][ndx]*epsbar[0][0]
                            + sig0[p][1][1][ndx]*epsbar[1][1]
                            + sig0[p][0][1][ndx]*epsbar[0][1]
                            + sig0[p][1][0][ndx]*epsbar[1][0] );
            f_homo[p] = EelAppl[p]*eta[p][ndx];
        }

        // heterogenous, local strain part of free energy
        for (int p=0; p<3; p++)
            f_hetero[p] = -2*eta[p][ndx]* ( sig0[p][0][0][ndx]*(eps[0][0][ndx] + dw[0][ndx]*dw[0][ndx])
                                          + sig0[p][0][1][ndx]*(eps[0][1][ndx] + dw[0][ndx]*dw[1][ndx])
                                          + sig0[p][1][0][ndx]*(eps[1][0][ndx] + dw[1][ndx]*dw[0][ndx])
                                          + sig0[p][1][1][ndx]*(eps[1][1][ndx] + dw[1][ndx]*dw[1][ndx]) );

        for (int p=0; p<3; p++)
        {
            chem[p][ndx]  = f_bulk[p]; 
            chem[p][ndx] -= ip.beta*lap[p][ndx];
            chem[p][ndx] += f_squeeze[p] + f_homo[p] + f_hetero[p];
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////
void calc_uxy(double * ux, double * uy, fftw_complex ** ku, 
              double *** G, double ** kxy, fftw_complex *** ksig00, 
              ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0)
// calculate the displacements in k-space and then inverse fourier tranform to real-space
// F{u} = G*k*F{sig0*eta^2}
//////////////////////////////////////////////////////////////////////////////////////////////
{
    const int N1c = N1/2+1;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1c; j++)
    {
        int ndx = i*N1c + j;

        for (int ii=0; ii<2; ii++)
        {
            ku[ii][ndx][Re] = 0;
            ku[ii][ndx][Im] = 0;

            for (int jj=0; jj<2; jj++)
            for (int kk=0; kk<2; kk++)
            {
                ku[ii][ndx][Re] += G[ii][jj][ndx]*kxy[kk][ndx]*ksig00[jj][kk][ndx][Im];
                ku[ii][ndx][Im] -= G[ii][jj][ndx]*kxy[kk][ndx]*ksig00[jj][kk][ndx][Re];
            }
        }
    }

    // ku -> (ux, uy)
    fftw_execute(planB_ux);
    fftw_execute(planB_uy);

    // nomalize ux,uy - necessary after fftw
    normalize(ux, N0, N1, local_n0);
    normalize(uy, N0, N1, local_n0);
}

void calc_uxy_bending(double * ux, double * uy, fftw_complex ** ku, double ** dw, double ** ddw,
                      double **** lam, double *** G, double ** kxy, fftw_complex *** ksig00,
                      ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0)
{
    const int N1c = N1/2+1;
    const int N1r = 2*(N1/2+1);

    double * N_klm = fftw_alloc_real(local_n0*N1r);
    fftw_complex * kN_klm = fftw_alloc_complex(local_n0*N1c);
    fftw_plan planF_N = fftw_mpi_plan_dft_r2c_2d(N0, N1, N_klm, kN_klm, MPI_COMM_WORLD, FFTW_ESTIMATE);

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1c; j++)
    {
        int ndx = i*N1c + j;

        for (int ii=0; ii<2; ii++)
        {
            ku[ii][ndx][Re] = 0;
            ku[ii][ndx][Im] = 0;

            for (int jj=0; jj<2; jj++)
            for (int kk=0; kk<2; kk++)
            {
                ku[ii][ndx][Re] += G[ii][jj][ndx]*kxy[kk][ndx]*ksig00[jj][kk][ndx][Im];
                ku[ii][ndx][Im] -= G[ii][jj][ndx]*kxy[kk][ndx]*ksig00[jj][kk][ndx][Re];
            }
        }
    }

    for (int kk=0; kk<2; kk++)
    for (int ll=0; ll<2; ll++)
    for (int mm=0; mm<2; mm++)
    {
        for (int i=0; i<local_n0; i++)
        for (int j=0; j<N1r; j++)
        {
            int ndx = i*N1r + j;
            N_klm[ndx] = dw[kk][ndx] * ddw[ll+mm][ndx];
        }

        fftw_execute(planF_N);

        for (int i=0; i<local_n0; i++)
        for (int j=0; j<N1c; j++)
        {
            int ndx = i*N1c + j;

            for (int ii=0; ii<2; ii++)
            for (int jj=0; jj<2; jj++)
            {
                ku[ii][ndx][Re] += G[ii][jj][ndx]*lam[jj][kk][ll][mm]*kN_klm[ndx][Re];
                ku[ii][ndx][Im] += G[ii][jj][ndx]*lam[jj][kk][ll][mm]*kN_klm[ndx][Im];
            }
        }
    }

    fftw_execute(planB_ux);
    fftw_execute(planB_uy);

    normalize(ux, N0, N1, local_n0);
    normalize(uy, N0, N1, local_n0);
}

double update_eta(double ** eta, double ** eta_old, double ** eta_new, double ** chem, ptrdiff_t local_n0, ptrdiff_t N1, struct input_parameters ip)
{
    const int N1r = 2*(N1/2+1);

    double dtg = 0.5*ip.dt*ip.gamma;
    double dtg2 = 1.0/(1.0+dtg);
    double dta2 = ip.dt*ip.dt*ip.alpha*ip.alpha;
    double change_etap_max = 0;

    for (int p=0; p<3; p++)
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;

        eta_new[p][ndx] = dtg2*(2*eta[p][ndx] + (dtg-1)*eta_old[p][ndx] - dta2*chem[p][ndx]);

        double delta = fabs( eta_new[p][ndx] - eta[p][ndx] );
        change_etap_max = std::max(change_etap_max, delta);

        eta_old[p][ndx] = eta[p][ndx];
        eta[p][ndx] = eta_new[p][ndx];
    }

    return change_etap_max;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void calc_ks0n2(double *** s0n2, double **** sig0, double ** eta, ptrdiff_t local_n0, ptrdiff_t N1)
// calculate and transform the nonlinear terms in the displacement equation 
// s0n2 = sig0_{jk}(p,r) * eta^2(p)
// sig0 = the tranformation stresses - lambda * eps0
// eta  = orientation order parameters
// local_n0 = size of local process in x-direction
// N1 = size of local process in y-direction
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
    const int X = 0;
    const int Y = 1;

    const int XX = 0;
    const int XY = 1;
    const int YY = 2;

    const int N1r = 2*(N1/2+1);

    for (int p=0; p<3; p++)
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        double eta_sq = eta[p][ndx] * eta[p][ndx];
        s0n2[p][XX][ndx] = sig0[p][X][X][ndx] * eta_sq;
        s0n2[p][XY][ndx] = sig0[p][X][Y][ndx] * eta_sq;
        s0n2[p][YY][ndx] = sig0[p][Y][Y][ndx] * eta_sq;
    }

    // s0n2 -> ks0n2
    for (int p=0; p<3; p++)
    for (int i=0; i<3; i++)
        fftw_execute(planF_s0n2[p][i]);
}

////////////////////////////////////////////////////////////////////////////////////////////
void calc_eps(double *** eps, fftw_complex ** keps, 
              double ** kxy, fftw_complex ** ku, 
              ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0)
// calculate the heterogeneous strain (delta-epsilon) in k-space and inverse tranform
////////////////////////////////////////////////////////////////////////////////////////////
{
    const int X = 0;
    const int Y = 1;

    const int XX = 0;
    const int YY = 1;
    const int XY = 2;

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<(N1/2+1); j++)
    {
        int ndx = i*(N1/2+1) + j;
        keps[XX][ndx][Re] = -kxy[X][ndx]*ku[X][ndx][Im];
        keps[XX][ndx][Im] =  kxy[X][ndx]*ku[X][ndx][Re];

        keps[YY][ndx][Re] = -kxy[Y][ndx]*ku[Y][ndx][Im];
        keps[YY][ndx][Im] =  kxy[Y][ndx]*ku[Y][ndx][Re];

        keps[XY][ndx][Re] = -0.5*(kxy[Y][ndx]*ku[X][ndx][Im] + kxy[X][ndx]*ku[Y][ndx][Im]);
        keps[XY][ndx][Im] =  0.5*(kxy[Y][ndx]*ku[X][ndx][Re] + kxy[X][ndx]*ku[Y][ndx][Re]);
    }

    // keps -> eps
    fftw_execute(plan_strain_xx);
    fftw_execute(plan_strain_yy);
    fftw_execute(plan_strain_xy);

    normalize(eps[0][0], N0, N1, local_n0);
    normalize(eps[1][1], N0, N1, local_n0);
    normalize(eps[0][1], N0, N1, local_n0);
    std::memcpy(eps[1][0], eps[0][1], sizeof(double)*local_n0*2*(N1/2+1));

    
}

////////////////////////////////////////////////////////////////////////////////////////////
void calc_lap(double ** lap, fftw_complex ** klap, fftw_complex ** keta, double ** kxy, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0)
// calculate the laplacian of the eta order paramters in k-space and inverse transform
// the laplacian comes from the gradient squared energy term
// it will be used to calculate the eta parameter chemical potential
////////////////////////////////////////////////////////////////////////////////////////////
{
    const int X = 0;
    const int Y = 1;

    for (int p=0; p<3; p++)
    {
        // eta -> keta
        fftw_execute(planF_eta[p]);

        for (int i=0; i<local_n0; i++)
        for (int j=0; j<(N1/2+1); j++)
        {
            int ndx = i*(N1/2+1) + j;
            double k2 = kxy[X][ndx]*kxy[X][ndx] + kxy[Y][ndx]*kxy[Y][ndx];
            //klap[p][ndx][Re] = -k2 * keta[p][ndx][Re];
            //klap[p][ndx][Im] = -k2 * keta[p][ndx][Im];

            double rk = (k2 >= 0.0) ? sqrt(k2) : 0.0;
            double kmod = 2.0*(1.0-cos(rk));
            klap[p][ndx][Re] = -kmod * keta[p][ndx][Re];
            klap[p][ndx][Im] = -kmod * keta[p][ndx][Im];
        }

        // klap -> lap
        fftw_execute(planB_lap[p]);
        normalize(lap[p], N0, N1, local_n0);
    }
}


void output(std::string path, double * data, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0)
{
    int np, rank;
    double * buffer;
    int alloc_local = local_n0 * (N1/2+1);
    int tag = 0;
    MPI_Status status;
    int dims[2] = {N0, 2*(N1/2+1)};

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if ( rank == 0 ) {
        buffer = new double [N0*2*(N1/2+1)];
        memcpy(buffer, data, 2*alloc_local*sizeof(double));

        for (int i=1; i<np; i++)
            MPI_Recv(buffer + i*2*alloc_local, 2*alloc_local, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);

        H5File h5;
        h5.open("out.h5", "a");
        h5.write_dataset(path, buffer, dims, 2);
        h5.close();

        delete [] buffer;
    } else {
        MPI_Send(data, 2*alloc_local, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }


}

std::string zeroFill(int x)
{
    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << x;
    return ss.str();
}

////////////////////////////////////////////////////////////////////////////////////////////
void interpolate(double * data, double m0, double m1, double * phi, 
                 ptrdiff_t local_n0, ptrdiff_t N1)
// iterpolate between values for the heterogenous composition
////////////////////////////////////////////////////////////////////////////////////////////
{
    const int N1r = 2*(N1/2 + 1);

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        data[ndx] = (1-phi[ndx])*m0 + phi[ndx]*m1;
    }
}

double calc_area(double ** eta, ptrdiff_t local_n0, ptrdiff_t N0, ptrdiff_t N1, double norm)
{
    const int N1r = 2*(N1/2+1);
    double sum = 0;
    double threshold = 0.5*norm;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        sum += std::abs(eta[0][ndx]) > threshold ? 1 : 0;
        sum += std::abs(eta[1][ndx]) > threshold ? 1 : 0;
        sum += std::abs(eta[2][ndx]) > threshold ? 1 : 0;
    }

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum/(N0*N1);
}

/*
///////////////////////////////////////////////////////////////////////////////////////////////
void calc_dw(double ** dw, double * w, ptrdiff_t local_n0, ptrdiff_t N1, double dx)
// calculate the derivatives of the out-of-plane displacement
// dw[i][ndx] is the first derivative of the out-of-plane displacement (w) in the direction "i = (X,Y)" at grid location "ndx"
// local_n0 is the size of the local process in the x direction
// N1 is the size of the local process in the y direction
// dx is the grid spacing
///////////////////////////////////////////////////////////////////////////////////////////////
{
    const int N1r = 2*(N1/2+1);
    const int X = 0;
    const int Y = 1;

    // first communicate edge data for parallel processes using MPI
    double * left   = new double [(int)N1];
    double * right  = new double [(int)N1];
    MPI_Request request;
    MPI_Status status;
    int dest, source, count, offset, tag=0;
    int rank, np;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // send right, recieve left
    dest = rank+1 ? rank<np-1 : 0;
    source = rank-1 ? rank>0 : np-1;
    offset = (local_n0-1)*N1r;
    count = (int) N1;
    MPI_Isend(w+offset, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &request);     
    MPI_Recv(left, count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status); 

    // send left, recieve right
    dest = rank-1 ? rank>0 : np-1;
    source = rank+1 ? rank<np-1 : 0;
    MPI_Isend(w, count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &request);
    MPI_Recv(right, count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);

    // loop through the local grid and calculate first derivatives

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;

        int xm = (i-1)*N1r + j  ? i>0           : left[j];
        int xp = (i+1)*N1r + j  ? i+1<local_n0  : right[j];
        int ym = ndx - 1        ? j>0           : i*N1r + (N1-1);
        int yp = ndx + 1        ? j+1<N1        : i*N1r + (0);

        dw[X][ndx] = (w[xp] - w[xm])/(2*dx);
        dw[Y][ndx] = (w[yp] - w[ym])/(2*dx);
    }

    delete [] left;
    delete [] right;
}
*/

double max_w ( double * data, int local_n0, int N1 )
{
    double m = 0;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        if (fabs(data[ndx]) > m) m = fabs(data[ndx]);
    }
    return m;
}

double max_eta ( double ** data, int local_n0, int N1)
{
      double m = 0;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<3; ii++)
        if (fabs(data[ii][ndx]) > m) m = fabs(data[ii][ndx]);
    }
    return m;
}


double max_kbend ( fftw_complex * data, int local_n0, int N1)
{
      double m = 0;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1/2+1; j++)
    {
        int ndx = i*(N1/2+1) + j;
        if (fabs(data[ndx][Re]) > m) m = fabs(data[ndx][Re]);
        if (fabs(data[ndx][Im]) > m) m = fabs(data[ndx][Im]);
     }
    return m;
}

double max_ksig00 ( fftw_complex *** data, int local_n0, int N1)
{
      double m = 0;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1/2+1; j++)
    {
        int ndx = i*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
       { if (fabs(data[ii][jj][ndx][Re]) > m) m = fabs(data[ii][jj][ndx][Re]);
         if (fabs(data[ii][jj][ndx][Im]) > m) m = fabs(data[ii][jj][ndx][Im]);
     }}
    return m;
}

double max_kxy ( double ** data, int local_n0, int N1)
{
      double m = 0;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1/2+1; j++)
    {
        int ndx = i*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
       { if (fabs(data[ii][ndx]) > m) m = fabs(data[ii][ndx]);
         if (fabs(data[ii][ndx]) > m) m = fabs(data[ii][ndx]);
     }}
    return m;
}

double max_G ( double *** data, int local_n0, int N1)
{
      double m = 0;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1/2+1; j++)
    {
        int ndx = i*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
       { if (fabs(data[ii][jj][ndx]) > m) m = fabs(data[ii][jj][ndx]);
         if (fabs(data[ii][jj][ndx]) > m) m = fabs(data[ii][jj][ndx]);
     }}
    return m;
}

double avg_eta ( double ** data, int local_n0, int N1)
{
    double n = 0;    

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<3; ii++)
        
         n += data[ii][ndx];  
}
         n = n/local_n0/N1;
    
    return n;
}

double max_eps00 ( double *** data, int local_n0, int N1 )
{
    double m = 0;
    double n = 0;
    double k = 0;

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
        {if (fabs(data[ii][jj][ndx]) > m) m = fabs(data[ii][jj][ndx]);
         n += data[ii][jj][ndx];
         k += fabs(data[ii][jj][ndx]);
}
    }
         n = n/local_n0;
         k = k/local_n0;
    
    return m;
}

double avg_eps00 ( double *** data, int local_n0, int N1 )
{
    double n = 0;    

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
        
         n += data[ii][jj][ndx];  
}
         n = n/local_n0/N1;
    
    return n;
}

double avg_abs_eps00 ( double *** data, int local_n0, int N1 )
{
    double n = 0;    

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
        
         n += fabs(data[ii][jj][ndx]);  
}
         n = n/local_n0/N1;
    
    return n;
}

double max_delta_S ( double ***** data, int local_n0, int N1)
{
    double m = 0;
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
        for (int kk=0; kk<2; kk++)
        for (int ll=0; ll<2; ll++)
        if (fabs(data[ii][jj][kk][ll][ndx]) > m) m = fabs(data[ii][jj][kk][ll][ndx]);
    }
    return m;
}


void calc_dw(double * w, double ** dw, double ** ddw, double ** kxy, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0)
{
    const int X = 0;
    const int Y = 1;
    const int XX = 0;
    const int XY = 1;
    const int YY = 2;

    const int N1c = N1/2+1;

    fftw_complex * kw = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * kdwx = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * kdwy = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * kddwxx = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * kddwxy = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * kddwyy = fftw_alloc_complex(local_n0*N1c);

    fftw_plan planF_w = fftw_mpi_plan_dft_r2c_2d(N0, N1, w, kw, MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planB_kdwx = fftw_mpi_plan_dft_c2r_2d(N0, N1, kdwx, dw[X], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planB_kdwy = fftw_mpi_plan_dft_c2r_2d(N0, N1, kdwy, dw[Y], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planB_kddwxx = fftw_mpi_plan_dft_c2r_2d(N0, N1, kddwxx, ddw[XX], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planB_kddwxy = fftw_mpi_plan_dft_c2r_2d(N0, N1, kddwxy, ddw[XY], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planB_kddwyy = fftw_mpi_plan_dft_c2r_2d(N0, N1, kddwyy, ddw[YY], MPI_COMM_WORLD, FFTW_ESTIMATE);

    fftw_execute(planF_w);

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1c; j++)
    {
        int ndx = i*N1c + j;

        kdwx[ndx][Re] = -kxy[X][ndx] * kw[ndx][Im];
        kdwx[ndx][Im] =  kxy[X][ndx] * kw[ndx][Re];

        kdwy[ndx][Re] = -kxy[Y][ndx] * kw[ndx][Im];
        kdwy[ndx][Im] =  kxy[Y][ndx] * kw[ndx][Re];

        kddwxx[ndx][Re] = -kxy[X][ndx]*kxy[X][ndx] * kw[ndx][Re];
        kddwxx[ndx][Im] = -kxy[X][ndx]*kxy[X][ndx] * kw[ndx][Im];
 
        kddwyy[ndx][Re] = -kxy[Y][ndx]*kxy[Y][ndx] * kw[ndx][Re];
        kddwyy[ndx][Im] = -kxy[Y][ndx]*kxy[Y][ndx] * kw[ndx][Im];

        kddwxy[ndx][Re] = -kxy[X][ndx]*kxy[Y][ndx] * kw[ndx][Re];
        kddwxy[ndx][Im] = -kxy[X][ndx]*kxy[Y][ndx] * kw[ndx][Im];
    }

    fftw_execute(planB_kdwx);
    fftw_execute(planB_kdwy);
    fftw_execute(planB_kddwxx);
    fftw_execute(planB_kddwxy);
    fftw_execute(planB_kddwyy);

    normalize(dw[X], N0, N1, local_n0);
    normalize(dw[Y], N0, N1, local_n0);
    normalize(ddw[XX], N0, N1, local_n0);
    normalize(ddw[XY], N0, N1, local_n0);
    normalize(ddw[YY], N0, N1, local_n0);

    fftw_destroy_plan(planF_w);
    fftw_destroy_plan(planB_kdwx);
    fftw_destroy_plan(planB_kdwy);
    fftw_destroy_plan(planB_kddwxx);
    fftw_destroy_plan(planB_kddwxy);
    fftw_destroy_plan(planB_kddwyy);

    fftw_free(kw);
    fftw_free(kdwx);
    fftw_free(kdwy);
    fftw_free(kddwxx);
    fftw_free(kddwxy);
    fftw_free(kddwyy);
}

void calc_dNuterm(double * kappa, double * nu, double * Nuterm, double ** ddNuterm, double ** kxy, double * phi, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0, struct input_parameters ip)
{
    const int N1r = 2*(N1/2 + 1);

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        kappa[ndx] = (1-phi[ndx])*ip.kappa_el0 + phi[ndx]*ip.kappa_el1;
        nu[ndx] = (1-phi[ndx])*ip.nu_el0 + phi[ndx]*ip.nu_el1;
        Nuterm[ndx] = (1-nu[ndx]);
    } 
   
    const int X = 0;
    const int Y = 1;
    const int XX = 0;
    const int XY = 1;
    const int YY = 2;

    const int N1c = N1/2+1;

    fftw_complex * kNuterm = fftw_alloc_complex(local_n0*N1c);

    fftw_complex * kddNutermxx = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * kddNutermxy = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * kddNutermyy = fftw_alloc_complex(local_n0*N1c);

    fftw_plan planF_Nuterm = fftw_mpi_plan_dft_r2c_2d(N0, N1, Nuterm, kNuterm, MPI_COMM_WORLD, FFTW_ESTIMATE);
    
    fftw_plan planB_kddNutermxx = fftw_mpi_plan_dft_c2r_2d(N0, N1, kddNutermxx, ddNuterm[XX], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planB_kddNutermxy = fftw_mpi_plan_dft_c2r_2d(N0, N1, kddNutermxy, ddNuterm[XY], MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planB_kddNutermyy = fftw_mpi_plan_dft_c2r_2d(N0, N1, kddNutermyy, ddNuterm[YY], MPI_COMM_WORLD, FFTW_ESTIMATE);

    fftw_execute(planF_Nuterm);

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1c; j++)
    {
        int ndx = i*N1c + j;

        kddNutermxx[ndx][Re] = -kxy[X][ndx]*kxy[X][ndx] * kNuterm[ndx][Re];
        kddNutermxx[ndx][Im] = -kxy[X][ndx]*kxy[X][ndx] * kNuterm[ndx][Im];
 
        kddNutermyy[ndx][Re] = -kxy[Y][ndx]*kxy[Y][ndx] * kNuterm[ndx][Re];
        kddNutermyy[ndx][Im] = -kxy[Y][ndx]*kxy[Y][ndx] * kNuterm[ndx][Im];

        kddNutermxy[ndx][Re] = -kxy[X][ndx]*kxy[Y][ndx] * kNuterm[ndx][Re];
        kddNutermxy[ndx][Im] = -kxy[X][ndx]*kxy[Y][ndx] * kNuterm[ndx][Im];
    }

    fftw_execute(planB_kddNutermxx);
    fftw_execute(planB_kddNutermxy);
    fftw_execute(planB_kddNutermyy);

    normalize(ddNuterm[XX], N0, N1, local_n0);
    normalize(ddNuterm[XY], N0, N1, local_n0);
    normalize(ddNuterm[YY], N0, N1, local_n0);

    fftw_destroy_plan(planF_Nuterm);
    fftw_destroy_plan(planB_kddNutermxx);
    fftw_destroy_plan(planB_kddNutermxy);
    fftw_destroy_plan(planB_kddNutermyy);

    fftw_free(kNuterm);
    fftw_free(kddNutermxx);
    fftw_free(kddNutermxy);
    fftw_free(kddNutermyy);
}

void calc_kbend(double * kappa, double ** ddw, double ** ddNuterm, double * bend_1, double * bend_2, fftw_complex * kbend_1, fftw_complex * kbend_2, ptrdiff_t local_n0, ptrdiff_t N0, ptrdiff_t N1)
{
    const int N1r = 2*(N1/2+1);
    const int X = 0;
    const int Y = 1;
    const int XX = 0;
    const int XY = 1;
    const int YY = 2;

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;

        bend_1[ndx] = kappa[ndx]*(ddw[XX][ndx]+ddw[YY][ndx]); 
        bend_2[ndx] = 2*ddNuterm[XY][ndx]*ddw[XY][ndx] - ddNuterm[XX][ndx]*ddw[YY][ndx] - ddNuterm[YY][ndx]*ddw[XX][ndx];
}
      fftw_execute(planF_bend_1);
      fftw_execute(planF_bend_2);
    
}

///////////////////////////////////////////////////////////////////////////////////////////////
void calc_dFdw(double * dFdw, fftw_complex * kdFdw_1, fftw_complex * kdFdw_2, fftw_complex * kdFdw_3, fftw_complex * kdFdw_4, fftw_complex * ktemp0_1, fftw_complex * ktemp0_2, fftw_complex * ktemp1_1, fftw_complex * ktemp1_2, double * kappa, double * w, double ** dw, fftw_complex * kbend_1, fftw_complex * kbend_2, double **** lam, double *** eps, double ** sigbar, double *** sig00, double ** kxy, ptrdiff_t local_n0, ptrdiff_t N0, ptrdiff_t N1)
// Calculate the chemical potentail of the out-of-plane displacement that will be used for evolution

// dFdw[ndx] is the variational derivative (chemical potential) of the out-of-plane displacement
// w[ndx], dw[i][ndx] are the out-of-plane displacements and its first derivatives
// lam[i][j][k][l] is the elastic stiffness tensor (lambda)
// eps[i][j][ndx] is the heterogeneous strain 0.5(u_{ij} + u_{ji})
// epsbar[i][ndx] is the homogeneous strain on the system
// s0n2[p][ij][ndx] is the product sig0(p,r) * eta(p) 
// kxy[i][ndx] are the k-vectors for calculating derivatives in k-space
// kappa is the bending modulus
///////////////////////////////////////////////////////////////////////////////////////////////
{
    const int N1c = N1/2 + 1;
    const int N1r = 2*(N1/2+1);

    const int X = 0;
    const int Y = 1;

    // allocate temporary memory
    double * temp0 = fftw_alloc_real(local_n0*N1r);
    double * temp1 = fftw_alloc_real(local_n0*N1r);

    fftw_complex * ktemp0 = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * ktemp1 = fftw_alloc_complex(local_n0*N1c);

    fftw_complex * kdFdw  = fftw_alloc_complex(local_n0*N1c);
    fftw_complex * kw     = fftw_alloc_complex(local_n0*N1c);

    // initialize fast fourier transforms

    fftw_plan planF_temp0 = fftw_mpi_plan_dft_r2c_2d(N0, N1, temp0, ktemp0, MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planF_temp1 = fftw_mpi_plan_dft_r2c_2d(N0, N1, temp1, ktemp1, MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planF_w    = fftw_mpi_plan_dft_r2c_2d(N0, N1, w, kw, MPI_COMM_WORLD, FFTW_ESTIMATE);

    fftw_plan planB_dFdw = fftw_mpi_plan_dft_c2r_2d(N0, N1, kdFdw, dFdw, MPI_COMM_WORLD, FFTW_ESTIMATE);

   // just for test
    double * temp0_1 = fftw_alloc_real(local_n0*N1r);
    double * temp1_1 = fftw_alloc_real(local_n0*N1r);
    double * temp0_2 = fftw_alloc_real(local_n0*N1r);
    double * temp1_2 = fftw_alloc_real(local_n0*N1r);
    fftw_plan planF_temp0_1 = fftw_mpi_plan_dft_r2c_2d(N0, N1, temp0_1, ktemp0_1, MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planF_temp0_2 = fftw_mpi_plan_dft_r2c_2d(N0, N1, temp0_2, ktemp0_2, MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planF_temp1_1 = fftw_mpi_plan_dft_r2c_2d(N0, N1, temp1_1, ktemp1_1, MPI_COMM_WORLD, FFTW_ESTIMATE);
    fftw_plan planF_temp1_2 = fftw_mpi_plan_dft_r2c_2d(N0, N1, temp1_2, ktemp1_2, MPI_COMM_WORLD, FFTW_ESTIMATE);

    // do some of the calculations in real space before taking derivatives
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        temp0[ndx] = 0;
        temp1[ndx] = 0;

        temp0_1[ndx] = 0;    //test
        temp0_2[ndx] = 0;
        temp1_1[ndx] = 0;
        temp1_2[ndx] = 0;

        for (int ii=0; ii<2; ii++)
        {
            temp0[ndx] += (sigbar[ii][0] - sig00[ii][0][ndx])*dw[ii][ndx];
            temp1[ndx] += (sigbar[ii][1] - sig00[ii][1][ndx])*dw[ii][ndx];

            temp0_1[ndx] += (sigbar[ii][0] - sig00[ii][0][ndx])*dw[ii][ndx];  //test
            temp1_1[ndx] += (sigbar[ii][1] - sig00[ii][1][ndx])*dw[ii][ndx];  
        }

        for (int ii=0; ii<2; ii++)
        for (int kk=0; kk<2; kk++)
        for (int ll=0; ll<2; ll++)
        {
            temp0[ndx] += lam[ii][0][kk][ll]*dw[ii][ndx]*(eps[ii][0][ndx] + 0.5*dw[kk][ndx]*dw[ll][ndx]);
            temp1[ndx] += lam[ii][1][kk][ll]*dw[ii][ndx]*(eps[ii][1][ndx] + 0.5*dw[kk][ndx]*dw[ll][ndx]);

            temp0_2[ndx] += lam[ii][0][kk][ll]*dw[ii][ndx]*(eps[ii][0][ndx] + 0.5*dw[kk][ndx]*dw[ll][ndx]);  //test
            temp1_2[ndx] += lam[ii][1][kk][ll]*dw[ii][ndx]*(eps[ii][1][ndx] + 0.5*dw[kk][ndx]*dw[ll][ndx]);
        }
        
    }

    // forward tranform to k-space
    fftw_execute(planF_temp0);
    fftw_execute(planF_temp1);
    fftw_execute(planF_w);

    fftw_execute(planF_temp0_1);     //test
    fftw_execute(planF_temp1_2);
    fftw_execute(planF_temp0_1);
    fftw_execute(planF_temp1_2);

    // calculate the derivatives in k-space
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1c; j++)
    {
        int ndx = i*N1c + j;
        double k4x = kxy[X][ndx]*kxy[X][ndx]*kxy[X][ndx]*kxy[X][ndx];
        double k4y = kxy[Y][ndx]*kxy[Y][ndx]*kxy[Y][ndx]*kxy[Y][ndx];

        kdFdw_1[ndx][Re] = kxy[X][ndx]*ktemp0[ndx][Im];
        kdFdw_2[ndx][Re] = kxy[Y][ndx]*ktemp1[ndx][Im];
        kdFdw_3[ndx][Re] = (k4x + k4y)*kw[ndx][Re]*0.22;
        //kdFdw_4[ndx][Re] = kbend_2[ndx][Re];
        kdFdw_1[ndx][Im] = -kxy[X][ndx]*ktemp0[ndx][Re];
        kdFdw_2[ndx][Im] = -kxy[Y][ndx]*ktemp1[ndx][Re];
        kdFdw_3[ndx][Im] = (k4x + k4y)*kw[ndx][Im]*0.22;
        //kdFdw_4[ndx][Im] = kbend_2[ndx][Im];
        kdFdw[ndx][Re] = kdFdw_1[ndx][Re] + kdFdw_2[ndx][Re] + kdFdw_3[ndx][Re];// + kdFdw_4[ndx][Re];
        kdFdw[ndx][Im] = kdFdw_1[ndx][Im] + kdFdw_2[ndx][Im] + kdFdw_3[ndx][Im];// + kdFdw_4[ndx][Im];


/*        kdFdw[ndx][Re] =  kxy[X][ndx]*ktemp0[ndx][Im] + kxy[Y][ndx]*ktemp1[ndx][Im] + (k4x + k4y)*kappa[ndx]*kw[ndx][Re] + kbend_2[ndx][Re];
        kdFdw[ndx][Im] = -kxy[X][ndx]*ktemp0[ndx][Re] - kxy[Y][ndx]*ktemp1[ndx][Re] + (k4x + k4y)*kappa[ndx]*kw[ndx][Im] + kbend_2[ndx][Im]; */
    }

    // inverse fourier transform kdFdw -> dFdw
    fftw_execute(planB_dFdw); 

    normalize(dFdw, N0, N1, local_n0);

    // free memory
    fftw_free(temp0);
    fftw_free(temp1);
    fftw_free(ktemp0);
    fftw_free(ktemp1);
    fftw_free(kdFdw);
    fftw_free(kw);

    fftw_destroy_plan(planF_w);
    fftw_destroy_plan(planF_temp0);
    fftw_destroy_plan(planF_temp1);
    fftw_destroy_plan(planB_dFdw);

    fftw_free(temp0_1);    //test
    fftw_free(temp0_2);
    fftw_free(temp1_1);
    fftw_free(temp1_2);
    fftw_destroy_plan(planF_temp0_1);
    fftw_destroy_plan(planF_temp0_2);
    fftw_destroy_plan(planF_temp1_1);
    fftw_destroy_plan(planF_temp1_2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void update_w(double * w, double * w_old, double * w_new, double * dFdw, ptrdiff_t local_n0, ptrdiff_t N1, struct input_parameters ip)
// step the out-of-plane displacement in time using the evolution wave equation
// HERE I SET ALL w=0 FIRST
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
    const int N1r = 2*(N1/2+1);
    double dtw = ip.dt/20.0;

    double dtg = 0.5*dtw*ip.gamma;
    double dtg2 = 1.0/(1.0+dtg);
    double dta2 = dtw*dtw*ip.alpha*ip.alpha;

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;

        w_new[ndx] = dtg2*(2*w[ndx] + (dtg-1)*w_old[ndx] - dta2*dFdw[ndx]);

        w_old[ndx] = w[ndx];
        w[ndx] = w_new[ndx];
    }
}

void add_w_noise(double * w, ptrdiff_t local_n0, ptrdiff_t N1)
{
    const int N1r = 2*(N1/2+1);
    double epdt2 = 0.00004;

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;
        double rnum = rand() / (double) RAND_MAX;
        rnum = 2*rnum - 1;
        w[ndx] = w[ndx] + epdt2*rnum;
    }
}

void calc_ksig00(double *** eps00, double *** sig00,double **** lam, ptrdiff_t local_n0, ptrdiff_t N1)
{    
    const int N1r = 2*(N1/2+1);
     
    for (int i=0; i<local_n0; i++)
        for (int j=0; j<N1; j++)
        {
            int ndx = i*N1r + j;

            for (int ii=0; ii<2; ii++)
            for (int jj=0; jj<2; jj++)
            {
                sig00[ii][jj][ndx] = 0.0;

                for (int kk=0; kk<2; kk++)
                for (int ll=0; ll<2; ll++)
                    sig00[ii][jj][ndx] += lam[ii][jj][kk][ll] * eps00[kk][ll][ndx];
            }
        }
     // sig00 -> ksig00
    for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
        fftw_execute(planF_sig00[i][j]);

}

void calc_homo_eps00(double **** lam, double *** homo_eps00, double ** epsbar, double *** eps00, double *** ift_eps00, fftw_complex *** ft_eps00, double *** G, double ** kxy, fftw_complex *** ksig00, ptrdiff_t N0, ptrdiff_t N1, ptrdiff_t local_n0)
{   const int N1c = N1/2+1;
    const int N1r = 2*(N1/2+1);
  
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1c; j++)
  {
    int ndx = i*N1c + j;

     for (int kk=0; kk<2; kk++)
     for (int ll=0; ll<2; ll++)
    { 
      ft_eps00[kk][ll][ndx][Re] = 0;
      ft_eps00[kk][ll][ndx][Im] = 0;
       
       for (int nn=0; nn<2; nn++)
       for (int mm=0; mm<2; mm++)
  
      {  
         ft_eps00[kk][ll][ndx][Re] += (kxy[kk][ndx] * G[ll][mm][ndx] + kxy[ll][ndx] * G[kk][mm][ndx]) * kxy[nn][ndx] * ksig00[mm][nn][ndx][Re];
         ft_eps00[kk][ll][ndx][Im] += (kxy[kk][ndx] * G[ll][mm][ndx] + kxy[ll][ndx] * G[kk][mm][ndx]) * kxy[nn][ndx] * ksig00[mm][nn][ndx][Im];
   }
  }
 }
    for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
     {
      fftw_execute(planB_ft_eps00[i][j]);
      normalize(ift_eps00[i][j], N0, N1, local_n0);
}
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
      int ndx = i*N1r + j;
      for (int ii=0; ii<2; ii++)
      for (int jj=0; jj<2; jj++) 
      {
      homo_eps00[ii][jj][ndx] = 0;
      for (int kk=0; kk<2; kk++)
      for (int ll=0; ll<2; ll++)  
      homo_eps00[ii][jj][ndx] += -0.5*lam[ii][jj][kk][ll]*ift_eps00[kk][ll][ndx] - lam[ii][jj][kk][ll]*epsbar[kk][ll] + lam[ii][jj][kk][ll]*eps00[kk][ll][ndx];
}}
}


void calc_inhom_eps00(double *** inhom_eps00, double *** inhom_eps00_1, double *** inhom_eps00_2, double *** inhom_eps00_3, double *** inhom_eps00_4, double **** lam, double ***** delta_S, double *** eps00, double **** eps0, double ** eta, ptrdiff_t local_n0, ptrdiff_t N1)
{    
    const int N1r = 2*(N1/2+1);     

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
        {
            int ndx = i*N1r + j;
     for (int ii=0; ii<2; ii++)
     for (int jj=0; jj<2; jj++)
     { 
        inhom_eps00[ii][jj][ndx] = 0;
        inhom_eps00_1[ii][jj][ndx] = 0;
        inhom_eps00_2[ii][jj][ndx] = 0;
        inhom_eps00_3[ii][jj][ndx] = 0;
        inhom_eps00_4[ii][jj][ndx] = 0;

       for (int kk=0; kk<2; kk++)
       for (int ll=0; ll<2; ll++)
          inhom_eps00_1[ii][jj][ndx] += (-1)*lam[ii][jj][kk][ll]*eps00[kk][ll][ndx];

       for (int kk=0; kk<2; kk++)
       for (int ll=0; ll<2; ll++)
       for (int p=0; p<3; p++)
             inhom_eps00_2[ii][jj][ndx] += (1)*lam[ii][jj][kk][ll]*eps0[p][kk][ll][ndx]*eta[p][ndx]*eta[p][ndx];

       for (int kk=0; kk<2; kk++)
       for (int ll=0; ll<2; ll++)
       for (int mm=0; mm<2; mm++)
       for (int nn=0; nn<2; nn++)
       for (int ss=0; ss<2; ss++)
       for (int tt=0; tt<2; tt++) 
          inhom_eps00_3[ii][jj][ndx] += (1)*lam[ii][jj][kk][ll]*fabs(delta_S[kk][ll][mm][nn][ndx])*lam[mm][nn][ss][tt]*eps00[ss][tt][ndx];
           
       for (int kk=0; kk<2; kk++)
       for (int ll=0; ll<2; ll++)
       for (int mm=0; mm<2; mm++)
       for (int nn=0; nn<2; nn++)
       for (int ss=0; ss<2; ss++)
       for (int tt=0; tt<2; tt++) 
       for (int pp=0; pp<3; pp++)
            inhom_eps00_4[ii][jj][ndx] += (-1)*lam[ii][jj][kk][ll]*fabs(delta_S[kk][ll][mm][nn][ndx])*lam[mm][nn][ss][tt]*eps0[pp][ss][tt][ndx]*eta[pp][ndx]*eta[pp][ndx];
       
       inhom_eps00[ii][jj][ndx] = inhom_eps00_1[ii][jj][ndx] + inhom_eps00_2[ii][jj][ndx] + inhom_eps00_3[ii][jj][ndx] + inhom_eps00_4[ii][jj][ndx];
}
}
}



void calc_dFdeps00(double *** dFdeps00, double **** lam, double *** homo_eps00, double *** inhom_eps00, ptrdiff_t local_n0, ptrdiff_t N1)
{
     const int N1r = 2*(N1/2+1);
     
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
      {
        int ndx = i*N1r + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
            dFdeps00[ii][jj][ndx] = homo_eps00[ii][jj][ndx] + inhom_eps00[ii][jj][ndx];
       
      }
}



double update_eps00(double *** eps00, double *** eps00_new, double *** eps00_old, double *** dFdeps00, ptrdiff_t local_n0, ptrdiff_t N1, struct input_parameters ip)
{
     double dtg = 0.5*ip.dt*ip.gamma_eps;
    double dtg2 = 1.0/(1.0+dtg);
    double dta2 = ip.dt*ip.dt*ip.alpha_eps*ip.alpha_eps;
     const int N1r = 2*(N1/2+1);
    double change_eps00_max = 0;

    for (int ii=0; ii<2; ii++)
    for (int jj=0; jj<2; jj++)
    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*N1r + j;

        //eps00_new[ii][jj][ndx] = dtg2*(2*eps00[ii][jj][ndx] + (dtg-1)*eps00_old[ii][jj][ndx] - dta2*dFdeps00[ii][jj][ndx]);
        eps00_new[ii][jj][ndx] = (-1)*ip.L_eps00*ip.dt*dFdeps00[ii][jj][ndx] + eps00[ii][jj][ndx];

        double delta = fabs( eps00_new[ii][jj][ndx] - eps00[ii][jj][ndx] );
        change_eps00_max = std::max(change_eps00_max, delta);

        //eps00_old[ii][jj][ndx] = eps00[ii][jj][ndx];
        eps00[ii][jj][ndx] = eps00_new[ii][jj][ndx];
    }

    return change_eps00_max;
}

void make_it_2D(double * w, double ** dw, double ** ddw, ptrdiff_t alloc_local)
{
     for (int i=0; i<2*alloc_local; i++)
     {  w[i] = 0;
        for (int ii=0; ii<2; ii++)
        dw[ii][i] = 0;
        for (int jj=0; jj<3; jj++)
        ddw[jj][i] = 0;  
}
}

double check(double **** lam, double **** eps0, double ** eta, double *** eps00, int local_n0, int N1)
{
    double m = 0;
    double delta;

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
{
        for (int kk=0; kk<2; kk++)
        for (int ll=0; ll<2; ll++)
        {   
            delta = 0;
            for(int p=0; p<3; p++)
            delta += lam[ii][jj][kk][ll]*eps0[p][kk][ll][ndx]*eta[p][ndx]*eta[p][ndx];
            delta += (-1)*lam[ii][jj][kk][ll]*eps00[kk][ll][ndx];
         if (fabs(delta) > m) m = fabs(delta);
        
}}
    }
      
    
    return m;
    
}

double avg_check ( double **** lam, double **** eps0, double ** eta, double *** eps00, int local_n0, int N1 )
{
    double m = 0;
    double delta;

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
{
        for (int kk=0; kk<2; kk++)
        for (int ll=0; ll<2; ll++)
        {   
            delta = 0;
            for(int p=0; p<3; p++)
            delta += lam[ii][jj][kk][ll]*eps0[p][kk][ll][ndx]*eta[p][ndx]*eta[p][ndx];
            delta += (-1)*lam[ii][jj][kk][ll]*eps00[kk][ll][ndx];
         m += delta;
        
}}
        m = m/local_n0/N1;
    }
      
    
    return m;
}

double avg_abs_check ( double **** lam, double **** eps0, double ** eta, double *** eps00, int local_n0, int N1 )
{
    double m = 0;
    double delta;

    for (int i=0; i<local_n0; i++)
    for (int j=0; j<N1; j++)
    {
        int ndx = i*2*(N1/2+1) + j;
        for (int ii=0; ii<2; ii++)
        for (int jj=0; jj<2; jj++)
{
        for (int kk=0; kk<2; kk++)
        for (int ll=0; ll<2; ll++)
        {   

            delta = 0;
            for(int p=0; p<3; p++)
            delta += lam[ii][jj][kk][ll]*eps0[p][kk][ll][ndx]*eta[p][ndx]*eta[p][ndx];
            delta += (-1)*lam[ii][jj][kk][ll]*eps00[kk][ll][ndx];
         m += fabs(delta);
        
}}
        m = m/local_n0/N1;
    }
      
    
    return m;
}

void calc_sigbar(double ** sigbar, double ** epsbar, double **** lam)
{
      for(int ii=0; ii<2; ii++)
      for(int jj=0; jj<2; jj++)
      { sigbar[ii][jj] = 0;
      for(int kk=0; kk<2; kk++)
      for(int ll=0; ll<2; ll++)
        sigbar[ii][jj] += lam[ii][jj][kk][ll]*epsbar[kk][ll];

}
}

int main(int argc, char ** argv)
{

    MPI_Init(&argc, &argv);
    fftw_mpi_init();

    srand(time(NULL));

    struct input_parameters ip;

    ParameterFile pf;
    pf.readParameters("input.txt");
    pf.unpack("Nx", ip.Nx);
    pf.unpack("Ny", ip.Ny);
    pf.unpack("dx", ip.dx);
    pf.unpack("dt", ip.dt);

    pf.unpack("mu_el0", ip.mu_el0);
    pf.unpack("nu_el0", ip.nu_el0);
    pf.unpack("mu_el1", ip.mu_el1);
    pf.unpack("nu_el1", ip.nu_el1);
    pf.unpack("nsteps", ip.nsteps);
    pf.unpack("out_freq", ip.out_freq);

    pf.unpack("epsx", ip.epsx);
    pf.unpack("epsy", ip.epsy);

    pf.unpack("gamma", ip.gamma);
    pf.unpack("gamma", ip.gamma_eps);
    pf.unpack("kappa_el0", ip.kappa_el0);
    pf.unpack("kappa_el1", ip.kappa_el1);
    pf.unpack("alpha", ip.alpha);
    pf.unpack("alpha", ip.alpha_eps);
    pf.unpack("L_eps00", ip.L_eps00);
    pf.unpack("change_etap_thresh", ip.change_etap_thresh);
    pf.unpack("change_eps00_thresh", ip.change_eps00_thresh);
    pf.unpack("beta", ip.beta);

    pf.unpack("M0_chem_a", ip.M0_chem_a);
    pf.unpack("M0_chem_b", ip.M0_chem_b);
    pf.unpack("M0_chem_c", ip.M0_chem_c);

    pf.unpack("M1_chem_a", ip.M1_chem_a);
    pf.unpack("M1_chem_b", ip.M1_chem_b);
    pf.unpack("M1_chem_c", ip.M1_chem_c);

    pf.unpack("M0_2H_a", ip.M0_2H_a);
    pf.unpack("M0_2H_b", ip.M0_2H_b);
    pf.unpack("M0_Tp_a", ip.M0_Tp_a);
    pf.unpack("M0_Tp_b", ip.M0_Tp_b);

    pf.unpack("M1_2H_a", ip.M1_2H_a);
    pf.unpack("M1_2H_b", ip.M1_2H_b);
    pf.unpack("M1_Tp_a", ip.M1_Tp_a);
    pf.unpack("M1_Tp_b", ip.M1_Tp_b);

    pf.unpack("M0_norm", ip.M0_norm);
    pf.unpack("M1_norm", ip.M1_norm);

    ptrdiff_t local_n0;
    ptrdiff_t local_0_start;
    ptrdiff_t N0 = (ptrdiff_t) ip.Nx;
    ptrdiff_t N1 = (ptrdiff_t) ip.Ny;

    ptrdiff_t alloc_local = fftw_mpi_local_size_2d(N0, N1/2+1, MPI_COMM_WORLD, &local_n0, &local_0_start);

    double * ux;        // in-plane displacement in the x-direction     ux[ndx]
    double * uy;        // in-plane displacement in the y-direction     uy[ndx]
    double * phi;       // material composition                         phi[ndx]
    double * lsf;       // level set function - used for initialization lsf[ndx]
    double ** eta;      // orientation order parmater                   eta[p][ndx]
    double ** eta_old;  // previous step order parmater                 eta[p][ndx]
    double ** eta_new;  // next step order parmater                     eta[p][ndx]    

    double **** lam;    // stiffness tensor         lam[i][j][k][l]
    double *** G;       // greens functions         G[i][j][ndx]
    double ** kxy;      // fourier k-vectors        kxy[i][ndx] 
    double **** epsT;   // transformation strains   epsT[M][p][i][j]
    double **** eps0;   // transformation strains   eps0[p][i][j][ndx]
    double **** sig0;   // transformation stress    sig0[p][i][j][ndx]
    double *** sigeps;  // sig0*eps0                sigeps[p][q][ndx]
    double ** epsbar;   // homogeneous strain       epsbar[i][j]
    double ** sigbar;
    double *** eps;     // heterogeneous strain     eps[p][i][j][ndx]
    double *** s0n2;    // sig0*eta                 s0n2[p][i+j][ndx]
    double ** chem;     // chemical potential       chem[p][ndx]
    double * w;         // out-of-plane bending     w[ndx]
    double ** dw;       // first derivatives        dw[i][ndx]
    double ** ddw;      // second derivatives       ddw[i+j][ndx]
    double * dFdw;      // delta F / delta w        dFdw[ndx]                    
    double * w_old;
    double * w_new;
    double * Nuterm;     // (1-nu)                 Nuterm[ndx]
    double ** ddNuterm;  // second derivatives       ddNuterm[i+j][ndx]
    double * kappa;     // interpolate K            kappa[ndx]
    double * nu;        // interpolate nu           nu[ndx] 

    double *** eps00;                    //eps00[i][j][ndx]
    double *** sig00;                    //sig00[i][j][ndx] 
    double *** ift_eps00;                //ift_eps00[i][j][ndx]
    double *** homo_eps00;                //homo_eps00[i][j][ndx]
    double *** inhom_eps00;              //inhom_eps00[i][j][ndx] 
    double *** inhom_eps00_1;     
    double *** inhom_eps00_2; 
    double *** inhom_eps00_3; 
    double *** inhom_eps00_4; 

    double *** dFdeps00;                 //dFdeps00[i][j][ndx]
    double *** eps00_new;                //eps00_new[i][j][ndx] 
    double *** eps00_old;                //eps00_old[i][j][ndx]     
    
    double ***** delta_S;                //delta_S[i][j][k][l][ndx]
    double * bend_1;
    double * bend_2;

    fftw_complex *** ft_eps00;           //ft_eps00[i][j][ndx]
    fftw_complex *** ksig00;             //ksig00[i][j][ndx]
    fftw_complex *** ks0n2; // fourier transform of s0n2, ks0n2[p][j+k][ndx]
    fftw_complex * kbend_1; // FT of first term of bending energy
    fftw_complex * kbend_2;

    // allocate all of the memory needed for the simulation

    eps00   = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    sig00   = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    ksig00  = (fftw_complex ***) kd_alloc2(sizeof(fftw_complex), 3, 2, 2, alloc_local);
    ift_eps00   = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    homo_eps00   = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    ft_eps00  = (fftw_complex ***) kd_alloc2(sizeof(fftw_complex), 3, 2, 2, alloc_local);
    inhom_eps00 = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    inhom_eps00_1 = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    inhom_eps00_2 = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    inhom_eps00_3 = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    inhom_eps00_4 = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    delta_S = (double *****) kd_alloc2(sizeof(double), 5, 2, 2, 2, 2, 2*alloc_local);
    dFdeps00 = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    eps00_new = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    eps00_old = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);

    lam     = (double ****) kd_alloc2(sizeof(double), 4, 2, 2, 2, 2);
    G       = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, alloc_local);
    kxy     = (double **)   kd_alloc2(sizeof(double), 2, 2, alloc_local);
    epsT    = (double ****) kd_alloc2(sizeof(double), 4, 2, 3, 2, 2);
    eps0    = (double ****) kd_alloc2(sizeof(double), 4, 3, 2, 2, 2*alloc_local);
    sig0    = (double ****) kd_alloc2(sizeof(double), 4, 3, 2, 2, 2*alloc_local);
    sigeps  = (double ***)  kd_alloc2(sizeof(double), 3, 3, 3, 2*alloc_local);
    epsbar  = (double **)   kd_alloc2(sizeof(double), 2, 2, 2);
    sigbar  = (double **)   kd_alloc2(sizeof(double), 2, 2, 2);
    eps     = (double ***)  kd_alloc2(sizeof(double), 3, 2, 2, 2*alloc_local);
    s0n2    = (double ***)  kd_alloc2(sizeof(double), 3, 3, 3, 2*alloc_local);
    ks0n2   = (fftw_complex ***) kd_alloc2(sizeof(fftw_complex), 3, 3, 3, alloc_local);
    
        

    w       = fftw_alloc_real(2*alloc_local);
    w_old   = fftw_alloc_real(2*alloc_local);
    w_new   = fftw_alloc_real(2*alloc_local);
    dFdw    = fftw_alloc_real(2*alloc_local);
    dw      = (double **) kd_alloc2(sizeof(double), 2, 2, 2*alloc_local);
    ddw     = (double **) kd_alloc2(sizeof(double), 2, 3, 2*alloc_local);
    Nuterm   = fftw_alloc_real(2*alloc_local);
    ddNuterm = (double **) kd_alloc2(sizeof(double), 2, 3, 2*alloc_local);
    kappa   = fftw_alloc_real(2*alloc_local);
    nu      = fftw_alloc_real(2*alloc_local);
    bend_1  = fftw_alloc_real(2*alloc_local);         
    bend_2  = fftw_alloc_real(2*alloc_local);
    kbend_1 = fftw_alloc_complex(alloc_local);
    kbend_2 = fftw_alloc_complex(alloc_local); 

    fftw_complex * kdFdw_1;                               //test
    fftw_complex * kdFdw_2;
    fftw_complex * kdFdw_3;
    fftw_complex * kdFdw_4; 
    kdFdw_1 = fftw_alloc_complex(alloc_local);
    kdFdw_2 = fftw_alloc_complex(alloc_local);  
    kdFdw_3 = fftw_alloc_complex(alloc_local);
    kdFdw_4 = fftw_alloc_complex(alloc_local);

    fftw_complex * ktemp0_1;
    ktemp0_1 = fftw_alloc_complex(alloc_local);
    fftw_complex * ktemp0_2;
    ktemp0_2 = fftw_alloc_complex(alloc_local);
    fftw_complex * ktemp1_1;
    ktemp1_1 = fftw_alloc_complex(alloc_local);
    fftw_complex * ktemp1_2;
    ktemp1_2 = fftw_alloc_complex(alloc_local);

    chem = new double * [3];
    chem[0] = fftw_alloc_real(2*alloc_local);
    chem[1] = fftw_alloc_real(2*alloc_local);
    chem[2] = fftw_alloc_real(2*alloc_local);

    eta = new double * [3];
    eta[0] = fftw_alloc_real(2*alloc_local);
    eta[1] = fftw_alloc_real(2*alloc_local);
    eta[2] = fftw_alloc_real(2*alloc_local);

    eta_old = new double * [3];
    eta_old[0] = fftw_alloc_real(2*alloc_local);
    eta_old[1] = fftw_alloc_real(2*alloc_local);
    eta_old[2] = fftw_alloc_real(2*alloc_local);

    eta_new = new double * [3];
    eta_new[0] = fftw_alloc_real(2*alloc_local);
    eta_new[1] = fftw_alloc_real(2*alloc_local);
    eta_new[2] = fftw_alloc_real(2*alloc_local);

    double ** lap = new double * [3];
    lap[0] = fftw_alloc_real(2*alloc_local);
    lap[1] = fftw_alloc_real(2*alloc_local);
    lap[2] = fftw_alloc_real(2*alloc_local);

    fftw_complex ** ketapsq = new fftw_complex * [3];
    ketapsq[0] = fftw_alloc_complex(alloc_local);
    ketapsq[1] = fftw_alloc_complex(alloc_local);
    ketapsq[2] = fftw_alloc_complex(alloc_local);

    fftw_complex ** keta = new fftw_complex * [3];
    keta[0] = fftw_alloc_complex(alloc_local);
    keta[1] = fftw_alloc_complex(alloc_local);
    keta[2] = fftw_alloc_complex(alloc_local);

    fftw_complex ** klap = new fftw_complex * [3];
    klap[0] = fftw_alloc_complex(alloc_local);
    klap[1] = fftw_alloc_complex(alloc_local);
    klap[2] = fftw_alloc_complex(alloc_local);

    fftw_complex ** keps = new fftw_complex * [3];
    keps[0] = fftw_alloc_complex(alloc_local);
    keps[1] = fftw_alloc_complex(alloc_local);
    keps[2] = fftw_alloc_complex(alloc_local);

    ux = fftw_alloc_real(2*alloc_local);
    uy = fftw_alloc_real(2*alloc_local);
    phi = fftw_alloc_real(2*alloc_local);
    lsf = fftw_alloc_real(local_n0*N1);
    

    fftw_complex ** ku = new fftw_complex * [2];
    ku[0] = fftw_alloc_complex(alloc_local);
    ku[1] = fftw_alloc_complex(alloc_local);


    // initialize the necessary fourier transforms
    for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
        planF_sig00[i][j] = fftw_mpi_plan_dft_r2c_2d(N0, N1, sig00[i][j], ksig00[i][j], MPI_COMM_WORLD, FFTW_MEASURE);

    for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
        planB_ft_eps00[i][j] = fftw_mpi_plan_dft_c2r_2d(N0, N1, ft_eps00[i][j], ift_eps00[i][j], MPI_COMM_WORLD, FFTW_MEASURE);

    planF_eta[0] = fftw_mpi_plan_dft_r2c_2d(N0, N1, eta[0], keta[0], MPI_COMM_WORLD, FFTW_MEASURE);
    planF_eta[1] = fftw_mpi_plan_dft_r2c_2d(N0, N1, eta[1], keta[1], MPI_COMM_WORLD, FFTW_MEASURE);
    planF_eta[2] = fftw_mpi_plan_dft_r2c_2d(N0, N1, eta[2], keta[2], MPI_COMM_WORLD, FFTW_MEASURE);

    planB_lap[0] = fftw_mpi_plan_dft_c2r_2d(N0, N1, klap[0], lap[0], MPI_COMM_WORLD, FFTW_MEASURE);
    planB_lap[1] = fftw_mpi_plan_dft_c2r_2d(N0, N1, klap[1], lap[1], MPI_COMM_WORLD, FFTW_MEASURE);
    planB_lap[2] = fftw_mpi_plan_dft_c2r_2d(N0, N1, klap[2], lap[2], MPI_COMM_WORLD, FFTW_MEASURE);

    for (int p=0; p<3; p++)
    for (int i=0; i<3; i++)
        planF_s0n2[p][i] = fftw_mpi_plan_dft_r2c_2d(N0, N1, s0n2[p][i], ks0n2[p][i], MPI_COMM_WORLD, FFTW_MEASURE);

    planB_ux = fftw_mpi_plan_dft_c2r_2d(N0, N1, ku[0], ux, MPI_COMM_WORLD, FFTW_MEASURE);
    planB_uy = fftw_mpi_plan_dft_c2r_2d(N0, N1, ku[1], uy, MPI_COMM_WORLD, FFTW_MEASURE);

    plan_strain_xx = fftw_mpi_plan_dft_c2r_2d(N0, N1, keps[0], eps[0][0], MPI_COMM_WORLD, FFTW_MEASURE);
    plan_strain_yy = fftw_mpi_plan_dft_c2r_2d(N0, N1, keps[1], eps[1][1], MPI_COMM_WORLD, FFTW_MEASURE);
    plan_strain_xy = fftw_mpi_plan_dft_c2r_2d(N0, N1, keps[2], eps[0][1], MPI_COMM_WORLD, FFTW_MEASURE);

     planF_bend_1 = fftw_mpi_plan_dft_r2c_2d(N0, N1, bend_1, kbend_1, MPI_COMM_WORLD, FFTW_ESTIMATE);
     planF_bend_2 = fftw_mpi_plan_dft_r2c_2d(N0, N1, bend_2, kbend_2, MPI_COMM_WORLD, FFTW_ESTIMATE);

    // calculate the elastic parameters

    calc_greens_function(G, kxy, local_n0, local_0_start, N1, ip);
    calc_transformation_strains(epsT, ip);

    // initialize the system with in-plane heterogeneity

    initialize_phi_1(phi,local_n0,N1);
    //initialize_lsf_circle(lsf, local_n0, local_0_start, N1);
    //initialize_lsf_stripe(lsf, local_n0, local_0_start, N1);
    //initialize_lsf_zigzag(lsf, local_n0, local_0_start, N1);
    //diffuse_lsf(lsf, local_n0, N1);
    //copy_lsf(lsf, phi, local_n0, N1);

    for (int p=0; p<3; p++)
    for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
        interpolate(eps0[p][i][j], epsT[0][p][i][j], epsT[1][p][i][j], phi, local_n0, N1);

    calc_elastic_tensors(lam, eps0, sig0, sigeps, ip.mu_el0, ip.nu_el0, local_n0, N1, delta_S, ip.mu_el1, ip.nu_el1, phi);

    log_greens_function(G, kxy, local_n0, N1);
    log_elastic_tensors(lam, epsT);

    // initialize eta parameters to zero
    initialize(eta, eta_old, local_n0, N1);

    
    // initialize displacements to zero
    for (int i=0; i<2*alloc_local; i++) { ux[i]=0; uy[i]=0; w_old[i]=0; w[i]=0; }
    
    // initialize eps00 to 0
    for (int i=0; i<2*alloc_local; i++) 
    for (int ii=0; ii<2; ii++)
    for (int jj=0; jj<2; jj++)
    { eps00[ii][jj][i] = 0; 
      eps00_old[ii][jj][i] = 0;
      homo_eps00[ii][jj][i] = 0;
      inhom_eps00[ii][jj][i] = 0;
      dFdeps00[ii][jj][i] = 0;
}    

    /*initialize_lsf_circle(lsf, local_n0, local_0_start,N1);
    diffuse_lsf(lsf, local_n0, N1);
    copy_lsf(lsf, eta[0], local_n0, N1);
    copy_lsf(lsf, eta_old[0], local_n0, N1);*/



    // begin writing the output file
    H5File h5;
    h5.open("out.h5", "w");
    output("phi", phi, N0, N1, local_n0);
    h5.close();

    FILE * fp = fopen("area_fraction.dat", "w");
    fclose(fp);   

    std::ofstream myfile ("test.log", std::ofstream::out);

    // begin the simulation loop
    int frame = 0;
    for (int step=1; step<=ip.nsteps; step++)
    {
        // increase load on system
        epsbar[0][0] = ip.epsx * (step/(double)ip.nsteps);
        epsbar[1][1] = ip.epsy * (step/(double)ip.nsteps);
        epsbar[0][1] = 0;
        epsbar[1][0] = 0;
        calc_sigbar(sigbar,epsbar, lam);

        // iterative relaxation loop for eta_p parameters
        double change_etap_max = 1;
        double change_eps00_max = 1;
        double area_fraction;
        
        while (change_etap_max > ip.change_etap_thresh)
        {
            change_etap_max = 0;
            change_eps00_max = 1;
            //make_it_2D( w, dw, ddw, alloc_local);
          while (change_eps00_max > ip.change_eps00_thresh)
           { 
             change_eps00_max = 0;
             //introduce_noise_eps00(eps00, local_n0, N1);
             //Calculate fourier term in misfit strain equation
            calc_ksig00(eps00, sig00, lam, local_n0, N1);
            // Calculate homo term in misfit strain equation
            calc_homo_eps00(lam, homo_eps00, epsbar, eps00, ift_eps00, ft_eps00, G, kxy, ksig00, N0, N1, local_n0);

            // Calculate inhom terms in misfir strain equation
            calc_inhom_eps00(inhom_eps00, inhom_eps00_1, inhom_eps00_2, inhom_eps00_3, inhom_eps00_4, lam, delta_S, eps00, eps0, eta, local_n0, N1);

            // Calculate the chemical potential of misfit strain
            calc_dFdeps00(dFdeps00, lam, homo_eps00, inhom_eps00, local_n0, N1);

            // Step the eps00 parameters in time using evolution equation
            change_eps00_max = update_eps00(eps00, eps00_new, eps00_old, dFdeps00, local_n0, N1, ip);

            // share convergence info with all processes for parallel computation
            MPI_Allreduce(MPI_IN_PLACE, &change_eps00_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            area_fraction = calc_area(eta, local_n0, N0, N1, ip.M1_norm);
            printf("%8d cepmax=%12.10f, ceps00max=%12.10f, Af=%12.10f\n",step,change_etap_max,change_eps00_max,area_fraction);

            myfile << " " << std::endl;
            myfile << "step =" << step << std::endl;
            //myfile << "eps00 = " << max_eps00(eps00, local_n0, N1) << std::endl;
            //myfile << "eps00_avg = " << avg_eps00(eps00, local_n0, N1) << std::endl;
            //myfile << "sig00 = " << max_eps00(sig00, local_n0, N1) << std::endl;
            //myfile << "sig00_avg = " << avg_eps00(sig00, local_n0, N1) << std::endl;
            //myfile << "eps00_avg_abs = " << avg_abs_eps00(eps00, local_n0, N1) << std::endl;
            //myfile << "homo_eps00 = " << max_eps00(homo_eps00, local_n0, N1) << std::endl;
            //myfile << "inhom_eps00 = " << max_eps00(inhom_eps00, local_n0, N1) << std::endl;
            //myfile << "homo_eps00_avg = " << avg_eps00(homo_eps00, local_n0, N1) << std::endl;
            //myfile << "inhomo_eps00_avg = " << avg_eps00(inhom_eps00, local_n0, N1) << std::endl;
            //myfile << "dFdeps00_avg = " << max_eps00(dFdeps00, local_n0, N1) << std::endl;
            //myfile << "ceps00max = " << change_eps00_max<< std::endl;*/
             }

            // fourier transform the nonlinear term in displacement equation sig0_{jk}*eta_p^2
            calc_ks0n2(s0n2, sig0, eta, local_n0, N1);

            // calculate the displacement in k-space using greens function
            //calc_uxy(ux, uy, ku, G, kxy, ksig00, N0, N1, local_n0);
            calc_uxy_bending(ux, uy, ku, dw, ddw, lam, G, kxy, ksig00, N0, N1, local_n0);

            // calculate the heterogeneous strain (delta-epsilon) in k-space
            calc_eps(eps, keps, kxy, ku, N0, N1, local_n0);

           //myfile << "eps_avg = " << avg_eps00(eps, local_n0, N1) << std::endl;
           // myfile << "eps_max = " << max_eps00(eps, local_n0, N1) << std::endl;
            //myfile << "kbend_1_max = " << max_kbend(kbend_1, local_n0, N1) << std::endl;
            //myfile << "kbend_2_max = " << max_kbend(kbend_2, local_n0, N1) << std::endl;
            //myfile << "kdFdw_1_max = " << max_kbend(kdFdw_1, local_n0, N1) << std::endl;
            //myfile << "kdFdw_2_max = " << max_kbend(kdFdw_2, local_n0, N1) << std::endl;
            //myfile << "kdFdw_3_max = " << max_kbend(kdFdw_3, local_n0, N1) << std::endl;
            //myfile << "kdFdw_4_max = " << max_kbend(kdFdw_4, local_n0, N1) << std::endl;
            myfile << "ktemp0_1_max = " << max_kbend(ktemp0_1, local_n0, N1) << std::endl;
            myfile << "ktemp0_2_max = " << max_kbend(ktemp0_2, local_n0, N1) << std::endl;
            myfile << "ktemp1_1_max = " << max_kbend(ktemp1_1, local_n0, N1) << std::endl;
            myfile << "ktemp1_2_max = " << max_kbend(ktemp1_2, local_n0, N1) << std::endl;
          
           // myfile << "ksig00_max = " << max_ksig00(ksig00, local_n0, N1) << std::endl;
           // myfile << "G_max = " << max_G(G, local_n0, N1) << std::endl;
           // myfile << "kxy_max = " << max_kxy(kxy, local_n0, N1) << std::endl;            


            // introduce random noise into the eta parameters
            introduce_noise(eta, local_n0, N1);

            // calculate the laplacian of the eta parameters (for the gradient squared energy term)
            calc_lap(lap, klap, keta, kxy, N0, N1, local_n0);

            // calculate the chemical potential for the eta parameters 
            calc_chemical_potential(sig00, eps0, chem, eta, eps00, epsbar, sig0, eps, lap, phi, dw, local_n0, N1, ip);

            // step the eta parameters in time using the evolution wave equation
            change_etap_max = update_eta(eta, eta_old, eta_new, chem, local_n0, N1, ip);

            // share convergence info with all processes for parallel computation
            MPI_Allreduce(MPI_IN_PLACE, &change_etap_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

            //myfile << "eta = " << max_eta(eta, local_n0, N1) << std::endl;
            //myfile << "eta_avg = " << avg_eta(eta, local_n0, N1) << std::endl;
            //myfile << "dFdeta = " << max_eta(chem, local_n0, N1) << std::endl;
            //myfile << "dFdeta = " << avg_eta(chem, local_n0, N1) << std::endl;
            //std::cout << "delta_S = " << max_delta_S(delta_S, local_n0, N1) << std::endl;
            // the rest for out-of-plane displacements - in progress

            // calculate first derivatives of the out-of-plane displacement
            calc_dw(w, dw, ddw, kxy, N0, N1, local_n0);

            // calculate derivatives of (1-nu) and fourier transform of bending term
            //calc_dNuterm( kappa, nu, Nuterm, ddNuterm, kxy, phi, N0, N1, local_n0, ip);
            //calc_kbend( kappa, ddw, ddNuterm, bend_1, bend_2, kbend_1, kbend_2, local_n0, N0, N1);

            // calculate the chemical potential of out-of-plane displacement
            calc_dFdw(dFdw, kdFdw_1, kdFdw_2, kdFdw_3, kdFdw_4, ktemp0_1, ktemp0_2, ktemp1_1, ktemp1_2, kappa, w, dw, kbend_1, kbend_2, lam, eps, epsbar, sig00, kxy, local_n0, N0, N1);

            // step w in time using evolution wave equation  
            update_w(w, w_old, w_new, dFdw, local_n0, N1, ip);
            add_w_noise(w, local_n0, N1); 

            std::cout << "w = " << max_w(w, local_n0, N1) << std::endl;

            // calculate and output area - will change in future versions 
            area_fraction = calc_area(eta, local_n0, N0, N1, ip.M1_norm);
            printf("%8d cepmax=%12.10f, ceps00max=%12.10f, Af=%12.10f\n",step,change_etap_max,change_eps00_max,area_fraction);     
        }

        // eta_p parameters have reached a thermodynamic and mechanical equilibrium
        //myfile << "check = " << check(lam, eps0, eta, eps00, local_n0, N1) << std::endl;
        //myfile << "check_avg = " << avg_check(lam, eps0, eta, eps00, local_n0, N1) << std::endl;
        //myfile << "check_abs_avg = " << avg_abs_check(lam, eps0, eta, eps00, local_n0, N1) << std::endl;

        fp = fopen("area_fraction.dat", "a");
        fprintf(fp, "%10d %12.10f\n\n", step, area_fraction);
        fclose(fp);

        // output eta_p data
        if (step % ip.out_freq == 0) {
            frame++;
            output("eta0/"+zeroFill(frame), eta[0], N0, N1, local_n0);
            output("eta1/"+zeroFill(frame), eta[1], N0, N1, local_n0);
            output("eta2/"+zeroFill(frame), eta[2], N0, N1, local_n0);
            output("sig00_00/"+zeroFill(frame), sig00[0][0], N0, N1, local_n0);
            output("sig00_01/"+zeroFill(frame), sig00[0][1], N0, N1, local_n0);
            output("sig00_10/"+zeroFill(frame), sig00[1][0], N0, N1, local_n0);
            output("sig00_11/"+zeroFill(frame), sig00[1][1], N0, N1, local_n0);
            output("dFdw/"+zeroFill(frame), dFdw, N0, N1, local_n0);
            output("bend_2/"+zeroFill(frame), bend_2, N0, N1, local_n0);
            output("ddNuterm0/"+zeroFill(frame), ddNuterm[0], N0, N1, local_n0);
            output("ddNuterm1/"+zeroFill(frame), ddNuterm[1], N0, N1, local_n0);
            output("ddNuterm2/"+zeroFill(frame), ddNuterm[2], N0, N1, local_n0);
            output("ddw0/"+zeroFill(frame), ddw[0], N0, N1, local_n0);
            output("ddw1/"+zeroFill(frame), ddw[1], N0, N1, local_n0);
            output("ddw2/"+zeroFill(frame), ddw[2], N0, N1, local_n0);
            output("inhom_eps00_3/"+zeroFill(frame), inhom_eps00_3[1][1], N0, N1, local_n0);
            output("inhom_eps00_4/"+zeroFill(frame), inhom_eps00_4[1][1], N0, N1, local_n0);
            output("w/"+zeroFill(frame), w, N0, N1, local_n0);
            output("ux/"+zeroFill(frame), ux, N0, N1, local_n0);
            output("uy/"+zeroFill(frame), uy, N0, N1, local_n0);
 
       }
    }

    myfile.close();
    
    fftw_mpi_cleanup();
    MPI_Finalize();

    return 0;
}
