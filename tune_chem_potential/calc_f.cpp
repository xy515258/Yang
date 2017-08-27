
#include <stdio.h>

#include <fstream>
#include <iostream>
#include <cmath>

int main()
{
double a = 0.055;
double b = 0.04;
//double c = b-a+0.0651018; //MoS2
double d = 0.32;
double c = b-a-d-0.00046124; //MoTe2
double mu = 1.;
double nu = 0.24;

double A = -1;
double D = 0;
double B = -1-2*A+D;
double C = 1+A-2*D;

double epsbar[2][2];
epsbar[0][0] = 0.;
epsbar[1][1] = 0.05;
epsbar[0][1] = 0.0;
epsbar[1][0] = 0.0;

double epsbar_star[2][2];
epsbar_star[0][0] = 0.;
epsbar_star[1][1] = 0.05;
epsbar_star[0][1] = 0.0;
epsbar_star[1][0] = 0.0;


double eps0[2][2];

//MoTe2
/*eps0[0][0] = 0.021485;
eps0[0][1] = 0.027855;
eps0[1][0] = 0.027855;
eps0[1][1] = -0.010679;
*/

eps0[0][0] = -0.0267605634;
eps0[0][1] = 0;
eps0[1][0] = 0;
eps0[1][1] = 0.0375670841;

//MoS2
//eps0[0][0] = -0.0012566761;
//eps0[0][1] = 0;
//eps0[1][0] = 0;
//eps0[1][1] = 0.0370034464;

double eta;
double f;
double f_elas;
double f_bulk;
double f_norm;
double df;
double df_elas;
double df_bulk;
double df_norm;
double lam[2][2][2][2];

double E = 2*mu*(1+nu);
    lam[0][0][0][0] = E/(1-nu*nu);
    lam[0][0][1][1] = E*nu/(1-nu*nu);
    lam[0][0][1][0] = 0.0;
    lam[0][0][0][1] = 0.0;
    lam[1][1][0][0] = E*nu/(1-nu*nu);
    lam[1][1][1][1] = E/(1-nu*nu);
    lam[1][1][0][1] = 0;
    lam[1][1][1][0] = 0;
    lam[0][1][0][0] = 0;
    lam[0][1][1][1] = 0;
    lam[0][1][0][1] = mu;
    lam[0][1][1][0] = mu;
    lam[1][0][0][0] = 0;
    lam[1][0][1][1] = 0;
    lam[1][0][0][1] = mu;
    lam[1][0][1][0] = mu;

std::ofstream myfile ("output", std::ofstream::out);

for(int i=0; i<=20000; i++)
{
    eta = -2+i*0.0002;
    f_bulk = 0.5*a*eta*eta - 0.25*b*eta*eta*eta*eta + c/6*eta*eta*eta*eta*eta*eta + d/8*eta*eta*eta*eta*eta*eta*eta*eta;
    df_bulk = a*eta - b*eta*eta*eta + c*eta*eta*eta*eta*eta + d*eta*eta*eta*eta*eta*eta*eta;
  
    f_elas = 0;
    df_elas = 0;
    for(int ii=0; ii<2; ii++)
    for(int jj=0; jj<2; jj++)
    for(int kk=0; kk<2; kk++)
    for(int ll=0; ll<2; ll++)
   { f_elas += 0.5*lam[ii][jj][kk][ll]*eps0[ii][jj]*eps0[kk][ll]*eta*eta*eta*eta - lam[ii][jj][kk][ll]*eps0[ii][jj]*epsbar[kk][ll]*eta*eta + 0.5*lam[ii][jj][kk][ll]*epsbar[ii][jj]*epsbar[kk][ll];
    df_elas += 2*lam[ii][jj][kk][ll]*eps0[ii][jj]*eps0[kk][ll]*eta*eta*eta - 2*lam[ii][jj][kk][ll]*eps0[ii][jj]*epsbar[kk][ll]*eta;
}

    f_norm = 0;
    df_norm = 0;
    for(int ii=0; ii<2; ii++)
    for(int jj=0; jj<2; jj++)
    for(int kk=0; kk<2; kk++)
    for(int ll=0; ll<2; ll++)
{    f_norm += lam[ii][jj][kk][ll]*eps0[kk][ll]*(epsbar[ii][jj]-epsbar_star[ii][jj])*(A*eta*eta + B*eta*eta*eta*eta + C*eta*eta*eta*eta*eta*eta + D*eta*eta*eta*eta*eta*eta*eta*eta);
    df_norm += lam[ii][jj][kk][ll]*eps0[kk][ll]*(epsbar[ii][jj]-epsbar_star[ii][jj])*(-2*eta+4*eta*eta*eta);
}
 
    f = f_bulk + f_elas + f_norm;
    df = df_bulk + df_elas + df_norm;
    myfile << eta << "    " << f << "    " << df <<std::endl;
   
    if(eta==-1||eta==1)
   std::cout<< eta << "    df= " << df << "   df_bulk= "<< df_bulk << "   df_elas= "<< df_elas << "    df_norm= " << df_norm << std::endl;

}

   std::cout<< lam[0][0][0][0] << "	" << lam[0][0][1][1] << "	" << std::endl;

   std::cout<< "a= " << a << "    b= " << b << "   c= " << c << std::endl;   

myfile.close();
return 0;
}
