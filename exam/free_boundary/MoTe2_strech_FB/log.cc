

#include "log.h"

void log_greens_function(double *** G, double ** kxy, ptrdiff_t local_n0, ptrdiff_t N1)
{
    FILE * fp = fopen("greens_function.dat", "w");
    for (ptrdiff_t i=0; i<local_n0; i++)
    for (ptrdiff_t j=0; j<N1/2+1; j++)
    {
        int ndx = i*(N1/2+1) + j;
        fprintf(fp, "%10ld", i);
        fprintf(fp, "%10ld", j);
        fprintf(fp, "%10.3f", kxy[0][ndx]);
        fprintf(fp, "%10.3f", kxy[1][ndx]);
        fprintf(fp, "%10.3f", G[0][0][ndx]);
        fprintf(fp, "%10.3f", G[0][1][ndx]);
        fprintf(fp, "%10.3f", G[1][0][ndx]);
        fprintf(fp, "%10.3f", G[1][1][ndx]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void log_elastic_tensors(double **** lam, double **** epsT)
{
    FILE * fp = fopen("elastic_tensors.dat","w");
    const char * f1 = "{%10.4f,%10.4f,%10.4f,%10.4f}\n";
    const char * f2 = "epsT[m=%1d][p=%1d] = {{%8.4f,%8.4f}, {%8.4f,%8.4f}}\n";

    fprintf(fp, "%10s\n", "lambda_ijkl = ");
    fprintf(fp, f1,  lam[0][0][0][0], lam[0][0][1][1], lam[0][0][0][1], lam[0][0][1][0]);
    fprintf(fp, f1,  lam[1][1][0][0], lam[1][1][1][1], lam[1][1][0][1], lam[1][1][1][0]);
    fprintf(fp, f1,  lam[0][1][0][0], lam[0][1][1][1], lam[0][1][0][1], lam[0][1][1][0]);
    fprintf(fp, f1,  lam[1][0][0][0], lam[1][0][1][1], lam[1][0][0][1], lam[1][0][1][0]);

    fprintf(fp, "\n");
    for (int p=0; p<3; p++)
        fprintf(fp, f2, 0, p, epsT[0][p][0][0], epsT[0][p][0][1], epsT[0][p][1][0], epsT[0][p][1][1]);

    fprintf(fp, "\n");
    for (int p=0; p<3; p++)
        fprintf(fp, f2, 1, p, epsT[1][p][0][0], epsT[1][p][0][1], epsT[1][p][1][0], epsT[1][p][1][1]);

    fclose(fp);
}
