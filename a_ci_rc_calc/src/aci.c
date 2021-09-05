#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* [a.u.] energije izrazena u [eV] */
#define eV_au (2 * 13.60569172)
/* Borov radius */
#define CGS_a0 5.2918e-9
/* Brzina svetlosti        */
#define CGS_c 2.9979e10
/* Naelektrisanje elektrona        */
#define CGS_e 4.802e-10
/* Masa elektrona       */
#define CGS_me 9.108e-28
/* Plankova konstanta   */
#define CGS_h 6.625e-27
#define CGS_ch 1.054e-27
/* Konstanta fine strukture     */
#define CGS_alpha (1 / 137.0)
/* Borov magneton       */
#define CGS_mu_0 0.9273e-20
/* Ridbergova konstanta */
#define CGS_R 109737
/* Bolcmanova konstanta */
#define CGS_k 1.380e-16
/* Avogadrov broj       */
#define CGS_N_0 6.02e23
/* Loshmidt-ov broj     */
#define CGS_L_0 2.687e19
/* Komptonovka talasna duzina elektrona */
#define CGS_lambda_hash_e 3.861e-11
/* Klasicni radius elektrona    */
#define CGS_r_0 2.818e-13

double r_n_OK(double Ne)
{
    double a, b;

    a = 1 / Ne; /* cm^3 / jonu */
    b = pow((3 * a / (4 * M_PI)), (1.0 / 3.0));

    return (b);
}

double debaj_OK(double Ne, double T)
{
    double f;

    f = pow((CGS_k * T / (4 * M_PI * Ne * pow(CGS_e, 2))), 0.5);
    return (f);
}

double aci(double xx)
{
    double rv;
    double a0, a1, a2, a3, a4, a5, a6, a7;
    a0 = 1.00147;    //+/- 6.949e-05    (0.006939%)
    a1 = -0.371242;  //+/- 0.001344     (0.3621%)
    a2 = 0.236772;   //+/- 0.008196     (3.462%)
    a3 = 0.183652;   //+/- 0.02203      (12%)
    a4 = -0.25849;   //+/- 0.03011      (11.65%)
    a5 = 0.113241;   //+/- 0.02181      (19.26%)
    a6 = -0.0215963; //+/- 0.00798      (36.95%)
    a7 = 0.00142737; //+/- 0.001159     (81.23%)

    rv = a0;
    rv += a1 * pow(xx, 1.0);
    rv += a2 * pow(xx, 2.0);
    rv += a3 * pow(xx, 3.0);
    rv += a4 * pow(xx, 4.0);
    rv += a5 * pow(xx, 5.0);
    rv += a6 * pow(xx, 6.0);
    rv += a7 * pow(xx, 7.0);

    return (rv);
}

double rc_aci(double r_ws, double r_D)
{
    double rc;
    double x;

    if (r_ws < 0.0)
    {
        return (0.0);
    }

    if (r_D < 0.0)
    {
        return (0.0);
    }

    x = r_ws / r_D;

    if (x < 0)
    {
        return (0.0);
    }
    if (x > 2.0)
    {
        return (0.0);
    }

    rc = r_ws * aci(x);
    return (rc);
}

double rc_od_Ne_i_T(double Ne, double T)
{
    double rD, rs;

    rD = debaj_OK(Ne, T) / CGS_a0;
    rs = r_n_OK(Ne) / CGS_a0;

    /* printf("r_D = %g\n", rD); */
    /* printf("r_n = %g\n", rs); */

    return (rc_aci(rs, rD));
}

int main(int argc, char **argv)
{
    double Ne, T;
    double rd;
    double rws;
    double rc;

    if (argc != 3)
    {
        fprintf(stderr, "Usage:\n%s Ne[cm^{-3}] T[K]\n", argv[0]);
        return(-1);
    }

    Ne = atof(argv[1]);
    T  = atof(argv[2]);
    
    printf("Input:\tNe = %g cm^-3\tT = %g K\n", Ne, T);
    rd = debaj_OK(Ne, T);
    rws = r_n_OK(Ne);

    printf("Output:\tr Deb. = %g a.u.\tr W.S.= %g a.u.\n", 
        rd / CGS_a0, rws / CGS_a0);
    
    rc = rc_od_Ne_i_T(Ne, T);

    printf("\t\tr cut = %g a.u.\n", rc);

    return (0);
}
