/*! \file madau_absorption.c
 *  \brief Function definition for IGM absorption following Madau 1995.
 *
 *  See Madau 1995, ApJ 441, 18
 *  http://adsabs.harvard.edu/abs/1995ApJ...441...18M
 */
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"madau_absorption.h"
/*! \fn double *madau_absorption(double *lambda, int n_lambda, double z);
 *  \brief Wavelength-dependent absorption owing to neutral H in the IGM.
 *
 *   This code is ported from a routine supplied by Alice Shapley, originally
 *   developed by Kurt Adelberger.
 *
 *   Provide an array of wavelengths in Angstrom, with the length of the 
 *   array n_lambda, at a redshift z and return an array of the attenuation 
 *   owing to neutral H in the IGM along the line of sight.
 *
 *   See Madau 1995, ApJ 441, 18
 *   http://adsabs.harvard.edu/abs/1995ApJ...441...18M
 */
double *madau_absorption(double *lambda, int n_lambda, double z)
{
	//lambda is in angstroms
	double *attenuation;

        double lylambda[32];
        double madauA[32];

	attenuation  = (double *) calloc(n_lambda, sizeof(double));

        int nseries = 4;

	for(int i=0;i<n_lambda;i++)
		attenuation[i] = 1.0;

        //ported from k. adelberger via a. shapley

	double l;

	
	double xc,xem;
	double teffline, teffcont, tefftot;
	double lmin, lmax;

        for(int i=2;i<=nseries+1;i++)
		lylambda[i-2]=912.*pow(i,2)/(pow(i,2)-1);

        madauA[0] = 0.0036;
        madauA[1] = 0.0017;
        madauA[2] = 0.0012;
        madauA[3] = 0.00093;

	//printf("here %d\n",n_lambda);
	//fflush(stdout);
        for(int j=0;j<n_lambda;j++)
	{
                //l = lambda[j]*1.0e4;	//(convert to angstroms)
                l = lambda[j]*(1.+z); //leave in angstroms, redshift
                teffline = 0.0;

                for(int i=0;i<nseries;i++)
		{
                        lmin = lylambda[i];
                        lmax = lylambda[i]*(1.0+z);
                        //;print,'series = ',i,' lmin = ',lmin,' lmax = ',lmax
                        if( (l>=lmin) && (l <= lmax))
			{
                                teffline += madauA[i]*exp(3.46*log(l/lylambda[i]));
			}
                }
                xc  = l/912.0;
                xem = 1.0+z;
                if (xc < 1.0) 
			xc=1.0;
                if (xc > xem) 
			xc=xem;
                teffcont = 0.25*xc*xc*xc*(exp(0.46*log(xem)) - exp(0.46*log(xc)));
                teffcont += 9.4*exp(1.5*log(xc))*(exp(0.18*log(xem)) - exp(0.18*log(xc)));
                teffcont += -0.7*xc*xc*xc*(exp(-1.32*log(xc)) - exp(-1.32*log(xem)));
                teffcont += -0.023*(exp(1.68*log(xem)) - exp(1.68*log(xc)));
                tefftot = teffline+teffcont;

		//printf("l%e t%e\n",lambda[j],tefftot);	
		//fflush(stdout);
                attenuation[j] *= exp(-tefftot);
	}

        return attenuation;
}

void print_madau_absorption(double *lambda, int n_lambda, double z)
{
  int i;
  FILE *fp;
  char fname[200];
  sprintf(fname,"madau_absorption.%4.2f.txt",z);
  double *absorption = madau_absorption(lambda, n_lambda, z);
  fp = fopen(fname,"w");
  for(i=0;i<n_lambda;i++)
    fprintf(fp,"%e\t%e\n",lambda[i],absorption[i]);
  fclose(fp);
  free(absorption);
}
