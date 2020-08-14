#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "multinest.h"
#include "bpass_models.h"
#include <gsl/gsl_interp.h>

#define SMC_DUST


#define A_MAX 100.0
#define A_MIN  0.0
#define F_ESC_MIN 0.0
#define F_ESC_MAX 1.0

#define FIT_F_ESC
#define FIT_LINES
#define FIT_CONTINUUM

double klambda_smc(double lambda)
{
	double lambda_smc[31] = {0.09,0.116,0.119,0.123,0.127,0.131,0.136,0.140,0.145,0.151,0.157,0.163,0.170,0.178,0.186,0.195,0.205,0.216,0.229,0.242,0.258,0.276,0.296,0.370,0.440,0.550,0.660,0.810,1.250,1.650,2.198};
	double k_smc[31] = {9.5,6.992,6.436,6.297,6.074,5.795,5.575,5.272,5.000,4.776,4.472,4.243,4.013,3.866,3.637,3.489,3.293,3.161,2.947,2.661,2.428,2.220,2.000,1.672,1.374,1.000,0.801,0.567,0.131,0.169,0.016};
	double Rv_smc = 2.74; // gordon et al. 2003
	int nlsmc = 31;
	int    ki = gsl_interp_bsearch(lambda_smc,lambda,0,nlsmc);
	double ksmc =  k_smc[ki] + (k_smc[ki+1] - k_smc[ki])*(lambda - lambda_smc[ki])/(lambda_smc[ki+1] - lambda_smc[ki]);
	return ksmc * Rv_smc;
}


double gaussian_func(double x, double sigma, double mu)
{
	double A = 1./sqrt(2.0*M_PI*sigma*sigma);
	double xx = (x-mu)/sigma;
	return A*exp(-0.5*xx*xx);
}


/******************************************** loglikelihood routine ****************************************************/
// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//	 
// Output arguments
// lnew 						= loglikelihood

void LogLikeExample(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{
	double chi = 1.0;
	int i;
	double sigma = 10.0*Cube[0];
	double mu    = 10.0*Cube[1] - 5.0;
	double l = 0;
	double x;
	double y;
	for(i=0;i<n_data;i++)
	{
		y = gaussian_func(data_x[i],sigma,mu);
		x = (y-data_y[i])/data_ye[i];
	}
	Cube[0] = sigma;
	Cube[1] = mu;
	*lnew = l;
}


void LogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{

	//we have the following parameters
	//A     == amplitude
	//f_bin == binarity
	//z_met == metallicity
	//age   == age

	//we perform an interpolation on our model grid
	//and then compare the model SED with the data


	double A = A_MIN + Cube[0]*(A_MAX-A_MIN);								 //SED amplitude
#ifdef FIT_F_ESC
	double f_esc = F_ESC_MIN + Cube[2]*(F_ESC_MAX-F_ESC_MIN);
#endif //FIT_F_ESC
	//limit age to 9.3081373786380386
	//double a_age   = age[0] + (age[n_age-1]-age[0])*Cube[1]; //log10 age
	double age_max = 9.3081373786380386;
	double a_age   = age[0] + (age_max-age[0])*Cube[1]; //log10 age
	double a_f_bin = 0;	//not using currently
	double a_z_met = 0;	//not using currently

	int m;
	double l = 0;
	double x;
	double y;
#ifdef FIT_LINES
	double yl;
#endif //FIT_LINES
#ifdef FIT_CONTINUUM
	double yc;
#endif //FIT_CONTINUUM
	//interpolation along age direction
	int     k_age = gsl_interp_bsearch(age,a_age,0,n_age-2);
	double dx_age = (a_age - age[k_age])/(age[k_age+1] - age[k_age]);

	printf("dx_age %e\n",dx_age);
	exit(0);

	for(m=0;m<n_data;m++)
	{
		//just interpolate along age direction to start
		y = (1.0-dx_age)*sed_model[0][0][k_age][m] + dx_age*sed_model[0][0][k_age+1][m];
#ifdef FIT_LINES
		yl = (1.0-dx_age)*line_model[0][0][k_age][m] + dx_age*line_model[0][0][k_age+1][m];
#endif //FIT_LINES
#ifdef FIT_CONTINUUM
		yc = (1.0-dx_age)*cont_model[0][0][k_age][m] + dx_age*cont_model[0][0][k_age+1][m];
#endif //FIT_CONTINUUM

#ifdef FIT_F_ESC

#ifdef FIT_LINES
		y += (1.0-f_esc)*yl;
#endif //FIT_LINES
#ifdef FIT_CONTINUUM
		y += (1.0-f_esc)*yc;
#endif //FIT_CONTINUUM


		if(m==0)
			y*=f_esc;

#endif //FIT_F_ESC
		x = (A*y-data_y[m])/data_ye[m];
		//if(((data_y[m]>1.0e-10)&&(fabs(data_x[m]-1216)>25.))||(m==0)) //avoid non-detections for the time being, and Lya
		if(((data_y[m]>1.0e-10)&&(fabs(data_x[m]-1215.67)>25.))) //avoid non-detections for the time being, and Lya

			l += -0.5*x*x;
	}
	Cube[0] = A;
	Cube[1] = a_age;
#ifdef FIT_F_ESC
	Cube[2] = f_esc;
#endif //FIT_F_ESC
	*lnew = l;
}

/***********************************************************************************************************************/




/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *logZerr, void *context)
{
	// convert the 2D Fortran arrays to C arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[*nSamples][*nPar + 2];
	for( i = 0; i < *nPar + 2; i++ )
		for( j = 0; j < *nSamples; j++ )
			postdist[j][i] = posterior[0][i * (*nSamples) + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[*nlive][*nPar + 1];
	for( i = 0; i < *nPar + 1; i++ )
		for( j = 0; j < *nlive; j++ )
			pLivePts[j][i] = physLive[0][i * (*nlive) + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
	int i;
	
	// set the MultiNest sampling parameters
	
	
	int mmodal = 0;					// do mode separation?
	
	int ceff = 0;					// run in constant efficiency mode?
	
	int nlive = 1000;				// number of live points
	
	double efr = 1.0;				// set the required efficiency
	
	double tol = 0.5;				// tol, defines the stopping criteria
	
#ifndef FIT_F_ESC
	int ndims = 2;					// dimensionality (no. of free parameters)
	
	int nPar = 2;					// total no. of parameters including free & derived parameters
	
	int nClsPar = 2;				// no. of parameters to do mode separation on
#else //FIT_F_ESC
	int ndims = 3;					// dimensionality (no. of free parameters)
	
	int nPar = 3;					// total no. of parameters including free & derived parameters
	
	int nClsPar = 3;				// no. of parameters to do mode separation on
#endif //FIT_F_ESC
	
	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;
	
	char root[100] = "chains/test_gaussian-";		// root for output files
	
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	
	int fb = 1;					// need feedback on standard output?
	
	int resume = 0;					// resume from a previous job?
	
	int outfile = 1;				// write output files?
	
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	
	double logZero = -DBL_MAX;			// points with loglike < logZero will be ignored by MultiNest
	
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass

	//load the data first
	char fname[200];
	sprintf(fname,"fnudata_total_02_90340");
	if(argc>1)
	{
		sprintf(fname,"%s",argv[1]);
	}

	printf("Data file name = %s.\n",fname);

	//read data
	read_data(fname);

	//initialize the bpass models
	//1) read the bpass files
	//2) based on the source redshift, apply IGM attenuation
	//3) remember only the sampled locations on the SED
	initialize_bpass_models();
	
	// calling MultiNest
	
	run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, 
	logZero, maxiter, LogLike, dumper, context);
}

/***********************************************************************************************************************/
