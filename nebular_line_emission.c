
/*! \file nebular_line_emission.c
 *  \brief Function defintions for routines that provide the nebular line emission
 *  from H, HeI, HeII, and metals. */
#include <math.h>
#include <gsl/gsl_interp.h>
#include "constants.h"
#include "routines.hpp"
#include "nebular_line_emission.h"

/*! \def FIXED_LINE_WIDTH
 *  \brief Fractional line width in delta lambda / lambda,
 *   If set, this is the fractional line width for all lines */
//#define FIXED_LINE_WIDTH 0.05 //in delta l / l
//#define FIXED_LINE_WIDTH 0.01 //in delta l / l
//#define FIXED_LINE_WIDTH 0.05 //in delta l / l
//#define FIXED_LINE_WIDTH 0.001 //in delta l / l
//#define FIXED_LINE_WIDTH 0.0001 //in delta l / l
#define FIXED_LINE_WIDTH 0.00001 //in delta l / l

/*! \var Constants C_lines
 *  \brief An instance of the Constants class for use 
 *  in calculating nebular line emission. */
Constants C_lines;

/*! \fn void add_all_spectra(double F_HI, double Z, double n_p, double n_Hep, double n_Hepp, double T, double *&lambdas, double *&spectrums, int *n_lambda)
 *  \brief Adds all nebular limes to the input spectrum, with a 
 *  strength based relative F_HI = F_H_beta.
 */
void add_all_spectra(double F_HI, double Z, double n_p, double n_Hep, double n_Hepp, double T, double *&lambdas, double *&spectrums, int *n_lambda)
{

	//HI first
	/*for(int i=0;i<*n_lambda;i++)
	{
		printf("test here l %e s %e\n",lambdas[i],spectrums[i]);
		fflush(stdout);
	}*/
	add_HI_spectrum(F_HI,T,lambdas,spectrums,n_lambda);
	//HeI second
	add_HeI_spectrum(F_HI,n_Hep,n_p,T,lambdas,spectrums,n_lambda);
	//HeII second
	add_HeII_spectrum(F_HI,n_Hepp,n_p,T,lambdas,spectrums,n_lambda);
	//metals third
	add_metal_spectrum(F_HI,Z,T,lambdas,spectrums,n_lambda);
}

/*! \fn void add_HI_spectrum(double F_HI, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Adds hydrogen recombination line spectrum, with a strength
 *  	  scaled to F_HI = F_H_beta.
 */
void add_HI_spectrum(double F_HI, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
{
	//add the entire HI spectrum
	//normalized such that Hbeta = F_HI
	//T is in K
	double l;
	double lw;
	double lwf = line_width_factor_H(T);
	double f;
	double f_sun = 3.826e33;

	//Balmer lines
	//H alpha
	l  = lambda_H_alpha();
	lw = l*lwf;
	f  = F_HI*F_H_alpha();
	/*for(int i=0;i<*n_lambda;i++)
	{
		printf("test l %e s %e\n",lambda_lines[i],spectrum_lines[i]);
		fflush(stdout);
	}*/
	//printf("f Halpha = %e, F_HI %e, F_H_alpha %e\n",f,F_HI,F_H_alpha());
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);
	


	//H beta
	l  = lambda_H_beta();
	lw = l*lwf;
	f  = F_HI*F_H_beta();
	//printf("H beta flux:  %lf\n",f/f_sun);
	//printf("f Hbeta = %e, F_HI %e, F_H_beta %e\n",f,F_HI,F_H_beta());
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);


	//H gamma
	l  = lambda_H_gamma();
	lw = l*lwf;
	f  = F_HI*F_H_gamma();
	//printf("f Hgamma = %e, F_HI %e, F_H_gamma %e lambda %e\n",f,F_HI,F_H_gamma(),l);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//H delta
	l  = lambda_H_delta();
	lw = l*lwf;
	f  = F_HI*F_H_delta();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//H 10
	l  = lambda_H_10();
	lw = l*lwf;
	f  = F_HI*F_H_10();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//H 15
	l  = lambda_H_15();
	lw = l*lwf;
	f  = F_HI*F_H_15();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);
	
	//H 20
	l  = lambda_H_20();
	lw = l*lwf;
	f  = F_HI*F_H_20();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);
	

	//Paschen lines
	//P alpha
	l  = lambda_H_Palpha();
	lw = l*lwf;
	f  = F_HI*F_H_Palpha();
	//printf("f Paschen alpha = %e, F_HI %e, F_H_Palpha %e lambda %e\n",f,F_HI,F_H_Palpha(),l);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//P beta
	l  = lambda_H_Pbeta();
	lw = l*lwf;
	f  = F_HI*F_H_Pbeta();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//P gamma
	l  = lambda_H_Pgamma();
	lw = l*lwf;
	f  = F_HI*F_H_Pgamma();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//P delta
	l  = lambda_H_Pdelta();
	lw = l*lwf;
	f  = F_HI*F_H_Pdelta();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//P 10
	l  = lambda_H_P10();
	lw = l*lwf;
	f  = F_HI*F_H_P10();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//P 15
	l  = lambda_H_P15();
	lw = l*lwf;
	f  = F_HI*F_H_P15();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);
	
	//P 20
	l  = lambda_H_P20();
	lw = l*lwf;
	f  = F_HI*F_H_P20();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//Brackett lines
	//B alpha
	l  = lambda_H_Balpha();
	lw = l*lwf;
	f  = F_HI*F_H_Balpha();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//B beta
	l  = lambda_H_Bbeta();
	lw = l*lwf;
	f  = F_HI*F_H_Bbeta();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//B gamma
	l  = lambda_H_Bgamma();
	lw = l*lwf;
	f  = F_HI*F_H_Bgamma();
	//printf("f Brackett gamma = %e, F_HI %e, F_H_Bgamma %e lambda %e\n",f,F_HI,F_H_Bgamma(),l);

	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//B delta
	l  = lambda_H_Bdelta();
	lw = l*lwf;
	f  = F_HI*F_H_Bdelta();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//B 10
	l  = lambda_H_B10();
	lw = l*lwf;
	f  = F_HI*F_H_B10();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//B 15
	l  = lambda_H_B15();
	lw = l*lwf;
	f  = F_HI*F_H_B15();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);
	
	//B 20
	l  = lambda_H_B20();
	lw = l*lwf;
	f  = F_HI*F_H_B20();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);


}

/*! \fn void add_HeI_spectrum(double F_Hbeta, double n_Hep, double n_p, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Adds HeI recombination line spectrum, with a strength
 *  	  scaled to F_HI = F_H_beta.
 */
void add_HeI_spectrum(double F_Hbeta, double n_Hep, double n_p, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
{
	//add the entire HeI spectrum
	//normalized such that Hbeta = F_Hbeta
	//T is in K
	double l;
	double lw;
	double lwf = line_width_factor_He(T);
	double f;
	double F_HeI = F_Hbeta*F_HeI_norm(n_Hep,n_p);

	//HeI 4471
	l  = lambda_HeI_4471();
	lw = l*lwf;
	f  = F_HeI*F_HeI_4471();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	
	//HeI 5876
	l  = lambda_HeI_5876();
	lw = l*lwf;
	f  = F_HeI*F_HeI_5876();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 4026
	l  = lambda_HeI_4026();
	lw = l*lwf;
	f  = F_HeI*F_HeI_4026();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 7065
	l  = lambda_HeI_7065();
	lw = l*lwf;
	f  = F_HeI*F_HeI_7065();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 10830
	l  = lambda_HeI_10830();
	lw = l*lwf;
	f  = F_HeI*F_HeI_10830();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 3889
	l  = lambda_HeI_3889();
	lw = l*lwf;
	f  = F_HeI*F_HeI_3889();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 3187
	l  = lambda_HeI_3187();
	lw = l*lwf;
	f  = F_HeI*F_HeI_3187();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 6678
	l  = lambda_HeI_6678();
	lw = l*lwf;
	f  = F_HeI*F_HeI_6678();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 4922
	l  = lambda_HeI_4922();
	lw = l*lwf;
	f  = F_HeI*F_HeI_4922();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 5016
	l  = lambda_HeI_5016();
	lw = l*lwf;
	f  = F_HeI*F_HeI_5016();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeI 3965
	l  = lambda_HeI_3965();
	lw = l*lwf;
	f  = F_HeI*F_HeI_3965();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	
}

/*! \fn void add_HeII_spectrum(double F_Hbeta, double n_Hepp, double n_p, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Adds HeII recombination line spectrum, with a strength
 *  	  scaled to F_HI = F_H_beta.
 */
void add_HeII_spectrum(double F_Hbeta, double n_Hepp, double n_p, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
{
	//add the entire HeII spectrum
	//normalized such that Hbeta = F_Hbeta
	//T is in K
	double l;
	double lw;
	double lwf = line_width_factor_He(T);
	double f;
	double F_HeII = F_Hbeta*F_HeII_norm(n_Hepp,n_p);
	//printf("F_HEII %e\n",F_HeII);
	//exit(-1);

	//return if evaluation is not necessary
	if(n_Hepp<1.0e-6*n_p)
		return;


	//HeII 32
	l  = lambda_HeII_32();
	lw = l*lwf;
	f  = F_HeII*F_HeII_32();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 42
	l  = lambda_HeII_42();
	lw = l*lwf;
	f  = F_HeII*F_HeII_42();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 52
	l  = lambda_HeII_52();
	lw = l*lwf;
	f  = F_HeII*F_HeII_52();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 62
	l  = lambda_HeII_62();
	lw = l*lwf;
	f  = F_HeII*F_HeII_62();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 102
	l  = lambda_HeII_102();
	lw = l*lwf;
	f  = F_HeII*F_HeII_102();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 152
	l  = lambda_HeII_152();
	lw = l*lwf;
	f  = F_HeII*F_HeII_152();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 202
	l  = lambda_HeII_202();
	lw = l*lwf;
	f  = F_HeII*F_HeII_202();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 43
	l  = lambda_HeII_43();
	lw = l*lwf;
	f  = F_HeII*F_HeII_43();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 53
	l  = lambda_HeII_53();
	lw = l*lwf;
	f  = F_HeII*F_HeII_53();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 63
	l  = lambda_HeII_63();
	lw = l*lwf;
	f  = F_HeII*F_HeII_63();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 73
	l  = lambda_HeII_73();
	lw = l*lwf;
	f  = F_HeII*F_HeII_73();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 103
	l  = lambda_HeII_103();
	lw = l*lwf;
	f  = F_HeII*F_HeII_103();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 153
	l  = lambda_HeII_153();
	lw = l*lwf;
	f  = F_HeII*F_HeII_153();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 203
	l  = lambda_HeII_203();
	lw = l*lwf;
	f  = F_HeII*F_HeII_203();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 54
	l  = lambda_HeII_54();
	lw = l*lwf;
	f  = F_HeII*F_HeII_54();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//printf("l %e\n",lambda_HeII_64());
	//exit(-1);

	//HeII 64 //causes problems with h alpha currently
	/*l  = lambda_HeII_64();
	lw = l*lwf;
	f  = F_HeII*F_HeII_64();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);*/


	//HeII 74
	l  = lambda_HeII_74();
	lw = l*lwf;
	f  = F_HeII*F_HeII_74();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 84
	l  = lambda_HeII_84();
	lw = l*lwf;
	f  = F_HeII*F_HeII_84();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);


	//causes problems with the spectrum
	//HeII 104
	l  = lambda_HeII_104();
	lw = l*lwf;
	f  = F_HeII*F_HeII_104();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 154
	l  = lambda_HeII_154();
	lw = l*lwf;
	f  = F_HeII*F_HeII_154();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//HeII 204
  	l  = lambda_HeII_204();
	lw = l*lwf;
	f  = F_HeII*F_HeII_204();
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

}

/*! \fn void add_metal_spectrum(double F_Hbeta, double Z, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Adds metal recombination line spectrum, with a strength
 *  	  scaled to F_HI = F_H_beta, according to 
 *
 *        See P. Anders & U. Fritze-v. Alvensleben 2003, AA, 401, 1063-1070. 
 *
 *        http://adsabs.harvard.edu/abs/2003A%26A...401.1063A
 *	  
 */
void add_metal_spectrum(double F_Hbeta, double Z, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
{
	//add the entire metal line spectrum
	//normalized such that Hbeta = F_Hbeta
	//T is in K
	double l;
	double lw;
	double lwf_He = line_width_factor_He(T);
	double lwf_Ar = line_width_factor_Ar(T);
	double lwf_C  = line_width_factor_C(T);
	double lwf_O  = line_width_factor_O(T);
	double lwf_N  = line_width_factor_N(T);
	double lwf_Ne = line_width_factor_Ne(T);
	double lwf_Mg = line_width_factor_Mg(T);
	double lwf_S  = line_width_factor_S(T);
	double f;

	//printf("F %e\n",F_Hbeta);
	//fflush(stdout);


	//H_epsilon + NeIII 3970 - H_e not added in hydrogen spectrum
	l  = lambda_H_epsilon_HeII_3970();
	lw = l*lwf_Ar;
	f = F_Hbeta * F_H_epsilon_HeII_3970(Z);
	//printf("l %e f %e\n",l,f);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//ArIII 7136
	l  = lambda_ArIII_7136();
	lw = l*lwf_Ar;
	f = F_Hbeta * F_ArIII_7136(Z);
	//printf("l %e f %e\n",l,f);

	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//ArIII 7751
	l  = lambda_ArIII_7751();
	lw = l*lwf_Ar;
	f = F_Hbeta * F_ArIII_7751(Z);
	//printf("l %e f %e\n",l,f);

	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//ArIV HeI 4711
	//we already have HeI 4711 from the helium spectrum
	//l  = lambda_ArIV_HeI_4711();
	//lw = l*lwf_He; //make it large
	//f = F_Hbeta * F_ArIV_HeI_4711(Z);
	//printf("l %e f %e Z %e\n",l,f,Z);
	//add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//CII 1335
	l  = lambda_CII_1335();
	lw = l*lwf_C; //make it large
	f = F_Hbeta * F_CII_1335(Z);
	//printf("l %e f %e Z %e\n",l,f,Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//CII 2326
	l  = lambda_CII_2326();
	lw = l*lwf_C; //make it large
	f = F_Hbeta * F_CII_2326(Z);
	//printf("l %e f %e\n",l,f);

	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//CIII 1909
	l  = lambda_CIII_1909();
	lw = l*lwf_C; //make it large
	f = F_Hbeta * F_CIII_1909(Z);
	//printf("l %e f %e\n",l,f);

	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//MgII 2798
	l  = lambda_MgII_2798();
	lw = l*lwf_Mg; //make it large
	f = F_Hbeta * F_MgII_2798(Z);
	/*printf("l %e lw %e f %e MgII %e\n",l,lw,f,F_MgII_2798(Z));

	FILE *fpt = fopen("line_test.txt","w");
	for(int i=0;i<*n_lambda;i++)
	{
		fprintf(fpt,"%10.9e\t%e\n",lambda_lines[i],spectrum_lines[i]);
		if(i>0)
			if(lambda_lines[i]<=lambda_lines[i-1])
				printf("ll %e %e i %d\n",lambda_lines[i-1],lambda_lines[i],i);
	}
	fclose(fpt);*/
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//NI 5199
	l  = lambda_NI_5199();
	lw = l*lwf_N; //make it large
	f = F_Hbeta * F_NI_5199(Z);
	//printf("F_NI_5199 %e\n",F_NI_5199(Z));
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//NII 2141
	l  = lambda_NII_2141();
	lw = l*lwf_N; //make it large
	f = F_Hbeta * F_NII_2141(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//NII 5755
	l  = lambda_NII_5755();
	lw = l*lwf_N; //make it large
	f = F_Hbeta * F_NII_5755(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//NII 6548
	l  = lambda_NII_6548();
	lw = l*lwf_N; //make it large
	f = F_Hbeta * F_NII_6548(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//NII 6583
	l  = lambda_NII_6583();
	lw = l*lwf_N; //make it large
	f = F_Hbeta * F_NII_6583(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//NeIII 5755
	l  = lambda_NeIII_3869();
	lw = l*lwf_Ne; //make it large
	f = F_Hbeta * F_NeIII_3869(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//OII 3727
	l  = lambda_OII_3727();
	lw = l*lwf_O; //make it large
	f = F_Hbeta * F_OII_3727(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//OIII 1663
	l  = lambda_OIII_1663();
	lw = l*lwf_O; //make it large
	f = F_Hbeta * F_OIII_1663(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//OIII 4363
	l  = lambda_OIII_4363();
	lw = l*lwf_O; //make it large
	f = F_Hbeta * F_OIII_4363(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//OIII 4959
	l  = lambda_OIII_4959();
	lw = l*lwf_O; //make it large
	f = F_Hbeta * F_OIII_4959(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//OIII 5007
	l  = lambda_OIII_5007();
	lw = l*lwf_O; //make it large
	f = F_Hbeta * F_OIII_5007(Z);
	//printf("f H OIII = %e, F_HI %e, F_OIII %e lambda %e Z %e\n",f,F_Hbeta,F_OIII_5007(Z),l,Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//OI 6300
	l  = lambda_OI_6300();
	lw = l*lwf_O; //make it large
	f = F_Hbeta * F_OI_6300(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//OII 7320
	l  = lambda_OII_7320();
	lw = l*lwf_O; //make it large
	f = F_Hbeta * F_OII_7320(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//OII 7331
	l  = lambda_OII_7331();
	lw = l*lwf_O; //make it large
	f = F_Hbeta * F_OII_7331(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SII 4069
	l  = lambda_SII_4069();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SII_4069(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SII 4076
	l  = lambda_SII_4076();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SII_4076(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SII 6312
	l  = lambda_SIII_6312();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SIII_6312(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SII 6716
	l  = lambda_SII_6716();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SII_6716(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SII 6730
	l  = lambda_SII_6730();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SII_6730(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SIII 9069
	l  = lambda_SIII_9069();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SIII_9069(Z);
	//printf("F_SIII_9069 %e\n",F_SIII_9069(Z));
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SIII 9531
	l  = lambda_SIII_9531();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SIII_9531(Z);
	//printf("F_SIII_9531 %e\n",F_SIII_9531(Z));
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SII 10287
	l  = lambda_SII_10287();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SII_10287(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SII 10320
	l  = lambda_SII_10320();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SII_10320(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

	//SII 10336
	l  = lambda_SII_10336();
	lw = l*lwf_S; //make it large
	f = F_Hbeta * F_SII_10336(Z);
	add_line_to_spectrum(l,lw,f,lambda_lines,spectrum_lines,n_lambda);

}

/*! \fn double line_profile(double lambda, double lambda_c, double lambda_width)
 *  \brief Gaussian line profile for a naturally broadened line. */
double line_profile(double lambda, double lambda_c, double lambda_width)
{
	double sigma = lambda_width;
	double mu    = lambda_c;
	double x     = lambda;

	return exp(-pow(x-mu,2)/(2*sigma*sigma))/sqrt(2*M_PI*sigma*sigma);
	
}


/*! \fn double line_profile_nu(double nu, double nu_c, double nu_width)
 *  \brief Gaussian line profile for a naturally broadened line. */
double line_profile_nu(double nu, double nu_c, double nu_width)
{
	double sigma = nu_width;
	double mu    = nu_c;
	double x     = nu;

	return exp(-pow(x-mu,2)/(2*sigma*sigma))/sqrt(2*M_PI*sigma*sigma);
	
}

/*! \fn void add_line_to_spectrum(double lambda, double lambda_width, double f, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Insers a single line into the spectrum. */
void add_line_to_spectrum(double lambda, double lambda_width, double f, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
{
	int n_profile = 250;
	double *lambda_profile   = calloc_double_array(n_profile);
	double *spectrum_profile = calloc_double_array(n_profile);

	double c = C_lines.c*1.0e6;
	int ioff = 0;
	int no = (*n_lambda);
	double nw = 5;
	size_t il = gsl_interp_bsearch(lambda_lines,lambda-nw*lambda_width,0,no-1);	//find where the line starts;
	int ip,ix,in;	
	double lambda_change = 1.0e-4;//angs
	double *lambda_save;
	double *spectrum_save;
	double *lambda_new;
	double *spectrum_new;
	double lp;
	double tol =1.0e-60;
	double nu;

	double lambda_max=0;
	double f_nu_max=0;

	double flux_sum=0;
	double nu_o;
	double line_width_nu;
	double nu_c;

	double ltol = 1.0e-12;

	gsl_spline *spline_spectrum_profile;
	gsl_interp *interp_spectrum_profile;
	gsl_interp_accel *acc_spectrum_profile;

	//printf("%d\n",*n_lambda);

	
	line_width_nu = c/lambda - c/(lambda + lambda_width);
	nu_c = c/lambda;
	//printf("%e\t%e\n",nu_c,line_width_nu);


	//limit the extent of where we add lines
	if((lambda-nw*lambda_width<lambda_lines[0])||(lambda+nw*lambda_width>lambda_lines[(*n_lambda)-1]))
	{
		//don't extend the spectrum unecessarily
		free(lambda_profile);
		free(spectrum_profile);
		return;
	}
	
	
	//make the blue half of the line profile 
	ioff = 0;
	for(int i=0;i<(n_profile-1)/2+1;i++)
	{
		lambda_profile[i]   = double_linear_index(i,(n_profile-1)/2+1,-nw*lambda_width,0) + lambda;
	}
	//make the red half of the line profile
	ioff = (n_profile-1)/2;
	for(int i=0;i<(n_profile-1)/2+2;i++)
	{
		lambda_profile[i+ioff]   = double_linear_index(i,(n_profile-1)/2+2,0,nw*lambda_width) + lambda;
	}

	for(int i=0;i<n_profile-1;i++)
	{
		if(lambda_profile[i]==lambda_profile[i+1])
		{
			printf("PROFILE ERROR i %d lp %e %e\n",i,lambda_profile[i],lambda_profile[i+1]);
		}
		//printf("i %d lp %e\n",i,lambda_profile[i]);

	}
	//the profile always seems OK


	//exit(0);


	//store the input spectra

	lambda_save   = calloc_double_array(no);
	spectrum_save = calloc_double_array(no);

	for(int i=0;i<no;i++)
	{
		lambda_save[i] = lambda_lines[i];
		spectrum_save[i] = spectrum_lines[i];

		//if(lambda_save[i]<0.2)
		//	printf("i %d ls %10.9e\n",i,lambda_save[i]);

		if(i>0)
			if(fabs(lambda_save[i]-lambda_save[i-1])<ltol)
			{
				printf("A failure\n");
				printf("i %d lln %e llnimo %e\n",i,lambda_save[i],lambda_save[i-1]);
				printf("i %d sln %e slnimo %e\n",i,spectrum_save[i],spectrum_save[i-1]);
				exit(0);

			}
	}

	//printf("test1\n");
	//allocate the new line+spectra
	lambda_new = calloc_double_array(no+n_profile);
	spectrum_new = calloc_double_array(no+n_profile);

	//save
	//old wavelengths blueward of the line profile
	for(int i=0;i<=il;i++)
	{
		lambda_new [i]     = lambda_save[i];

		if(i>0)
			if(fabs(lambda_new[i]-lambda_new[i-1])<ltol)
			{
				printf("D failure\n");
				printf("i %d lln %e llnimo %e\n",i,lambda_new[i],lambda_new[i-1]);
			}
	}
	int np_lim = no+n_profile;

	//printf("test2\n");
	//if the line is simply inserted between two points
	//then go ahead and insert it
	if(lambda_save[il+1]>lambda_profile[n_profile-1])
	{
		//in this case, we can just insert the new spectrum
		for(int i=0;i<n_profile;i++)
		{
			lambda_new[il+1+i] = lambda_profile[i];
		}
		for(int i=0;i<no-il;i++)
			lambda_new[il+i+1+n_profile] = lambda_save[il+1+i];
	}else{

		//in this case the new profile extends across previous
		//parts of the table; now we interleave the profiles

		ip = il+1;	//index for the new wavelength array
		ix = il+1;	//index for the saved wavelength array
		in = 0;		//index for the profile

		//
		while(ip<no+n_profile)
		{
			//if this entry in the line profile comes
			//first, go ahead and add it to the new
			//spectrum

			if(lambda_save[ix]-lambda_profile[in]>ltol)
			{
				//in this case the next wavelength sample is from the
				//line profile being added to the new wavelength array

				//go ahead an add the wavelength sample from the line profile
				//to the new wavelength array, and advance the line profile
				//counter in	

				lambda_new[ip] = lambda_profile[in];
				in++;
			}else{

				//in this case, either...

				if( fabs(lambda_profile[in]-lambda_save[ix])<ltol)
				{

					//1) the next lambda sample in line profile and the
					//previous wavelength array are equal, so it doesn't
					//matter which one we pick -- go ahead and take the
					//saved lambda sample
					in++;

					lambda_new[ip] = lambda_save[ix];
					ix++;

				}else{

					//2) in this case, the next profile sample is at a longer
					//wavelength than the next saved lambda sample

					//there are two further possibilities -- either the 
					//next saved lambda sample is greater than the current
					//new lambda sample or it's not.  if it is, we can just continue
					//as normal, adopting the saved lambda sample as the
					//next new sample
					if(ip==0)
					{

						//first wavelength sample, so go ahead
						lambda_new[ip] = lambda_save[ix];
						ix++;

					}else{

						if(lambda_save[ix]-lambda_new[ip-1] > ltol)
						{
							//OK, the next saved lambda sample is 
							//larger than the last sample in the new array
							//-- go ahead and proceed as normal
							lambda_new[ip] = lambda_save[ix];
							ix++;

						}else{
							//In this case lambda_save <= lambda_new[ip-1]
							//so we don't want to use it.

							//we advance to the next lambda_save sample
							ix++;

							//take the shorter wavelength of 
							//the line profile lambda or the saved lambda
							if(lambda_save[ix]-lambda_profile[in]>ltol)
							{
								lambda_new[ip] = lambda_profile[in];
								in++;
							}else{
								lambda_new[ip] = lambda_save[ix];
								ix++;
							}
						}
					}
				}
			}

			ip++;

			if(in==n_profile)
			{
				//entire profile has been included
				//now append the rest of the table
				//printf("finishing\n");
				//fflush(stdout);
				while((ip<no+n_profile)&&(ix<no))
				{
					lambda_new[ip] = lambda_save[ix];
					//printf("ip %d no+n_p %d ln %e, ix %d no %d ls %e %e\n",ip,no+n_profile,lambda_new[ip],ix,no,lambda_save[ix-1],lambda_save[ix]);

					ix++;
					ip++;
				}
				np_lim = ip;
				break;
			}

		}
		/*
		//in this case the new profile extends across previous
		//parts of the table; now we interleave the profiles

		ip = il+1;	
		ix = il+1;
		in = 0;
		while(ip<no+n_profile)
		{

			//if this entry in the line profile comes
			//first, go ahead and add it to the new
			//spectrum
			if(lambda_profile[in]<lambda_save[ix])
			{
				lambda_new[ip] = lambda_profile[in];
				in++;
			}else{
				//if there is a duplicate entry
				//just add another point in the table
				//if(lambda_profile[in]==lambda_save[ix])
				if( fabs(lambda_profile[in]-lambda_save[ix])<1.0e-9 )
				{
					printf("CASE E in %d no %d n_profile %d lp %e ls %e lpinpo %e\n",in,no,n_profile,lambda_profile[in],lambda_save[ix],lambda_profile[in+1]);
					//move lambda profile
					lambda_profile[in] = 0.5*(lambda_profile[in+1]+lambda_save[ix]);
				}

				//add the saved wavelength point to the table
				lambda_new[ip] = lambda_save[ix];
				ix++;
			}
			
			ip++;

			if(in==n_profile)
			{
				//entire profile has been included
				//now append the rest of the table
				//printf("finishing\n");
				//fflush(stdout);
				while(ip<no+n_profile)
				{
					lambda_new[ip] = lambda_save[ix];
					ix++;
					ip++;
				}
				break;
			}
		}
		*/
	}
	//printf("test3\n");
	//add new wavelengths	

	//allocate linear interpolation 
	interp_spectrum_profile = gsl_interp_alloc(gsl_interp_linear,no);
	acc_spectrum_profile    = gsl_interp_accel_alloc();

	for(int i=0;i<no;i++)
	{
	/*	if(lambda_save[i]>1.0 && lambda_max<1.0)
		{
			printf("lambda_max %e f_nu_max %e\n",lambda_max,f_nu_max);
			fflush(stdout);
		}
		if(lambda_save[i]>lambda_max)
		{
			lambda_max = lambda_save[i];
			f_nu_max  = spectrum_save[i];
		}*/
		lambda_save[i] = log10(lambda_save[i]);
		if(spectrum_save[i]>0)
		{
			spectrum_save[i] = log10(spectrum_save[i]);
		}else{
			spectrum_save[i] = -99.0;
		}
		if(i>0)
			//if(lambda_save[i]<=lambda_save[i-1])
			if(fabs(lambda_save[i]-lambda_save[i-1])<ltol)
			{
				printf("B failure\n");
				printf("i %d lln %e llnimo %e\n",i,pow(10,lambda_save[i]),pow(10,lambda_save[i-1]));
				printf("i %d sln %e slnimo %e\n",i,spectrum_save[i],spectrum_save[i-1]);

			}
	}
	//printf("lambda_max %e f_nu_max %e\n",lambda_max,f_nu_max);
	//fflush(stdout);
	//exit(-1);

	for(int i=0;i<no+n_profile;i++)
	{
		if(lambda_new[i]>lambda_max)
		{
			lambda_max = lambda_new[i];
			f_nu_max  = spectrum_new[i];
		}
	}
	/*printf("lambda_max %e f_nu_max %e\n",lambda_max,f_nu_max);
	fflush(stdout);
	exit(-1);*/
	//printf("test4\n");
	gsl_interp_init(interp_spectrum_profile,lambda_save,spectrum_save,no);

	//interpolate the previous spectrum
	//for(int i=0;i<no;i++)
//	for(int i=0;i<no+n_profile;i++)
	for(int i=0;i<np_lim;i++)
	{
		//spectrum_new[i] = gsl_interp_eval(interp_spectrum_profile,lambda_save,spectrum_save,lambda_new[i],acc_spectrum_profile);

		spectrum_new[i] = gsl_interp_eval(interp_spectrum_profile,lambda_save,spectrum_save,log10(lambda_new[i]),acc_spectrum_profile);
		spectrum_new[i] = pow(10.0,spectrum_new[i]);
	}

	//add the new line emission

	// f is in erg/s

//	for(int i=0;i<no+n_profile;i++)
	for(int i=0;i<np_lim;i++)
	{
		nu = c/lambda_new[i];
		//lp = f*line_profile(lambda_new[i],lambda,lambda_width)/nu;

		lp = f * line_profile_nu(nu,nu_c,line_width_nu);    // Normalizes

		// Make sure flux sum = 1

		if(lp>0)
		{
			//printf("%d",i);
			//flux_sum = flux_sum + lp * (c/lambda_new[i-1] - c/lambda_new[i-1]);

			//flux_sum = flux_sum + line_profile(lambda_new[i],lambda,lambda_width)*(lambda_new[i+1]-lambda_new[i]);
		
			nu_o = c/lambda_new[i-1];
			flux_sum = flux_sum + lp * (nu_o - nu);
			//flux_sum = line_profile_nu(nu,nu_c,line_width_nu);
			//printf("\n%e \t %e \t %e \n",nu,nu_c,line_width_nu);
			//printf("%e\t%e\n",lp,spectrum_new[i]);

		}


		//multiples line_profile by total flux - something wrong with profile?

		if(lp<0)
		{
			printf("l %e lc %e lw %e f %e lp %e\n",lambda_new[i],lambda,lambda_width,f,lp);
			exit(-1);
		}

		spectrum_new[i] += lp;

		//ensure no spline funny business
		if(fabs(spectrum_new[i])<tol)
			spectrum_new[i] = 0;
		
		if(i>0)
			if(fabs(lambda_new[i]-lambda_new[i-1])<ltol)
			{
				printf("F failure\n");
				printf("i %d lln %e llnimo %e\n",i,lambda_new[i],lambda_new[i-1]);
			}
	}

	//printf("Total line flux:   %e\t%e\n",flux_sum/3.826e33,flux_sum/f);
	//printf("Total line flux recovered:   f %e %e %e\n",f, flux_sum, flux_sum/f);
	if((f>0)&&(fabs(flux_sum/f-1.0)>1.0e-2))
	{
		printf("Error --- inaccurate line flux.\n");
		printf("Total line flux recovered:   f %e %e %e\n",f, flux_sum, flux_sum/f);
		exit(-1);
	}


	//printf("test5\n");

/*
	for(int i=0;i<no+n_profile-1;i++)
	{
		if(lambda_new[i+1]<lambda_new[i])
		{
			lambda_new[i+1]   = 0.5*(lambda_new[i+2] + lambda_new[i]);
			spectrum_new[i+1] = 0.5*(spectrum_new[i+2] + spectrum_new[i]);

		}
	}*/
		//printf("l %e\n",lambda_lines[i]);



	//adjust n_lambda
	*n_lambda = no+n_profile;
	*n_lambda = np_lim;


	free(lambda_lines);
	free(spectrum_lines);

	lambda_lines   = NULL;
	spectrum_lines = NULL;

	lambda_lines   = calloc_double_array(*n_lambda);
	spectrum_lines = calloc_double_array(*n_lambda);
	lp = 0;
	for(int i=0;i<*n_lambda;i++)
	{
	
		lambda_lines[i]   = lambda_new[i];	
		spectrum_lines[i] = spectrum_new[i];	
		if(fabs(lambda_lines[i]-lp)<ltol)
		{
			printf("error here i %d lp %e ll %e\n",i,lp,lambda_lines[i]);
			printf("lambda %e\n",lambda);
			fflush(stdout);
			exit(-1);
		}
		lp = lambda_lines[i];
	}

	//printf("test6\n");

	//gsl_interp_free(interp_spectrum_profile);
	//gsl_interp_accel_free(acc_spectrum_profile);
	//free(lambda_save);
	//free(spectrum_save);
	//free(lambda_new);
	//free(spectrum_new);
	//free(lambda_profile);
	//free(spectrum_profile);
	
	//return;//maybe need to comment this out?
	

	//gsl_spline_free(spline_spectrum_profile);
	gsl_interp_free(interp_spectrum_profile);
	gsl_interp_accel_free(acc_spectrum_profile);
	free(lambda_save);
	free(spectrum_save);
	free(lambda_new);
	free(spectrum_new);
	free(lambda_profile);
	free(spectrum_profile);

	return;
}

/*! \fn double *line_spectrum(int *n_lambda, double lambda_min, double lambda_max, double *&lambda_lines)
 *  \brief Allocates an empty array to use as a line spectrum. */
double *line_spectrum(int *n_lambda, double lambda_min, double lambda_max, double *&lambda_lines)
{
	double *spectrum;
	lambda_lines = calloc_double_array(*n_lambda);
	for(int i=0;i<(*n_lambda);i++)
		lambda_lines[i] = double_log10_index(i,(*n_lambda),lambda_min,lambda_max);
	spectrum = calloc_double_array((*n_lambda));
	for(int i=0;i<(*n_lambda);i++)
		spectrum[i] = 0;

	return spectrum;
}
/*! \fn double line_width_factor(double m, double T)
 *  \brief Returns the fractional line width in delta_nu/nu.  Either set by
 *  FIXED_LINE_WIDTH or the mass of the ion. */
double line_width_factor(double m, double T)
{
	//returns line width in delta_nu/nu 
	double k = C_lines.k_b_cgs;
	double c = C_lines.c_cgs;

#ifndef FIXED_LINE_WIDTH
	//return sqrt(8*log(2)*k/(c*c))*sqrt(T/m);
	return sqrt(k*T/(m*c*c));
#else //FIXED_LINE_WIDTH
	return (double) FIXED_LINE_WIDTH;
#endif //FIXED_LINE_WIDTH
}

/*! \fn double line_width_factor_H(double T)
 *  \brief Thermal line width for H */
double line_width_factor_H(double T)
{
	//returns line width in delta_nu/nu 
	double m = C_lines.m_proton_cgs;
	return line_width_factor(m,T);
}
/*! \fn double line_width_factor_He(double T)
 *  \brief Thermal line width for He */
double line_width_factor_He(double T)
{
	//returns line width in delta_nu/nu 
	double m = C_lines.m_alpha_cgs;
	return line_width_factor(m,T);
}

/*! \fn double line_width_factor_Ar(double T)
 *  \brief Thermal line width for Ar */
double line_width_factor_Ar(double T)
{
	//returns line width in delta_nu/nu 
	double m = 39.948 * C_lines.amu_to_cgs;
	return line_width_factor(m,T);
}
/*! \fn double line_width_factor_C(double T)
 *  \brief Thermal line width for C */
double line_width_factor_C(double T)
{
	//returns line width in delta_nu/nu 
	double m = 12.0956 * C_lines.amu_to_cgs;
	return line_width_factor(m,T);
}
/*! \fn double line_width_factor_Mg(double T)
 *  \brief Thermal line width for Mg */
double line_width_factor_Mg(double T)
{
	//returns line width in delta_nu/nu 
	double m = 24.3050 * C_lines.amu_to_cgs;
	return line_width_factor(m,T);
}

/*! \fn double line_width_factor_N(double T)
 *  \brief Thermal line width for N */
double line_width_factor_N(double T)
{
	//returns line width in delta_nu/nu 
	double m = 14.003241 * C_lines.amu_to_cgs;
	return line_width_factor(m,T);
}

/*! \fn double line_width_factor_Ne(double T)
 *  \brief Thermal line width for Ne */
double line_width_factor_Ne(double T)
{
	//returns line width in delta_nu/nu 
	double m = 20.1797 * C_lines.amu_to_cgs;
	return line_width_factor(m,T);
}

/*! \fn double line_width_factor_O(double T)
 *  \brief Thermal line width for O */
double line_width_factor_O(double T)
{
	//returns line width in delta_nu/nu 
	double m = 16.0 * C_lines.amu_to_cgs;
	return line_width_factor(m,T);
}


/*! \fn double line_width_factor_S(double T)
 *  \brief Thermal line width for S */
double line_width_factor_S(double T)
{
	//returns line width in delta_nu/nu 
	double m = 28.96880* C_lines.amu_to_cgs;
	return line_width_factor(m,T);
}

/*! \fn double lambda_H(int n1, int n2)
 *  \brief Wavelength of a transition in the H atom between
 *  orbital states n1 and n2. The Rydbeg formula.*/
double lambda_H(int n1, int n2)
{
	//double lambda_inv = C_lines.rydberg_cgs*(1.0/pow(n1,2) - 1.0/(pow(n2,2)));
	double lambda_inv = C_lines.rydberg_cgs*(1.0/n1/n1 - 1.0/n2/n2);
	
	return (1./lambda_inv)*1.0e6/1.0e2;
}
/*! \fn double lambda_He(int n1, int n2)
 *  \brief Wavelength of a transition in the He atom between
 *  orbital states n1 and n2. The Z-dependent Rydbeg formula.*/
double lambda_HeII(int n1, int n2)
{
	double Z = 2.0;
	double lambda_inv = Z*Z*C_lines.rydberg_cgs*(1./pow(n1,2) - 1.0/(pow(n2,2)));

	

	return (1./lambda_inv)*1.0e6/1.0e2;
}

/*  Assumes HI recombination lines at ne=100, 
 *  T=10^4, Case B recombination.
 *  Table 4.4 of Osterbrock & Ferland */


/*! \fn double lambda_H_alpha(void)
 *  \brief Wavelength of H-alpha */
double lambda_H_alpha(void)
{
	//in microns
	return lambda_H(2,3);
}

/*! \fn double lambda_H_beta(void)
 *  \brief Wavelength of H-beta */
double lambda_H_beta(void)
{
	//in microns
	return lambda_H(2,4);
}

/*! \fn double lambda_H_gamma(void)
 *  \brief Wavelength of H-gamma*/
double lambda_H_gamma(void)
{
	//in microns
	return lambda_H(2,5);
}

/*! \fn double lambda_H_delta(void)
 *  \brief Wavelength of H-delta*/
double lambda_H_delta(void)
{
	//in microns
	return lambda_H(2,6);
}
/*! \fn double lambda_H_10(void)
 *  \brief Wavelength of H-10*/
double lambda_H_10(void)
{
	//in microns
	return lambda_H(2,10);
}
/*! \fn double lambda_H_15(void)
 *  \brief Wavelength of H-15*/
double lambda_H_15(void)
{
	//in microns
	return lambda_H(2,15);
}
/*! \fn double lambda_H_20(void)
 *  \brief Wavelength of H-20*/
double lambda_H_20(void)
{
	//in microns
	return lambda_H(2,20);
}
/*! \fn double F_H_beta(void)
 *  \brief H-beta line emission strength relative to H-beta */
double F_H_beta(void)
{
	//relative to Hbeta
	return 1.0;
}
/*! \fn double F_H_alpha(void)
 *  \brief H-alpha line emission strength relative to H-beta */
double F_H_alpha(void)
{
	//relative to Hbeta
	return 2.863;
}
/*! \fn double F_H_gamma(void)
 *  \brief H-gamma line emission strength relative to H-beta */
double F_H_gamma(void)
{
	//relative to Hbeta
	return 0.468;
}
/*! \fn double F_H_delta(void)
 *  \brief H-delta line emission strength relative to H-beta */
double F_H_delta(void)
{
	//relative to Hbeta
	return 0.259;
}
/*! \fn double F_H_10(void)
 *  \brief H-10 line emission strength relative to H-beta */
double F_H_10(void)
{
	//relative to Hbeta
	return 0.0530;
}
/*! \fn double F_H_15(void)
 *  \brief H-15 line emission strength relative to H-beta */
double F_H_15(void)
{
	//relative to Hbeta
	return 0.01561;
}
/*! \fn double F_H_20(void)
 *  \brief H-20 line emission strength relative to H-beta */
double F_H_20(void)
{
	//relative to Hbeta
	return 0.00662;
}


//Paschen series

/*! \fn double lambda_H_Palpha(void)
 *  \brief Wavelength of Paschen-alpha line in microns. */
double lambda_H_Palpha(void)
{
	//in microns
	return lambda_H(3,4);
}
/*! \fn double lambda_H_Pbeta(void)
 *  \brief Wavelength of Paschen-beta line in microns. */
double lambda_H_Pbeta(void)
{
	//in microns
	return lambda_H(3,5);
}
/*! \fn double lambda_H_Pgamma(void)
 *  \brief Wavelength of Paschen-gamma line in microns. */
double lambda_H_Pgamma(void)
{
	//in microns
	return lambda_H(3,6);
}
/*! \fn double lambda_H_Pdelta(void)
 *  \brief Wavelength of Paschen-delta line in microns. */
double lambda_H_Pdelta(void)
{
	//in microns
	return lambda_H(3,7);
}
/*! \fn double lambda_H_P10(void)
 *  \brief Wavelength of Paschen-10 line in microns. */
double lambda_H_P10(void)
{
	//in microns
	return lambda_H(3,10);
}
/*! \fn double lambda_H_P15(void)
 *  \brief Wavelength of Paschen-15 line in microns. */
double lambda_H_P15(void)
{
	//in microns
	return lambda_H(3,15);
}
/*! \fn double lambda_H_P20(void)
 *  \brief Wavelength of Paschen-20 line in microns. */
double lambda_H_P20(void)
{
	//in microns
	return lambda_H(3,20);
}


/*! \fn double F_H_Palpha(void)
 *  \brief Emission strength of Paschen-alpha relative to H-beta */
double F_H_Palpha(void)
{
	//relative to Hbeta
	return 0.339;
}
/*! \fn double F_H_Pbeta(void)
 *  \brief Emission strength of Paschen-beta relative to H-beta */
double F_H_Pbeta(void)
{
	//relative to Hbeta
	return 0.163;
}
/*! \fn double F_H_Pgamma(void)
 *  \brief Emission strength of Paschen-gamma relative to H-beta */
double F_H_Pgamma(void)
{
	//relative to Hbeta
	return 0.0904;
}
/*! \fn double F_H_Pdelta(void)
 *  \brief Emission strength of Paschen-delta relative to H-beta */
double F_H_Pdelta(void)
{
	//relative to Hbeta
	return 0.0555;
}
/*! \fn double F_H_P10(void)
 *  \brief Emission strength of Paschen-10 relative to H-beta */
double F_H_P10(void)
{
	//relative to Hbeta
	return 0.0184;
}
/*! \fn double F_H_P15(void)
 *  \brief Emission strength of Paschen-15 relative to H-beta */
double F_H_P15(void)
{
	//relative to Hbeta
	return 0.00541;
}
/*! \fn double F_H_P20(void)
 *  \brief Emission strength of Paschen-20 relative to H-beta */
double F_H_P20(void)
{
	//relative to Hbeta
	return 0.00229;
}


/*! \fn double F_H_Balpha(void)
 *  \brief Emission strength of Brackett-alpha relative to H-beta */
double F_H_Balpha(void)
{
	//relative to Hbeta
	return 0.0802;
}

/*! \fn double F_H_Bbeta(void)
 *  \brief Emission strength of Brackett-beta relative to H-beta */
double F_H_Bbeta(void)
{
	//relative to Hbeta
	return 0.0455;
}

/*! \fn double F_H_Bgamma(void)
 *  \brief Emission strength of Brackett-gamma relative to H-beta */
double F_H_Bgamma(void)
{
	//relative to Hbeta
	return 0.0278;
}

/*! \fn double F_H_Bdelta(void)
 *  \brief Emission strength of Brackett-delta relative to H-beta */
double F_H_Bdelta(void)
{
	//relative to Hbeta
	return 0.0183;
}

/*! \fn double F_H_B10(void)
 *  \brief Emission strength of Brackett-10 relative to H-beta */
double F_H_B10(void)
{
	//relative to Hbeta
	return 0.00914;
}

/*! \fn double F_H_B15(void)
 *  \brief Emission strength of Brackett-15 relative to H-beta */
double F_H_B15(void)
{
	//relative to Hbeta
	return 0.00266;
}

/*! \fn double F_H_B20(void)
 *  \brief Emission strength of Brackett-20 relative to H-beta */
double F_H_B20(void)
{
	//relative to Hbeta
	return 0.00112;
}

/*! \fn double lambda_H_Balpha(void)
 *  \brief Wavelength of H Brackett-alpha in microns. */
double lambda_H_Balpha(void)
{
	//in microns
	return lambda_H(4,5);
}

/*! \fn double lambda_H_Bbeta(void)
 *  \brief Wavelength of H Brackett-beta in microns. */
double lambda_H_Bbeta(void)
{
	//in microns
	return lambda_H(4,6);
}

/*! \fn double lambda_H_Bgamma(void)
 *  \brief Wavelength of H Brackett-gamma in microns. */
double lambda_H_Bgamma(void)
{
	//in microns
	return lambda_H(4,7);
}

/*! \fn double lambda_H_Bdelta(void)
 *  \brief Wavelength of H Brackett-delta in microns. */
double lambda_H_Bdelta(void)
{
	//in microns
	return lambda_H(4,8);
}

/*! \fn double lambda_H_B10(void)
 *  \brief Wavelength of H Brackett-10 in microns. */
double lambda_H_B10(void)
{
	//in microns
	return lambda_H(4,10);
}

/*! \fn double lambda_H_B15(void)
 *  \brief Wavelength of H Brackett-15 in microns. */
double lambda_H_B15(void)
{
	//in microns
	return lambda_H(4,15);
}

/*! \fn double lambda_H_B20(void)
 *  \brief Wavelength of H Brackett-20 in microns. */
double lambda_H_B20(void)
{
	//in microns
	return lambda_H(4,20);
}




//HeII recombination lines at ne=100, T=10^4, Case B, Table 4.5 of Osterbrock & Ferland

/*! \fn double F_HeII_43(void)
 *  \brief Line strength relative to HeII 4686.*/
double F_HeII_43(void)
{	
	//same as 4686
	return 1.0;
}

/*! \fn double F_HeII_norm(double n_Hepp, double n_p)
 *  \brief Normalization of the HeII 4686 (n=4 -> n=3) line strength
 *  relative to H-beta */
double F_HeII_norm(double n_Hepp, double n_p)
{
	//n=4 -> n=3
	//Table 4.4 + Table 4.5 of osterbrock and ferland
	//relative to H BETA!!!
	//4*pi*j_lambda/(ne * n_He++) in 10^-24 ergs cm^3 s^-1
	//1.52;

	//Hbeta : 4*pi*j_HBeta/ (n_e * n_p) = 1.23e-25 ergs cm^3 s^-1

	//j_HeII_4686/j_Hbeta = (1.52e-24 / 1.23e-25) * n_Hepp / n_p

	return (1.52e-24/1.23e-25)*(n_Hepp/n_p);
}
/*! \fn double F_HeII_4686(void)
 *  \brief Line strength relative to HeII 4686.*/
double F_HeII_4686(void)
{
	//n=4 -> n=3
	return 1.0;
}


/*! \fn double lambda_HeII_32(void)
 *  \brief Wavelength of HeII 3->2 line, in microns. */
double lambda_HeII_32(void)
{
	//in microns
	return lambda_HeII(2,3);
}

/*! \fn double lambda_HeII_42(void)
 *  \brief Wavelength of HeII 4->2 line, in microns. */
double lambda_HeII_42(void)
{
	//in microns
	return lambda_HeII(2,4);
}

/*! \fn double lambda_HeII_52(void)
 *  \brief Wavelength of HeII 5->2 line, in microns. */
double lambda_HeII_52(void)
{
	//in microns
	return lambda_HeII(2,5);
}

/*! \fn double lambda_HeII_62(void)
 *  \brief Wavelength of HeII 6->2 line, in microns. */
double lambda_HeII_62(void)
{
	//in microns
	return lambda_HeII(2,6);
}

/*! \fn double lambda_HeII_102(void)
 *  \brief Wavelength of HeII 10->2 line, in microns. */
double lambda_HeII_102(void)
{
	//in microns
	return lambda_HeII(2,10);
}

/*! \fn double lambda_HeII_152(void)
 *  \brief Wavelength of HeII 15->2 line, in microns. */
double lambda_HeII_152(void)
{
	//in microns
	return lambda_HeII(2,15);
}

/*! \fn double lambda_HeII_202(void)
 *  \brief Wavelength of HeII 20->2 line, in microns. */
double lambda_HeII_202(void)
{
	//in microns
	return lambda_HeII(2,20);
}

/*! \fn double F_HeII_32(void)
 *  \brief Line strength of HeII 3->2, relative to HeII 4->3.*/
double F_HeII_32(void)
{
	//relative to 4686
	return 6.474;
}
/*! \fn double F_HeII_42(void)
 *  \brief Line strength of HeII 4->2, relative to HeII 4->3.*/
double F_HeII_42(void)
{
	//relative to 4686 
	return 1.958;
}
/*! \fn double F_HeII_52(void)
 *  \brief Line strength of HeII 5->2, relative to HeII 4->3.*/
double F_HeII_52(void)
{
	//relative to 4686 
	return 0.873;
}
/*! \fn double F_HeII_62(void)
 *  \brief Line strength of HeII 6->2, relative to HeII 4->3.*/
double F_HeII_62(void)
{
	//relative to 4686 
	return 0.474;
}
/*! \fn double F_HeII_102(void)
 *  \brief Line strength of HeII 10->2, relative to HeII 4->3.*/
double F_HeII_102(void)
{
	//relative to 4686 
	return 0.0970;
}
/*! \fn double F_HeII_152(void)
 *  \brief Line strength of HeII 15->2, relative to HeII 4->3.*/
double F_HeII_152(void)
{
	//relative to 4686 
	return 0.02925;
}
/*! \fn double F_HeII_202(void)
 *  \brief Line strength of HeII 20->2, relative to HeII 4->3.*/
double F_HeII_202(void)
{
	//relative to 4686 
	return 0.01260;
}
/*! \fn double lambda_HeII_43(void)
 *  \brief Wavelength of HeII 4->3 line, in microns. */
double lambda_HeII_43(void)
{
	//in microns
	//4686
	return lambda_HeII(3,4);
}
/*! \fn double lambda_HeII_53(void)
 *  \brief Wavelength of HeII 5->3 line, in microns. */
double lambda_HeII_53(void)
{
	//in microns
	return lambda_HeII(3,5);
}
/*! \fn double lambda_HeII_63(void)
 *  \brief Wavelength of HeII 6->3 line, in microns. */
double lambda_HeII_63(void)
{
	//in microns
	return lambda_HeII(3,6);
}
/*! \fn double lambda_HeII_73(void)
 *  \brief Wavelength of HeII 7->3 line, in microns. */
double lambda_HeII_73(void)
{
	//in microns
	return lambda_HeII(3,7);
}
/*! \fn double lambda_HeII_103(void)
 *  \brief Wavelength of HeII 10->3 line, in microns. */
double lambda_HeII_103(void)
{
	//in microns
	return lambda_HeII(3,10);
}
/*! \fn double lambda_HeII_153(void)
 *  \brief Wavelength of HeII 15->3 line, in microns. */
double lambda_HeII_153(void)
{
	//in microns
	return lambda_HeII(3,15);
}
/*! \fn double lambda_HeII_203(void)
 *  \brief Wavelength of HeII 20->3 line, in microns. */
double lambda_HeII_203(void)
{
	//in microns
	return lambda_HeII(3,20);
}
/*! \fn double F_HeII_53(void)
 *  \brief Line strength of HeII 5->3, relative to HeII 4->3.*/
double F_HeII_53(void)
{
	//relative to 4686 (4-3)
	return 0.405;
}
/*! \fn double F_HeII_63(void)
 *  \brief Line strength of HeII 6->3, relative to HeII 4->3.*/
double F_HeII_63(void)
{
	//relative to 4686 
	return 0.209;
}
/*! \fn double F_HeII_73(void)
 *  \brief Line strength of HeII 7->3, relative to HeII 4->3.*/
double F_HeII_73(void)
{
	//relative to 4686 
	return 0.123;
}
/*! \fn double F_HeII_103(void)
 *  \brief Line strength of HeII 10->3, relative to HeII 4->3.*/
double F_HeII_103(void)
{
	//relative to 4686 
	return 0.0397;
}
/*! \fn double F_HeII_103(void)
 *  \brief Line strength of HeII 10->3, relative to HeII 4->3.*/
double F_HeII_153(void)
{
	//relative to 4686 
	return 0.01166;
}
/*! \fn double F_HeII_203(void)
 *  \brief Line strength of HeII 20->3, relative to HeII 4->3.*/
double F_HeII_203(void)
{
	//relative to 4686 
	return 0.00498;
}
/*! \fn double F_HeII_54(void)
 *  \brief Line strength of HeII 5->4, relative to HeII 4->3.*/
double F_HeII_54(void)
{
	//relative to 4686 
	return 0.275;
}
/*! \fn double lambda_HeII_54(void)
 *  \brief Wavelength of HeII 5->4 line, in microns. */
double lambda_HeII_54(void)
{
	//in microns
	return lambda_HeII(4,5);
}
/*! \fn double lambda_HeII_64(void)
 *  \brief Wavelength of HeII 6->4 line, in microns. */
double lambda_HeII_64(void)
{
	//in microns
	return lambda_HeII(4,6);
}
/*! \fn double lambda_HeII_74(void)
 *  \brief Wavelength of HeII 7->4 line, in microns. */
double lambda_HeII_74(void)
{
	//in microns
	return lambda_HeII(4,7);
}
/*! \fn double lambda_HeII_84(void)
 *  \brief Wavelength of HeII 8->4 line, in microns. */
double lambda_HeII_84(void)
{
	//in microns
	return lambda_HeII(4,8);
}
/*! \fn double lambda_HeII_104(void)
 *  \brief Wavelength of HeII 10->4 line, in microns. */
double lambda_HeII_104(void)
{
	//in microns
	return lambda_HeII(4,10);
}
/*! \fn double lambda_HeII_154(void)
 *  \brief Wavelength of HeII 15->4 line, in microns. */
double lambda_HeII_154(void)
{
	//in microns
	return lambda_HeII(4,15);
}
/*! \fn double lambda_HeII_204(void)
 *  \brief Wavelength of HeII 20->4 line, in microns. */
double lambda_HeII_204(void)
{
	//in microns
	return lambda_HeII(4,20);
}

/*! \fn double F_HeII_64(void)
 *  \brief Line strength of HeII 6->4, relative to HeII 4->3.*/
double F_HeII_64(void)
{
	//relative to 4686 
	return 0.134;
}
/*! \fn double F_HeII_74(void)
 *  \brief Line strength of HeII 7->4, relative to HeII 4->3.*/
double F_HeII_74(void)
{
	//relative to 4686 
	return 0.0759;
}
/*! \fn double F_HeII_84(void)
 *  \brief Line strength of HeII 8->4, relative to HeII 4->3.*/
double F_HeII_84(void)
{
	//relative to 4686 
	return 0.0478;
}
/*! \fn double F_HeII_104(void)
 *  \brief Line strength of HeII 10->4, relative to HeII 4->3.*/
double F_HeII_104(void)
{
	//relative to 4686 
	return 0.0230;
}
/*! \fn double F_HeII_154(void)
 *  \brief Line strength of HeII 15->4, relative to HeII 4->3.*/
double F_HeII_154(void)
{
	//relative to 4686 
	return 0.00655;
}
/*! \fn double F_HeII_204(void)
 *  \brief Line strength of HeII 20->4, relative to HeII 4->3.*/
double F_HeII_204(void)
{
	//relative to 4686 
	return 0.00276;
}


//HeI Case B lines, ne=100, T = 10^4, table 4.6 of Osterbrock & Ferland

/*! \fn double F_HeI_norm(double n_Hep, double n_p)
 *  \brief Normalization of HeII 4471 line strength relative to H-beta. */
double F_HeI_norm(double n_Hep, double n_p)
{
	//tables 4.6 and 4.4 of osterbrock and ferland
	//relative to H-BETA!!!!
	//4*pi*j_HeI_4471/(n_e * n_Hep) = 0.612e-25 erg cm^3 s^-1

	//Hbeta : 4*pi*j_HBeta/ (n_e * n_p) = 1.23e-25 ergs cm^3 s^-1

	//j_HeI_4471/j_Hbeta = (0.612e-25 / 1.23e-25) * n_Hep / n_p

	return (0.612e-25/1.23e-25) * (n_Hep/n_p);
}
/*! \fn double F_HeI_4471(void)
 *  \brief Normalization of HeI 4471 line strength relative to HeI 4471. */
double F_HeI_4471(void)
{
	//relative to 4471
	return 1.0;
}
/*! \fn double F_HeI_5876(void)
 *  \brief Normalization of HeI 5876 line strength relative to HeI 4471. */
double F_HeI_5876(void)
{
	//relative to 4471
	return 2.67;
}
/*! \fn double F_HeI_4026(void)
 *  \brief Normalization of HeI 4026 line strength relative to HeI 4471. */
double F_HeI_4026(void)
{
	//relative to 4471
	return 0.476;
}
/*! \fn double F_HeI_7065(void)
 *  \brief Normalization of HeI 7065 line strength relative to HeI 4471. */
double F_HeI_7065(void)
{
	//relative to 4471
	return 0.489;
}
/*! \fn double F_HeI_10830(void)
 *  \brief Normalization of HeI 10830 line strength relative to HeI 4471. */
double F_HeI_10830(void)
{
	//relative to 4471
	return 5.41;
}
/*! \fn double F_HeI_3889(void)
 *  \brief Normalization of HeI 3889 line strength relative to HeI 4471. */
double F_HeI_3889(void)
{
	//relative to 4471
	return 2.31;
}
/*! \fn double F_HeI_3187(void)
 *  \brief Normalization of HeI 3187 line strength relative to HeI 4471. */
double F_HeI_3187(void)
{
	//relative to 4471
	return 0.917;
}
/*! \fn double F_HeI_6678(void)
 *  \brief Normalization of HeI 6678 line strength relative to HeI 4471. */
double F_HeI_6678(void)
{
	//relative to 4471
	return 0.756;
}
/*! \fn double F_HeI_4922(void)
 *  \brief Normalization of HeI 4922 line strength relative to HeI 4471. */
double F_HeI_4922(void)
{
	//relative to 4471
	return 0.270;
}
/*! \fn double F_HeI_5016(void)
 *  \brief Normalization of HeI 5016 line strength relative to HeI 4471. */
double F_HeI_5016(void)
{
	//relative to 4471
	return 0.578;
}
/*! \fn double F_HeI_3965(void)
 *  \brief Normalization of HeI 3965 line strength relative to HeI 4471. */
double F_HeI_3965(void)
{
	//relative to 4471
	return 0.230;
}

/*! \fn double lambda_HeI_5876(void)
 *  \brief Wavelength of HeI 5876 in microns. */
double lambda_HeI_5876(void)
{
	//microns
	return 0.5876;
}
/*! \fn double lambda_HeI_4026(void)
 *  \brief Wavelength of HeI 4026 in microns. */
double lambda_HeI_4026(void)
{
	//microns
	return 0.4026;
}
/*! \fn double lambda_HeI_7065(void)
 *  \brief Wavelength of HeI 7065 in microns. */
double lambda_HeI_7065(void)
{
	//microns
	return 0.7065;
}
/*! \fn double lambda_HeI_10830(void)
 *  \brief Wavelength of HeI 10830 in microns. */
double lambda_HeI_10830(void)
{
	//microns
	return 1.083;
}
/*! \fn double lambda_HeI_3889(void)
 *  \brief Wavelength of HeI 3889 in microns. */
double lambda_HeI_3889(void)
{
	//microns
	return 0.3889;
}
/*! \fn double lambda_HeI_3187(void)
 *  \brief Wavelength of HeI 3187 in microns. */
double lambda_HeI_3187(void)
{
	//microns
	return 0.3187;
}
/*! \fn double lambda_HeI_6678(void)
 *  \brief Wavelength of HeI 6678 in microns. */
double lambda_HeI_6678(void)
{
	//microns
	return 0.6678;
}
/*! \fn double lambda_HeI_4922(void)
 *  \brief Wavelength of HeI 4922 in microns. */
double lambda_HeI_4922(void)
{
	//microns
	return 0.4922;
}
/*! \fn double lambda_HeI_5016(void)
 *  \brief Wavelength of HeI 5016 in microns. */
double lambda_HeI_5016(void)
{
	//microns
	return 0.5016;
}
/*! \fn double lambda_HeI_3965(void)
 *  \brief Wavelength of HeI 3965 in microns. */
double lambda_HeI_3965(void)
{
	//microns
	return 0.3965;
}



//metals

/*! \var double Z_line[5]
 *  \brief Array of metallicities at which the metal line strengths are calculated.*/
double Z_line[5] = {log10(4e-4),log10(4e-3),log10(8e-3),log10(0.02),log10(0.05)};

/*! \fn double F_CII_1335(double Z)
 *  \brief Strength of CII 1335 line relative to H-beta. */
double F_CII_1335(double Z)
{
	//relative to Hbeta
	if(log10(Z)<Z_line[2])
		return 0;
	return 0.110;
}
/*! \fn double F_OIII_1663(double Z)
 *  \brief Strength of OIII 1663 line relative to H-beta. */
double F_OIII_1663(double Z)
{
	//relative to Hbeta
	if(log10(Z)<Z_line[1])
		return 0;
	if(log10(Z)>Z_line[2])
		return 0.010;
	return (0.010 - 0.058)*( log10(Z)-Z_line[1] )/(Z_line[2] - Z_line[1]) + 0.058;
}
/*! \fn double F_CIII_1909(double Z)
 *  \brief Strength of CIII 1909 line relative to H-beta. */
double F_CIII_1909(double Z)
{
	//relative to Hbeta
	if(log10(Z)<Z_line[2])
		return 0;
	return 0.180;
}
/*! \fn double F_NII_2141(double Z)
 *  \brief Strength of NII 2141 line relative to H-beta. */
double F_NII_2141(double Z)
{
	//relative to Hbeta
	if(log10(Z)<Z_line[2])
		return 0;
	return 0.010;
}
/*! \fn double F_CII_2326(double Z)
 *  \brief Strength of CII 2336 line relative to H-beta. */
double F_CII_2326(double Z)
{
	//relative to Hbeta
	if(log10(Z)<Z_line[2])
		return 0;
	return 0.290;
}
/*! \fn double F_MgII_2798(double Z)
 *  \brief Strength of MgII 2798 line relative to H-beta. */
double F_MgII_2798(double Z)
{
	//relative to Hbeta
	if(log10(Z)<Z_line[1])
		return 0;
	if(log10(Z)>Z_line[2])
		return 0.070;
	return (0.070 - 0.310)*( log10(Z)-Z_line[1] )/(Z_line[2] - Z_line[1]) + 0.310;
}
/*! \fn double F_OII_3727(double Z)
 *  \brief Strength of OII 3727 line relative to H-beta. */
double F_OII_3727(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.489,1.791,3.010,3.010,3.010};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_NeIII_3869(double Z)
 *  \brief Strength of NeIII 3869 line relative to H-beta. */
double F_NeIII_3869(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.295,0.416,0.3,0.3,0.3};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_H_zeta_HeII_3889(double Z)
 *  \brief Strength of Hzeta+HeII 3889 line relative to H-beta. */
double F_H_zeta_HeII_3889(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.203,0.192,0.107,0.107,0.107};
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_H_epsilon_HeII_3970(double Z)
 *  \brief Strength of Hepsilon+HeII 3970 line relative to H-beta. */
double F_H_epsilon_HeII_3970(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.270,0.283,0.159,0.159,0.159};
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_HeI_4026(double Z)
 *  \brief Strength of HeI 4026 line relative to H-beta. */
double F_HeI_4026(double Z)
{
	//relative to Hbeta
	return 0.015;
}
/*! \fn double F_SII_4069(double Z)
 *  \brief Strength of SII 4069 line relative to H-beta. */
double F_SII_4069(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.005,0.017,0.029,0.029,0.029};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SII_4076(double Z)
 *  \brief Strength of SII 4076 line relative to H-beta. */
double F_SII_4076(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.002,0.007,0.011,0.011,0.011};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_OIII_4363(double Z)
 *  \brief Strength of OIII 4363 line relative to H-beta. */
double F_OIII_4363(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.109,0.066,0.010,0.010,0.010};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*double F_HeI_4471(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.036,0.036,0.050,0.050,0.050};
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}*/
/*! \fn double F_ArIV_HeI_4711(double Z)
 *  \brief Strength of ArIV and HeI 4711 line relative to H-beta. */
double F_ArIV_HeI_4711(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.010,0.014,0.0,0.0,0.0};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_OIII_4959(double Z)
 *  \brief Strength of OIII 4959 line relative to H-beta. */
double F_OIII_4959(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {1.097,1.617,1.399,1.399,1.399};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_OIII_5007(double Z)
 *  \brief Strength of OIII 5007 line relative to H-beta. */
double F_OIII_5007(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {3.159,4.752,4.081,4.081,4.081};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_NI_5199(double Z)
 *  \brief Strength of NI 5199 line relative to H-beta. */
double F_NI_5199(double Z)
{
	//relative to Hbeta

	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.003,0.010,0.030,0.030,0.030};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	double f = (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
	//printf("i %d l %e %e f %e\n",i,l[i],l[i+1],f);
	return f;
}
/*! \fn double F_NII_5755(double Z)
 *  \brief Strength of NII 5755 line relative to H-beta. */
double F_NII_5755(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.000,0.000,0.010,0.010,0.010};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_HeI_5876(double Z)
 *  \brief Strength of HeI 5876 line relative to H-beta. */
double F_HeI_5876(double Z)
{
	//relative to Hbeta

	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.096,0.108,0.140,0.140,0.140};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_OI_6300(double Z)
 *  \brief Strength of OI 6300 line relative to H-beta. */
double F_OI_6300(double Z)
{
	//relative to Hbeta

	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.008,0.041,0.130,0.130,0.130};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SIII_6312(double Z)
 *  \brief Strength of SIII 6312 line relative to H-beta. */
double F_SIII_6312(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.009,0.017,0.030,0.030,0.030};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_NII_6548(double Z)
 *  \brief Strength of NII 6548 line relative to H-beta. */
double F_NII_6548(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.005,0.059,0.136,0.136,0.136};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_NII_6583(double Z)
 *  \brief Strength of NII 6583 line relative to H-beta. */
double F_NII_6583(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.015,0.175,0.404,0.404,0.404};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SII_6716(double Z)
 *  \brief Strength of SII 6716 line relative to H-beta. */
double F_SII_6716(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.037,0.188,0.300,0.300,0.300};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SII_6730(double Z)
 *  \brief Strength of SII 6730 line relative to H-beta. */
double F_SII_6730(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.029,0.138,0.210,0.210,0.210};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_HeI_7065(double Z)
 *  \brief Strength of HeI 7065 line relative to H-beta. */
double F_HeI_7065(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.028,0.023,0.040,0.040,0.040};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_ArIII_7136(double Z)
 *  \brief Strength of ArIII 7136 line relative to H-beta. */
double F_ArIII_7136(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.027,0.071,0.035,0.035,0.035};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_OII_7320(double Z)
 *  \brief Strength of OII 7320 line relative to H-beta. */
double F_OII_7320(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.012,0.027,0.026,0.026,0.026};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_OII_7331(double Z)
 *  \brief Strength of OII 7331 line relative to H-beta. */
double F_OII_7331(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.007,0.014,0.014,0.014,0.014};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_ArIII_7751(double Z)
 *  \brief Strength of ArIII 7751 line relative to H-beta. */
double F_ArIII_7751(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.067,0.176,0.086,0.086,0.086};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SIII_9069(double Z)
 *  \brief Strength of SIII 9069 line relative to H-beta. */
double F_SIII_9069(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.0,0.510,0.945,0.945,0.945};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SIII_9531(double Z)
 *  \brief Strength of SIII 9531 line relative to H-beta. */
double F_SIII_9531(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.0,0.0,0.365,0.365,0.365};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SII_10287(double Z)
 *  \brief Strength of SII 10287 line relative to H-beta. */
double F_SII_10287(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.000,0.000,0.048,0.048,0.048};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SII_10320(double Z)
 *  \brief Strength of SII 10320 line relative to H-beta. */
double F_SII_10320(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.000,0.000,0.058,0.058,0.058};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}
/*! \fn double F_SII_10336(double Z)
 *  \brief Strength of SII 10336 line relative to H-beta. */
double F_SII_10336(double Z)
{
	//relative to Hbeta
	size_t i = gsl_interp_bsearch(Z_line,log10(Z),0,3);
	double l[5] = {0.000,0.000,0.054,0.054,0.054};
	if(log10(Z)<Z_line[0])
		return l[0] * Z/pow(10,Z_line[0]);
	return (l[i+1]-l[i])*(log10(Z)-Z_line[i])/(Z_line[i+1]-Z_line[i]) + l[i];
}





/*! \fn double lambda_CII_1335(double Z)
 *  \brief Wavelength of CII 1335 line in microns.*/
double lambda_CII_1335(void)
{
	//microns
	return 0.1335;
}
/*! \fn double lambda_OIII_1663(double Z)
 *  \brief Wavelength of OIII 1663 line in microns.*/
double lambda_OIII_1663(void)
{
	//microns
	return 0.1663;
}
/*! \fn double lambda_CIII_1909(double Z)
 *  \brief Wavelength of CIII 1909 line in microns.*/
double lambda_CIII_1909(void)
{
	//microns
	return 0.1909;
}
/*! \fn double lambda_NII_2141(double Z)
 *  \brief Wavelength of NII 2141 line in microns.*/
double lambda_NII_2141(void)
{
	//microns
	return 0.2141;
}

/*! \fn double lambda_CII_2326(double Z)
 *  \brief Wavelength of CII 2326 line in microns.*/
double lambda_CII_2326(void)
{
	//microns
	return 0.2326;
}

/*! \fn double lambda_MgII_2798(double Z)
 *  \brief Wavelength of MgII 2798 line in microns.*/
double lambda_MgII_2798(void)
{
	//microns
	return 0.2798;
}


/*! \fn double lambda_OII_3727(double Z)
 *  \brief Wavelength of OII 3727 line in microns.*/
double lambda_OII_3727(void)
{
	//microns
	return 0.3727;
}

/*! \fn double lambda_NeIII_3869(double Z)
 *  \brief Wavelength of NeIII 3869 line in microns.*/
double lambda_NeIII_3869(void)
{
	//microns
	return 0.3869;
}

/*! \fn double lambda_H_zeta_HeII_3889(double Z)
 *  \brief Wavelength of Hzeta and HeII 3889 line in microns.*/
double lambda_H_zeta_HeII_3889(void)
{
	//microns
	return 0.3889;
}
/*! \fn double lambda_H_epsilon_HeII_3970(double Z)
 *  \brief Wavelength of Hepsilon and HeII 3970 line in microns.*/
double lambda_H_epsilon_HeII_3970(void)
{
	//microns
	return 0.3970;
}
/*double lambda_HeI_4026(void)
{
	//microns
	return 0.4026;
}*/
/*! \fn double lambda_SII_4069(double Z)
 *  \brief Wavelength of SII 4069 line in microns.*/
double lambda_SII_4069(void)
{
	//microns
	return 0.4069;
}
/*! \fn double lambda_SII_4076(double Z)
 *  \brief Wavelength of SII 4076 line in microns.*/
double lambda_SII_4076(void)
{
	//microns
	return 0.4076;
}
/*! \fn double lambda_OIII_4363(double Z)
 *  \brief Wavelength of OIII 4363 line in microns.*/
double lambda_OIII_4363(void)
{
	//microns
	return 0.4363;
}
/*! \fn double lambda_HeI_4471(double Z)
 *  \brief Wavelength of HeI 4471 line in microns.*/
double lambda_HeI_4471(void)
{
	//microns
	return 0.4471;
}
/*! \fn double lambda_ArIV_HeI_4711(double Z)
 *  \brief Wavelength of ArIV and HeI 4711 line in microns.*/
double lambda_ArIV_HeI_4711(void)
{
	//microns
	return 0.4711;
}
/*! \fn double lambda_OIII_4959(double Z)
 *  \brief Wavelength of OIII 4959 line in microns.*/
double lambda_OIII_4959(void)
{
	//microns
	return 0.4959;
}
/*! \fn double lambda_OIII_5007(double Z)
 *  \brief Wavelength of OIII 5007 line in microns.*/
double lambda_OIII_5007(void)
{
	//microns
	return 0.5007;
}
/*! \fn double lambda_NI_5199(double Z)
 *  \brief Wavelength of NI 5199 line in microns.*/
double lambda_NI_5199(void)
{
	//microns
	return 0.5199;
}

/*! \fn double lambda_NI_5755(double Z)
 *  \brief Wavelength of NI 5755 line in microns.*/
double lambda_NII_5755(void)
{
	//microns
	return 0.5755;
}
/*double lambda_HeI_5876(void)
{
	//microns
	return 0.5876;
}*/
/*! \fn double lambda_OI_6300(double Z)
 *  \brief Wavelength of OI 6300 line in microns.*/
double lambda_OI_6300(void)
{
	//microns
	return 0.6300;
}
/*! \fn double lambda_SIII_6312(double Z)
 *  \brief Wavelength of SIII 6312 line in microns.*/
double lambda_SIII_6312(void)
{
	//microns
	return 0.6312;
}
/*! \fn double lambda_NII_6548(double Z)
 *  \brief Wavelength of NII 6548 line in microns.*/
double lambda_NII_6548(void)
{
	//microns
	return 0.654805;
}
/*! \fn double lambda_NII_6583(double Z)
 *  \brief Wavelength of NII 6583 line in microns.*/
double lambda_NII_6583(void)
{
	//microns
	return 0.658345;
}
/*! \fn double lambda_SII_6716(double Z)
 *  \brief Wavelength of SII 6716 line in microns.*/
double lambda_SII_6716(void)
{
	//microns
	return 0.6716;
}
/*! \fn double lambda_SII_6730(double Z)
 *  \brief Wavelength of SII 6730 line in microns.*/
double lambda_SII_6730(void)
{
	//microns
	return 0.6730;
}
/*double lambda_HeI_7065(void)
{
	//microns
	return 0.7065;
}*/
/*! \fn double lambda_ArIII_7136(double Z)
 *  \brief Wavelength of ArIII 7136 line in microns.*/
double lambda_ArIII_7136(void)
{
	//microns
	return 0.7136;
}
/*! \fn double lambda_SIII_9069(double Z)
 *  \brief Wavelength of SIII 9069 line in microns.*/
double lambda_SIII_9069(void)
{
	//microns
	return 0.9069;
}
/*! \fn double lambda_OII_7320(double Z)
 *  \brief Wavelength of OII 7320 line in microns.*/
double lambda_OII_7320(void)
{
	//microns
	return 0.731999;
}
/*! \fn double lambda_OII_7331(double Z)
 *  \brief Wavelength of OII 7331 line in microns.*/
double lambda_OII_7331(void)
{
	//microns
	return 0.733073;
}
/*! \fn double lambda_ArIII_7751(double Z)
 *  \brief Wavelength of ArIII 7751 line in microns.*/
double lambda_ArIII_7751(void)
{
	//microns
	return 0.775111;
}
/*! \fn double lambda_SIII_9531(double Z)
 *  \brief Wavelength of SIII 9531 line in microns.*/
double lambda_SIII_9531(void)
{
	//microns
	return 0.9531;
}
/*! \fn double lambda_SII_10287(double Z)
 *  \brief Wavelength of SII 10287 line in microns.*/
double lambda_SII_10287(void)
{
	//microns
	return 1.028673;
}
/*! \fn double lambda_SII_10320(double Z)
 *  \brief Wavelength of SII 10320 line in microns.*/
double lambda_SII_10320(void)
{
	//microns
	return 1.032049;
}
/*! \fn double lambda_SII_10336(double Z)
 *  \brief Wavelength of SII 10336 line in microns.*/
double lambda_SII_10336(void)
{
	//microns
	return 1.033641;
}

