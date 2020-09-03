#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "bpass_models.h"
#include "igm-absorption.h"
#include "pass_bands.h"
#include "nebular_line_emission.h"
#include <gsl/gsl_interp.h>

#define CHECK_PHOT

#define FAKE_FESC

//#define FIT_LINE_STRENGTHS

#define LINES
#define CONT

char **load_z_met_filenames(void);
void allocate_sed_buffers(void);

double klambda_smc(double lambda)
{
  //note that lambda is in microns
  double lambda_smc[31] = {0.09,0.116,0.119,0.123,0.127,0.131,0.136,0.140,0.145,0.151,0.157,0.163,0.170,0.178,0.186,0.195,0.205,0.216,0.229,0.242,0.258,0.276,0.296,0.370,0.440,0.550,0.660,0.810,1.250,1.650,2.198};
  double k_smc[31] = {9.5,6.992,6.436,6.297,6.074,5.795,5.575,5.272,5.000,4.776,4.472,4.243,4.013,3.866,3.637,3.489,3.293,3.161,2.947,2.661,2.428,2.220,2.000,1.672,1.374,1.000,0.801,0.567,0.131,0.169,0.016};
  double Rv_smc = 2.74; // gordon et al. 2003
  int nlsmc = 31;
  int    ki = gsl_interp_bsearch(lambda_smc,lambda,0,nlsmc);
  double ksmc =  k_smc[ki] + (k_smc[ki+1] - k_smc[ki])*(lambda - lambda_smc[ki])/(lambda_smc[ki+1] - lambda_smc[ki]);
  return ksmc * Rv_smc;
}




int main(int argc, char **argv)
{
  double z     = 2.1724; //redshift
  double A     = 16.1;  //SFR in Msun/yr
  double agei  = 6.5;   //age in log10 years
  double f_esc = 0.66;  //escape fraction
  double dm;
  double Nlyc;
  double ebv = 0.0;

  double    AEBV, kp; //dust attenuation

  double lyman_limit = 912.0; //Ang

#ifdef FAKE_FESC
  lyman_limit = 970.;
#endif

  flag_binary = 1; //binary population flag

  if(argc>=2)
    z = atof(argv[1]);
	if(argc>=3)
		A = atof(argv[2]);
	if(argc>=4)
		agei = atof(argv[3]);
	if(argc>=5)
		f_esc = atof(argv[4]);
#ifndef FIT_LINE_STRENGTHS

  if(argc>=6)
    ebv = atof(argv[5]);
  if(argc>=7)
    flag_binary = atoi(argv[6]);


#else //FIT_LINE_STRENGTHS
  double B_line;
  if(argc>=6)
    B_line = atof(argv[5]);
  if(argc>=7)
    flag_binary = atoi(argv[6]);
#endif //FIT_LINE_STRENGTHS


  dm = distance_modulus(z);
  printf("z     = %f\n",z);
  printf("A     = %f [Msun/yr]\n",A);
  printf("age   = %f [log10 yr]\n",agei);
  printf("f_esc = %f\n",f_esc);
  printf("ebv   = %f\n",ebv);

#ifndef FIT_LINE_STRENGTHS
  if(argc>=6)
    printf("flag_binary = %d\n",flag_binary);
#else //FIT_LINE_STRENGTHS
  if(argc>=6)
    printf("B_line = %e\n",B_line);
  if(argc>=7)
    printf("flag_binary = %d\n",flag_binary);
#endif //FIT_LINE_STRENGTHS
  printf("dm    = %e\n",dm);

  //indices
  int i,j,k,l;

  //flag for switching binary fitting
  //flag_binary = 0 -- only fitting single
  //flag_binary = 1 -- only fitting binary
  //flag_binary = 2 -- fitting range of binarity

  //do we need to load the IGM absorption?
  int flag_igm_abs = 1;

  //filename properties
  char fname_single[200];
  char fname_binary[200];
  char fname_ion[200];
  char fname_base[200];
  char fname_sfr[200];
  char fname_suff[200];
  char fname[200];
  char fname_line[200];
  char fname_cont[200];

  initialize_bpass_models();


  char **fname_z_met;

  //load the metallicity file name components
  fname_z_met = load_z_met_filenames();

  //the directory containing the single
  //star SED models
  sprintf(fname_single,"single");

  //the directory containing the binary
  //star SED models
  if(flag_binary==0)
  {
    sprintf(fname_binary,"single"); 

  }else{
    sprintf(fname_binary,"binary"); 
  }

  //the base filename
  sprintf(fname_base,"fnu_cgs");

  //the base sfr extension
  sprintf(fname_sfr,"1Msun_per_year");

  //the base file suffix
  sprintf(fname_suff,"txt");

  //allocate the sed buffers
  allocate_sed_buffers();
  //allocate_sed_samples();

  //allocate interpolated buffers
  double *f_nu_sed_lo   = (double *) malloc(n_lambda*sizeof(double));
  double *f_nu_sed_hi   = (double *) malloc(n_lambda*sizeof(double));
  double *f_nu_line_lo   = (double *) malloc(n_lambda*sizeof(double));
  double *f_nu_line_hi   = (double *) malloc(n_lambda*sizeof(double));
  double *f_nu_cont_lo   = (double *) malloc(n_lambda*sizeof(double));
  double *f_nu_cont_hi   = (double *) malloc(n_lambda*sizeof(double));


  //Let's do Lyman continuum first
  sprintf(fname_ion,"ionizing");


  for(j=0;j<n_z_met;j++)
  {
    sprintf(fname,"%s/%s.%s.imf135_100.%s.age_%s.%s",fname_binary,fname_ion,fname_binary,fname_sfr,fname_z_met[j],fname_suff);
    printf("fname %s\n",fname);
    FILE *fp = fopen(fname,"r");
    double ab, nlb, cb, db, eb;
    for(k=0;k<n_age;k++)
    {
      fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",&ab,&nlb,&cb,&db,&eb);
      N_Lyc[j][k] = nlb;
      //printf("age %e nlb %e\n",age[k],N_Lyc[j][k]);
    }
  }
  int     k_age = gsl_interp_bsearch(age,agei,0,n_age-2);
  double dx_age = (agei - age[k_age])/(age[k_age+1] - age[k_age]);

#ifdef CHECK_PHOT
  printf("PHOT k_age %d dx_age %e\n",k_age,dx_age);
#endif //CHECK_PHOT
  //exit(0);
  Nlyc = log10(A)+((1.0-dx_age)*N_Lyc[0][k_age] + dx_age*N_Lyc[0][k_age+1]);
  printf("N_Lyc = %e\n",Nlyc);


  //sprintf(fname_cont,"nebular_continuum_emission.HI_only.txt");
  sprintf(fname_cont,"nebular_continuum_emission.HI_and_HeI.txt");
  sprintf(fname_line,"line_emission.%s.txt",fname_z_met[0]);

  //do low ages first
  //load_line_model(fname_line, z, log10(A)+N_Lyc[0][k_age]);
  //load_cont_model(fname_cont, z, log10(A)+N_Lyc[0][k_age]);
  load_line_model(fname_line, z, N_Lyc[0][k_age]);
  load_cont_model(fname_cont, z, N_Lyc[0][k_age]);
  sprintf(fname,"%s/%s.%s.%s.age_%3.1f.%s.%s",fname_binary,fname_base,fname_binary,fname_sfr,age[k_age],fname_z_met[0],fname_suff);
  load_sed_model(fname,z);

  //get absorption
  double *absorption = igm_absorption(lambda_line, n_lambda_line, z);

  print_igm_absorption(lambda_line, n_lambda_line, z);

  //interpolate
  for(i=0;i<n_lambda_line;i++)
  {
    //nebular continuum
    //l = gsl_interp_bsearch(lambda_cont,lambda[i],0,n_lambda_cont-2);
    //f_nu_cont_lo[i] = f_nu_cont[l] + (lambda[i]-lambda_cont[l])/(lambda_cont[l+1] - lambda_cont[l])*f_nu_cont[l+1];
    f_nu_cont_lo[i] = f_nu_cont[i];

    //Lines
    //l = gsl_interp_bsearch(lambda_line,lambda[i],0,n_lambda_line-2);
    //f_nu_line_lo[i] = f_nu_line[l] + (lambda[i]-lambda_line[l])/(lambda_line[l+1] - lambda_line[l])*f_nu_line[l+1];
    f_nu_line_lo[i] = f_nu_line[i];

    //SED
    //f_nu_sed_lo[i] = f_nu[i];

    l = gsl_interp_bsearch(lambda,lambda_line[i],0,n_lambda);
    //printf("lambda %e %e %e\n",lambda[l],lambda_line[i],lambda[l+1]);
    f_nu_sed_lo[i] = f_nu[l] + (lambda_line[i]-lambda[l])/(lambda[l+1] - lambda[l])*(f_nu[l+1]-f_nu[l]);
  }


  //do low ages first
  //do low ages first
  //load_line_model(fname_line, z, log10(A)+N_Lyc[0][k_age+1]);
  //load_cont_model(fname_cont, z, log10(A)+N_Lyc[0][k_age+1]);
  load_line_model(fname_line, z, N_Lyc[0][k_age+1]);
  load_cont_model(fname_cont, z, N_Lyc[0][k_age+1]);
  sprintf(fname,"%s/%s.%s.%s.age_%3.1f.%s.%s",fname_binary,fname_base,fname_binary,fname_sfr,age[k_age+1],fname_z_met[0],fname_suff);
  load_sed_model(fname,z);

  //interpolate
  for(i=0;i<n_lambda_line;i++)
  {
    //nebular continuum
    //l = gsl_interp_bsearch(lambda_cont,lambda[i],0,n_lambda_cont-2);
    //f_nu_cont_hi[i] = f_nu_cont[l] + (lambda[i]-lambda_cont[l])/(lambda_cont[l+1] - lambda_cont[l])*f_nu_cont[l+1];
    f_nu_cont_hi[i] = f_nu_cont[i];

    //Lines
    //l = gsl_interp_bsearch(lambda_line,lambda[i],0,n_lambda_line-2);
    //f_nu_line_hi[i] = f_nu_line[l] + (lambda[i]-lambda_line[l])/(lambda_line[l+1] - lambda_line[l])*f_nu_line[l+1];
    f_nu_line_hi[i] = f_nu_line[i];

    //SED
    //f_nu_sed_hi[i] = f_nu[i];
    l = gsl_interp_bsearch(lambda,lambda_line[i],0,n_lambda);
    f_nu_sed_hi[i] = f_nu[l] + (lambda_line[i]-lambda[l])/(lambda[l+1] - lambda[l])*(f_nu[l+1]-f_nu[l]);

  }

  double f_sed;
  double f_line;
  double f_cont;
  double f_tot;
  double f_tot_nl;
  double f_tot_bare;

  double *f_nu_tot = (double *) malloc(n_lambda_line*sizeof(double));
  double *f_nu_tot_nl = (double *) malloc(n_lambda_line*sizeof(double));


  FILE *fpsed, *fpsed_bare;
  if(flag_binary==1)
  {
    fpsed = fopen("sed.txt","w");
  }else{
    fpsed = fopen("sed.single.txt","w");
  }
  if(flag_binary==1)
  {
    fpsed_bare = fopen("sed.bare.txt","w");
  }else{
    fpsed_bare = fopen("sed.bare.single.txt","w");
  }
  for(i=0;i<n_lambda_line;i++)
  {
    f_sed   = A*((1.0-dx_age)*f_nu_sed_lo[i]  + dx_age*f_nu_sed_hi[i]);
    //f_sed   = f_nu_sed_hi[i];
    //f_sed   = f_nu_sed_lo[i];
    //f_sed = 0.0;

    f_line  = A*((1.0-dx_age)*f_nu_line_lo[i] + dx_age*f_nu_line_hi[i]);
#ifdef FIT_LINE_STRENGTHS
    f_line *= B_line;
#endif //FIT_LINE_STRENGTHS
    f_cont  = A*((1.0-dx_age)*f_nu_cont_lo[i] + dx_age*f_nu_cont_hi[i]);

    f_tot = f_sed;
    f_tot_nl = f_sed;


#ifdef LINES
    f_tot += (1.0-f_esc)*f_line;
#endif 
#ifdef CONT
    f_tot += (1.0-f_esc)*f_cont;
    f_tot_nl += (1.0-f_esc)*f_cont;

#endif 


    //apply dust 
    kp = klambda_smc(lambda_line[i] * 1.0e-4);
    if(lambda_line[i]  * 1.0e-4 < lyman_limit * 1.0e-4)
      kp = 0.0;
    AEBV = -0.4*ebv*kp;


    f_tot    *= pow(10.0, AEBV);
    f_tot_nl *= pow(10.0, AEBV);


    //apply absorption
    f_tot *= absorption[i];
    f_tot_bare = f_tot;
    f_tot_nl *= absorption[i];


    if(lambda_line[i]<=lyman_limit)
    {
      f_tot *= f_esc;
      f_tot_nl *= f_esc;
    }


    f_nu_tot[i] = f_tot;
    f_nu_tot_nl[i] = f_tot_nl;

    fprintf(fpsed,"%e\t%e\t%e\t%e\t%e\t%e\n",lambda_line[i],f_sed,(1.0-f_esc)*f_line,(1.0-f_esc)*f_cont,absorption[i],f_tot);
    fprintf(fpsed_bare,"%e\t%e\t%e\t%e\t%e\t%e\n",lambda_line[i],f_sed,(1.0-f_esc)*f_line,(1.0-f_esc)*f_cont,absorption[i],f_tot_bare);
  }
  fclose(fpsed);
  fclose(fpsed_bare);


  //all data now has 11 data points
  //output model photometry
  int n_data = 11;
  double fT, T;
  int i_lo, i_hi;
  double oopz = 1./(1.+z);
  double f_nu_lo, f_nu_hi;
  double dx_lo, dx_hi;
  double *f_model = (double *) malloc(n_data*sizeof(double));

  //double lambda_model[13] = {3.360000e+03,3.850000e+03,4.500000e+03,4.975000e+03,5.500000e+03,6.500000e+03,7.600000e+03,8.900000e+03,1.250000e+04,1.550000e+04,2.200000e+04,3.600000e+04,4.500000e+04};
  double lambda_model[11] = {2750.0, 3375.0, 4317.4, 5917.7, 7693.0, 8059.8, 9054.8, 10450.0, 12500.0, 14000.0, 15450.0};
  int m;
  double *x_pb_o;
  double p_lo, p_hi;
  double nu_lo, nu_hi;
  double nu_ave;
  double c = 2.99792458e18; //speed of light in Ang/sec

#ifdef ALTERNATE_LINE_EMISSION
  double y_pb_max = 0;      //maximum transmissivity of pass band
  double f_pb_define = 0.5; //fraction of y_pb_max that defines width of pass band
  double l_pb_lo;    //minimum wavelength of pass band
  double l_pb_hi;    //maximum wavelength of pass band
  double l_pb_max = 0;
  double dlam_lo;
  double dlam_hi;
  double dlam_lo_sum;
  double dlam_hi_sum;
  double l_pb_ave;
  double T_pb_ave;
  double y_T;

  rescale_line_strengths(z, Nlyc);

#endif // ALTERNATE_LINE_EMISSION

  fpsed = fopen("sed.model_photometry.txt","w");
  for(l=0;l<n_data;l++)
  {

#ifdef ALTERNATE_LINE_EMISSION
    y_pb_max = 0;
    dlam_lo = 0;
    dlam_hi = 0;
    dlam_lo_sum = 0;
    dlam_hi_sum = 0;
    l_pb_ave = 0.0;
    T_pb_ave = 0.0;
#endif // ALTERNATE_LINE_EMISSION


    //note that each data point corresponds to a pass band
    //so we load the pass band
    LoadPassBand(l);

    x_pb_o = (double *) malloc(n_pb*sizeof(double));
    for(m=0;m<n_pb;m++)
    {
      x_pb_o[m] = x_pb[m]*oopz;
#ifdef ALTERNATE_LINE_EMISSION
      //find the maximum passband emissivity
      T_pb_ave += y_pb[m];
      l_pb_ave += y_pb[m]*x_pb[m];
      if(y_pb[m]>y_pb_max)
      {
        l_pb_max = x_pb_o[m];
        y_pb_max = y_pb[m];
      }
#endif //ALTERNATE_LINE_EMISSION
    }

#ifdef ALTERNATE_LINE_EMISSION
    //define the bandwidth of the passband
    l_pb_lo = -1;
    l_pb_hi = -1;
    l_pb_ave /= T_pb_ave; //average wavelength in Angstrom
    for(m=0;m<n_pb;m++)
    {
      if(x_pb_o[m]<l_pb_max)
      {
        dlam_lo += y_pb[m]*(l_pb_max-x_pb_o[m]);
        dlam_lo_sum += y_pb[m];
      }
      if(x_pb_o[m]>=l_pb_max)
      {
        dlam_hi += y_pb[m]*(x_pb_o[m]-l_pb_max);
        dlam_hi_sum += y_pb[m];
      }
      //find short wavelength side of band
      if( (y_pb[m]>f_pb_define*y_pb_max)&&(l_pb_lo<0) )
      {
        l_pb_lo = x_pb_o[m];
      }
      //find long wavelength side of band
      if( (y_pb[m]<f_pb_define*y_pb_max)&&(l_pb_hi<0)&&(l_pb_lo>0)&&(x_pb_o[m]>l_pb_max) )
      {
        l_pb_hi = x_pb_o[m];
      }
      //if(l==8)
        //printf("BANDWIDTH l %e y %e y_pb_max %e\n",x_pb_o[m],y_pb[m],y_pb_max);
    }
    //compute the pass band in Hz
    //dnu_pb = c/l_pb_lo - c/l_pb_hi;
    dlam_lo /= dlam_lo_sum;
    dlam_hi /= dlam_hi_sum;


    //printf("BANDWIDTH for %d l_ave %f; width %f %f l_pb_lo %e l_pb_max %e l_pb_hi %e y_pb_max %e dn_pb %e check %e %e %e\n",l,l_pb_ave,(dnu_hi + dnu_lo)/oopz,(l_pb_hi - l_pb_lo)/oopz,l_pb_lo/oopz,l_pb_max/oopz,l_pb_hi/oopz,y_pb_max,dnu_pb,(l_pb_max-dnu_lo)/oopz,(dnu_hi+l_pb_max)/oopz,c*oopz/(l_pb_max-dnu_lo) + c*oopz/(l_pb_max+dnu_hi));
    //printf("BANDWIDTH for %d l_ave %f; width %f %f l_pb_lo %e l_pb_max %e l_pb_hi %e y_pb_max %e\n",l,l_pb_ave,(dlam_hi + dlam_lo)/oopz,(l_pb_hi - l_pb_lo)/oopz,l_pb_lo/oopz,l_pb_max/oopz,l_pb_hi/oopz,y_pb_max);

#endif // ALTERNATE_LINE_EMISSION

    //loop over band curve points and integrate
    fT = 0;
    T  = 0;
    /*
    for(int m=0;m<n_pb-1;m++)
    {
      
      //index where x[i_lo] <= x < x[i_lo+1]
      //i_lo = gsl_interp_bsearch(lambda_line,x_pb[m]*oopz,0,n_lambda_cont-2);
      //i_hi = gsl_interp_bsearch(lambda_line,x_pb[m+1]*oopz,0,n_lambda_cont-2);
      //f_nu_lo  = (1.0-dx_lo)*f_nu_tot[i_lo] + dx_lo*f_nu_tot[i_lo+1];
      //f_nu_hi  = (1.0-dx_hi)*f_nu_tot[i_hi] + dx_hi*f_nu_tot[i_hi+1];
      //fractional distance between model points in wavelength
      //dx_lo = (x_pb[m]*oopz-lambda_line[i_lo])/(lambda_line[i_lo+1]-lambda_line[i_lo]);
      //dx_hi = (x_pb[m+1]*oopz-lambda_line[i_hi])/(lambda_line[i_hi+1]-lambda_line[i_hi]);
      //fT += 0.5*(f_nu_lo*y_pb[m] + f_nu_hi*y_pb[m+1])*(x_pb[m+1]-x_pb[m]) / (0.5*(x_pb[m+1]+x_pb[m]));
      //T  += 0.5*(y_pb[m] + y_pb[m+1])*(x_pb[m+1]-x_pb[m])  / (0.5*(x_pb[m+1]+x_pb[m]));
      

      //index where x[i_lo] <= x < x[i_lo+1]
      i_lo = gsl_interp_bsearch(lambda_line,x_pb[m]*oopz,0,n_lambda_cont-2);
      i_hi = gsl_interp_bsearch(lambda_line,x_pb[m+1]*oopz,0,n_lambda_cont-2);
      f_nu_lo  = (1.0-dx_lo)*f_nu_tot[i_lo] + dx_lo*f_nu_tot[i_lo+1];
      f_nu_hi  = (1.0-dx_hi)*f_nu_tot[i_hi] + dx_hi*f_nu_tot[i_hi+1];
      //fractional distance between model points in wavelength
      dx_lo = (x_pb[m]*oopz-lambda_line[i_lo])/(lambda_line[i_lo+1]-lambda_line[i_lo]);
      dx_hi = (x_pb[m+1]*oopz-lambda_line[i_hi])/(lambda_line[i_hi+1]-lambda_line[i_hi]);
      fT += 0.5*(f_nu_lo*y_pb[m] + f_nu_hi*y_pb[m+1])*(x_pb[m+1]-x_pb[m]) / (0.5*(x_pb[m+1]+x_pb[m]));
      T  += 0.5*(y_pb[m] + y_pb[m+1])*(x_pb[m+1]-x_pb[m])  / (0.5*(x_pb[m+1]+x_pb[m]));



    }*/

    for(m=0;m<n_lambda_cont-1;m++)
    {
      if( (lambda_cont[m]>=x_pb_o[0])&&(lambda_cont[m+1]<=x_pb_o[n_pb-1]))
      {
      //index where x[i_lo] <= x < x[i_lo+1]
      i_lo = gsl_interp_bsearch(x_pb_o,lambda_cont[m],0,n_pb);
      i_hi = gsl_interp_bsearch(x_pb_o,lambda_cont[m+1],0,n_pb);

      //fractional distance between model points in wavelength
      dx_lo = (lambda_cont[m]-x_pb_o[i_lo])/(x_pb_o[i_lo+1]-x_pb_o[i_lo]);
      dx_hi = (lambda_cont[m+1]-x_pb_o[i_hi])/(x_pb_o[i_hi+1]-x_pb_o[i_hi]);

      p_lo  = (1.0-dx_lo)*y_pb[i_lo] + dx_lo*y_pb[i_lo+1];
      p_hi  = (1.0-dx_hi)*y_pb[i_hi] + dx_hi*y_pb[i_hi+1];

      //printf("m %d p_lo %e p_hi %e dx_lo %e dx_hi %e oopz %e i_lo %d i_hi %d n_pb %d x_pb_o_hi %e %e\n",m,p_lo,p_hi,dx_lo,dx_hi,oopz,i_lo,i_hi,n_pb,x_pb_o[i_hi],x_pb_o[i_hi+1]);

      nu_hi = c/lambda_cont[m];
      nu_lo = c/lambda_cont[m+1];
      nu_ave = 0.5*(nu_hi + nu_lo);

#ifndef ALTERNATE_LINE_EMISSION
      fT += 0.5*(p_lo*f_nu_tot[m] + p_hi*f_nu_tot[m+1])*(nu_hi-nu_lo)/nu_ave;
#else  //ALTERNATE_LINE_EMISSION
      fT += 0.5*(p_lo*f_nu_tot_nl[m] + p_hi*f_nu_tot_nl[m+1])*(nu_hi-nu_lo)/nu_ave;
#endif //ALTERNATE_LINE_EMISSION

      T  += 0.5*(p_hi + p_lo)*(nu_hi-nu_lo)/nu_ave;    


      }
    }
    fT/=T;

    //printf("Before lines %e\n",fT);

#ifdef ALTERNATE_LINE_EMISSION


    //find integral over pass band in observed frame
    double Tl  = 0;
    double T_norm = 0;
    for(m=0;m<n_pb-1;m++)
    {
      nu_lo = c/x_pb[m+1];
      nu_hi = c/x_pb[m];
      nu_ave = 0.5*(nu_hi + nu_lo);

      //only worry about the transmissivity above 0.001
      //when estimating the average transmissivity
      if(y_pb[m]>1.0e-3)
      {
        Tl  += 0.5*(y_pb[m+1] + y_pb[m])*(nu_hi-nu_lo)/nu_ave;
        T_norm  += (nu_hi-nu_lo)/nu_ave;
      }

    }

    //first compute f_lambda HERE
    double dlambda_pb;
    int il;
    dlambda_pb = (dlam_hi + dlam_lo)*(1+z); //band width in Ang, in observed frame
    double lambda_line_obs;
    double fac_pb;
    //printf("n_emission_lines %d\n",n_emission_lines);
    for(il=0;il<n_emission_lines;il++)
    {
      //find observed wavelength of line
      lambda_line_obs = lambda_emission_lines[il] * (1+z);

      //printf("llo %e\n",lambda_line_obs);

      //check whether line falls in this pass band
      if((lambda_line_obs>=x_pb[0])&&(lambda_line_obs<=x_pb[n_pb-1]) )
      {

        //this line falls in the pass band

        //find f_lambda in Jy * Hz/Ang
        y_T = f_nu_emission_lines[il] / dlambda_pb;

        //moderate by the pass band
        i_lo = gsl_interp_bsearch(x_pb,lambda_line_obs,0,n_pb);
        fac_pb = y_pb[i_lo] + (lambda_line_obs - x_pb[i_lo])*(y_pb[i_lo+1] - y_pb[i_lo])/(x_pb[i_lo+1] - x_pb[i_lo]);
        fac_pb /= (Tl/T_norm);

        y_T *= pow(lambda_line_obs,2) / c; // to uJy in observed frame


        //printf("BANDWIDTH lambda_lines[%d] %e y_T %e\n",il,lambda_lines[il],y_T*fac_pb);
        //apply dust 
        //printf("lambda_emission_lines %e lambda_line %e\n",lambda_emission_lines[il],lambda_line[il]);
        kp = klambda_smc(lambda_emission_lines[il] * 1.0e-4);
        if(lambda_emission_lines[il]  * 1.0e-4 < lyman_limit * 1.0e-4)
          kp = 0.0;
        AEBV = -0.4*ebv*kp;

        //if( 0.5*(y_pb[i_lo]+y_pb[i_lo+1])>1.0e-3)
        fT +=  y_T * fac_pb * (1 - f_esc) * pow(10.0, AEBV);

        //#error check factors of 1+z
      }
    }

#endif //ALTERNATE_LINE_EMISSION

    //f_model[l] = fT/T;
    f_model[l] = fT;


    //now we need to add the contribution from the lines
    //to the pass band

    fprintf(fpsed,"%e\t%e\n",lambda_model[l],f_model[l]);

    free(x_pb_o);

    //printf("After lines %e\n",fT);
  }
  fclose(fpsed);


  //double line_strength =     
  double j_nu_hbeta = 4.78e-13*(1-f_esc)*pow(10.0,Nlyc); //eqn 2 of ono's 2010 paper

  j_nu_hbeta *= areal_factor(z);
#ifdef FIT_LINE_STRENGTHS
  j_nu_hbeta *= B_line;
#endif //FIT_LINE_STRENGTHS

  fpsed = fopen("sed.lines.txt","w");
  fprintf(fpsed,"%e\n",Nlyc);
  fprintf(fpsed,"%e\n",z);
  fprintf(fpsed,"%e\n",f_esc);
  fprintf(fpsed,"%e\n",areal_factor(z));
  fprintf(fpsed,"%e\n",j_nu_hbeta);
  fprintf(fpsed,"%e\n",j_nu_hbeta*F_OIII_4959((0.1*0.02)));
  fprintf(fpsed,"%e\n",j_nu_hbeta*F_OIII_5007((0.1*0.02)));

  printf("Hbeta = %e\n",j_nu_hbeta);
  printf("F_OIII_4959 = %e fac %e\n",j_nu_hbeta*F_OIII_4959((0.1*0.02)),F_OIII_4959(0.1*0.02));
  printf("F_OIII_5007 = %e fac %e\n",j_nu_hbeta*F_OIII_5007((0.1*0.02)),F_OIII_5007(0.1*0.02));

  fclose(fpsed);


  free(f_nu_tot);

  free(f_nu_sed_lo);
  free(f_nu_sed_hi);
  free(f_nu_line_lo);
  free(f_nu_line_hi);
  free(f_nu_cont_lo);
  free(f_nu_cont_hi);

	return 0;
}
