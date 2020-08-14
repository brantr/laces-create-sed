#include "pass_bands.h"
#include <stdio.h>
#include <stdlib.h>

int n_pb;      //number of samples on the curve
double *x_pb;  //wavelength
double *y_pb;  //transmissivity


//load the pass band from a file
void LoadPassBand(int l)
{
  FILE *fp;
  char fname[200];
  int i;
  double xb, yb;
  double fac = 1.0;

  //select the pass band to read in
  switch(l)
  {
    case 0:
      sprintf(fname,"pass_bands/f275w_trans.dat");
      break;
    case 1:
      sprintf(fname,"pass_bands/f336w_trans.dat");
      break;
    case 2:
      sprintf(fname,"pass_bands/f435w_trans.dat");
      break;
    case 3:
      sprintf(fname,"pass_bands/f606w_trans.dat");
      break;
    case 4:
      sprintf(fname,"pass_bands/f775w_trans.dat");
      break;
    case 5:
      sprintf(fname,"pass_bands/f814w_trans.dat");
      break;
    case 6:
      sprintf(fname,"pass_bands/f850lp_trans.dat");
      break;
    case 7:
      sprintf(fname,"pass_bands/f105w_trans.dat");
      break;
    case 8:
      sprintf(fname,"pass_bands/f125w_trans.dat");
      break;
    case 9:
      sprintf(fname,"pass_bands/f140w_trans.dat");
      break;
    case 10:
      sprintf(fname,"pass_bands/f160w_trans.dat");
      break;
  }

  //load the passband transmissivity curve
  if(!(fp = fopen(fname,"r")))
  {
    printf("Error loading file %s.\n",fname);
    exit(0);
  }
  fscanf(fp,"%d\n",&n_pb);

  //allocate the pass band array
  x_pb = (double *) calloc(n_pb,sizeof(double));
  y_pb = (double *) calloc(n_pb,sizeof(double));

  for(i=0;i<n_pb;i++)
  {
    fscanf(fp,"%lf %lf\n",&xb,&yb);
    x_pb[i] = xb*fac;
    y_pb[i] = yb;
  }

  //close the passband file
  fclose(fp);
}
void FreePassBand(void)
{
  free(x_pb);
  free(y_pb);
  n_pb = 0;
}