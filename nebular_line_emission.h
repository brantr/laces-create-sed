/*! \file nebular_line_emission.h
 *  \brief Function declarations for routines that provide the nebular line emission
 *  from H, HeI, HeII, and metals. */
#ifndef NEBULAR_LINE_EMISSION
#define NEBULAR_LINE_EMISSION

/*! \fn double *line_spectrum(int *n_lambda, double lambda_min, double lambda_max, double *&lambda_lines)
 *  \brief Allocates an empty array to use as a line spectrum. */
double *line_spectrum(int *n_lambda, double lambda_min, double lambda_max, double *&lambda_lines);

/*! \fn double line_width_factor(double m, double T)
 *  \brief Returns the fractional line width in delta_nu/nu.  Either set by
 *  FIXED_LINE_WIDTH or the mass of the ion. */
double line_width_factor(double m, double T);

/*! \fn double line_width_factor_H(double T)
 *  \brief Thermal line width for H */
double line_width_factor_H(double T);

/*! \fn double line_width_factor_He(double T)
 *  \brief Thermal line width for He */
double line_width_factor_He(double T);

/*! \fn double line_width_factor_Ar(double T)
 *  \brief Thermal line width for Ar */
double line_width_factor_Ar(double T);

/*! \fn double line_width_factor_C(double T)
 *  \brief Thermal line width for C */
double line_width_factor_C(double T);

/*! \fn double line_width_factor_Mg(double T)
 *  \brief Thermal line width for Mg */
double line_width_factor_Mg(double T);

/*! \fn double line_width_factor_N(double T)
 *  \brief Thermal line width for N */
double line_width_factor_N(double T);

/*! \fn double line_width_factor_Ne(double T)
 *  \brief Thermal line width for Ne */
double line_width_factor_Ne(double T);

/*! \fn double line_width_factor_O(double T)
 *  \brief Thermal line width for O */
double line_width_factor_O(double T);

/*! \fn double line_width_factor_S(double T)
 *  \brief Thermal line width for S */
double line_width_factor_S(double T);

/*! \fn double line_profile(double lambda, double lambda_c, double lambda_width)
 *  \brief Gaussian line profile for a naturally broadened line. */
double line_profile(double lambda, double lambda_c, double lambda_width);

/*! \fn double line_profile_nu(double nu, double nu_c, double nu_width)
 *  \brief Gaussian line profile for a naturally broadened line. */
double line_profile_nu(double nu, double nu_c, double nu_width);


/*! \fn void add_line_to_spectrum(double lambda, double lambda_width, double f, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Insers a single line into the spectrum. */
void add_line_to_spectrum(double lambda, double lambda_width, double f, double *&lambda_lines, double *&spectrum_lines, int *n_lambda);

/*! \fn void add_HI_spectrum(double F_HI, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Adds hydrogen recombination line spectrum, with a strength
 *  	  scaled to F_HI = F_H_beta.
 */
void add_HI_spectrum(double F_HI, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda);
/*! \fn void add_HeI_spectrum(double F_Hbeta, double n_Hep, double n_p, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Adds HeI recombination line spectrum, with a strength
 *  	  scaled to F_HI = F_H_beta.
 */
void add_HeI_spectrum(double F_Hbeta, double n_Hep, double n_p, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda);
/*! \fn void add_HeII_spectrum(double F_Hbeta, double n_Hepp, double n_p, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Adds HeII recombination line spectrum, with a strength
 *  	  scaled to F_HI = F_H_beta.
 */
void add_HeII_spectrum(double F_Hbeta, double n_Hepp, double n_p, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda);
/*! \fn void add_metal_spectrum(double F_Hbeta, double Z, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda)
 *  \brief Adds metal recombination line spectrum, with a strength
 *  	  scaled to F_HI = F_H_beta */
void add_metal_spectrum(double F_HI, double Z, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda);

/*! \fn void add_all_spectra(double F_HI, double Z, double n_p, double n_Hep, double n_Hepp, double T, double *&lambdas, double *&spectrums, int *n_lambda)
 *  \brief Adds all nebular lines to the input spectrum, with a 
 *  strength based relative F_HI = F_H_beta.
 */
void add_all_spectra(double F_HI, double Z, double n_p, double n_Hep, double n_Hepp, double T, double *&lambda_lines, double *&spectrum_lines, int *n_lambda);


//metal lines from P. Anders & 
//U. Fritze-v. Alvensleben 2003, AA, 401, 1063-1070


/*! \var double Z_line[5]
 *  \brief Array of metallicities at which the metal line strengths are calculated.*/
extern double Z_line[5];

/*! \fn double F_ArIII_7136(double Z)
 *  \brief Strength of ArIII 7136 line relative to H-beta. */
double F_ArIII_7136(double Z);

/*! \fn double F_ArIII_7136(double Z)
 *  \brief Strength of ArIII 7136 line relative to H-beta. */
double F_ArIII_7751(double Z);

/*! \fn double F_ArIV_HeI_4711(double Z)
 *  \brief Strength of ArIV and HeI 4711 line relative to H-beta. */
double F_ArIV_HeI_4711(double Z);

/*! \fn double F_CII_1335(double Z)
 *  \brief Strength of CII 1335 line relative to H-beta. */
double F_CII_1335(double Z);

/*! \fn double F_CII_2326(double Z)
 *  \brief Strength of CII 2336 line relative to H-beta. */
double F_CII_2326(double Z);

/*! \fn double F_CIII_1909(double Z)
 *  \brief Strength of CIII 1909 line relative to H-beta. */
double F_CIII_1909(double Z);

/*! \fn double F_MgII_2798(double Z)
 *  \brief Strength of MgII 2798 line relative to H-beta. */
double F_MgII_2798(double Z);

/*! \fn double F_NI_5199(double Z)
 *  \brief Strength of NI 5199 line relative to H-beta. */
double F_NI_5199(double Z);

/*! \fn double F_NII_2141(double Z)
 *  \brief Strength of NII 2141 line relative to H-beta. */
double F_NII_2141(double Z);

/*! \fn double F_NII_5755(double Z)
 *  \brief Strength of NII 5755 line relative to H-beta. */
double F_NII_5755(double Z);

/*! \fn double F_NII_5755(double Z)
 *  \brief Strength of NII 5755 line relative to H-beta. */
double F_NII_6548(double Z);

/*! \fn double F_NII_5755(double Z)
 *  \brief Strength of NII 5755 line relative to H-beta. */
double F_NII_6583(double Z);

/*! \fn double F_NeIII_3869(double Z)
 *  \brief Strength of NeIII 3869 line relative to H-beta. */
double F_NeIII_3869(double Z);

/*! \fn double F_OII_3727(double Z)
 *  \brief Strength of OII 3727 line relative to H-beta. */
double F_OII_3727(double Z);

/*! \fn double F_OIII_1663(double Z)
 *  \brief Strength of OIII 1663 line relative to H-beta. */
double F_OIII_1663(double Z);

/*! \fn double F_OIII_4363(double Z)
 *  \brief Strength of OIII 4363 line relative to H-beta. */
double F_OIII_4363(double Z);

/*! \fn double F_OIII_4959(double Z)
 *  \brief Strength of OIII 4959 line relative to H-beta. */
double F_OIII_4959(double Z);

/*! \fn double F_OIII_5007(double Z)
 *  \brief Strength of OIII 5007 line relative to H-beta. */
double F_OIII_5007(double Z);

/*! \fn double F_OII_3727(double Z)
 *  \brief Strength of OII 3727 line relative to H-beta. */
double F_OI_6300(double Z);

/*! \fn double F_SIII(double Z)
 *  \brief Strength of SIII 6312 line relative to H-beta. */
double F_SIII_6312(double Z);

/*! \fn double F_OII_3727(double Z)
 *  \brief Strength of OII 3727 line relative to H-beta. */
double F_OII_7320(double Z);

/*! \fn double F_OII_3727(double Z)
 *  \brief Strength of OII 3727 line relative to H-beta. */
double F_OII_7331(double Z);

/*! \fn double F_SII_4069(double Z)
 *  \brief Strength of SII 4069 line relative to H-beta. */
double F_SII_4069(double Z);

/*! \fn double F_SII_4076(double Z)
 *  \brief Strength of SII 4076 line relative to H-beta. */
double F_SII_4076(double Z);

/*! \fn double F_SII_4076(double Z)
 *  \brief Strength of SII 4076 line relative to H-beta. */
double F_SIII_6312(double Z);

/*! \fn double F_SII_6716(double Z)
 *  \brief Strength of SII 6716 line relative to H-beta. */
double F_SII_6716(double Z);

/*! \fn double F_SII_6730(double Z)
 *  \brief Strength of SII 6730 line relative to H-beta. */
double F_SII_6730(double Z);

/*! \fn double F_SIII_9069(double Z)
 *  \brief Strength of SIII 9069 line relative to H-beta. */
double F_SIII_9069(double Z);

/*! \fn double F_SIII_9531(double Z)
 *  \brief Strength of SIII 9531 line relative to H-beta. */
double F_SIII_9531(double Z);

/*! \fn double F_SII_4076(double Z)
 *  \brief Strength of SII 4076 line relative to H-beta. */
double F_SII_10287(double Z);

/*! \fn double F_SII_4076(double Z)
 *  \brief Strength of SII 4076 line relative to H-beta. */
double F_SII_10320(double Z);

/*! \fn double F_SII_4076(double Z)
 *  \brief Strength of SII 4076 line relative to H-beta. */
double F_SII_10336(double Z);

/*! \fn double F_H_zeta_HeII_3889(double Z)
 *  \brief Strength of Hzeta+HeII 3889 line relative to H-beta. */
double F_H_zeta_HeII_3889(double Z);

/*! \fn double F_H_epsilon_HeII_3970(double Z)
 *  \brief Strength of Hepsilon+HeII 3970 line relative to H-beta. */
double F_H_epsilon_HeII_3970(double Z);

/*! \fn double F_HeI_4026(double Z)
 *  \brief Strength of HeI 4026 line relative to H-beta. */
double F_HeI_4026(double Z);
//double F_HeI_4471(double Z);

/*! \fn double F_HeI_5876(double Z)
 *  \brief Strength of HeI 5876 line relative to H-beta. */
double F_HeI_5876(double Z);

/*! \fn double F_HeI_7065(double Z)
 *  \brief Strength of HeI 7065 line relative to H-beta. */
double F_HeI_7065(double Z);


/*! \fn double lambda_CII_1335(double Z)
 *  \brief Wavelength of CII 1335 line in microns.*/
double lambda_CII_1335(void);

/*! \fn double lambda_OIII_1663(double Z)
 *  \brief Wavelength of OIII 1663 line in microns.*/
double lambda_OIII_1663(void);

/*! \fn double lambda_CIII_1909(double Z)
 *  \brief Wavelength of CIII 1909 line in microns.*/
double lambda_CIII_1909(void);

/*! \fn double lambda_NII_2141(double Z)
 *  \brief Wavelength of NII 2141 line in microns.*/
double lambda_NII_2141(void);


/*! \fn double lambda_CII_2326(double Z)
 *  \brief Wavelength of CII 2326 line in microns.*/
double lambda_CII_2326(void);

/*! \fn double lambda_MgII_2798(double Z)
 *  \brief Wavelength of MgII 2798 line in microns.*/
double lambda_MgII_2798(void);

/*! \fn double lambda_OII_3727(double Z)
 *  \brief Wavelength of OII 3727 line in microns.*/
double lambda_OII_3727(void);

/*! \fn double lambda_NeIII_3869(double Z)
 *  \brief Wavelength of NeIII 3869 line in microns.*/
double lambda_NeIII_3869(void);

/*! \fn double lambda_H_zeta_HeII_3889(double Z)
 *  \brief Wavelength of Hzeta and HeII 3889 line in microns.*/
double lambda_H_zeta_HeII_3889(void);

/*! \fn double lambda_H_epsilon_HeII_3970(double Z)
 *  \brief Wavelength of Hepsilon and HeII 3970 line in microns.*/
double lambda_H_epsilon_HeII_3970(void);

/*! \fn double lambda_SII_4069(double Z)
 *  \brief Wavelength of SII 4069 line in microns.*/
double lambda_HeI_4026(void);

/*! \fn double lambda_SII_4069(double Z)
 *  \brief Wavelength of SII 4069 line in microns.*/
double lambda_SII_4069(void);

/*! \fn double lambda_SII_4076(double Z)
 *  \brief Wavelength of SII 4076 line in microns.*/
double lambda_SII_4076(void);

/*! \fn double lambda_OIII_4363(double Z)
 *  \brief Wavelength of OIII 4363 line in microns.*/
double lambda_OIII_4363(void);

/*! \fn double lambda_ArIV_HeI_4711(double Z)
 *  \brief Wavelength of ArIV and HeI 4711 line in microns.*/
double lambda_ArIV_HeI_4711(void);

/*! \fn double lambda_OIII_4959(double Z)
 *  \brief Wavelength of OIII 4959 line in microns.*/
double lambda_OIII_4959(void);

/*! \fn double lambda_OIII_5007(double Z)
 *  \brief Wavelength of OIII 5007 line in microns.*/
double lambda_OIII_5007(void);


/*! \fn double lambda_NI_5199(double Z)
 *  \brief Wavelength of NI 5199 line in microns.*/
double lambda_NI_5199(void);

/*! \fn double lambda_NI_5755(double Z)
 *  \brief Wavelength of NI 5755 line in microns.*/
double lambda_NII_5755(void);

/*! \fn double lambda_OIII_5007(double Z)
 *  \brief Wavelength of OIII 5007 line in microns.*/
double lambda_OI_6300(void);

/*! \fn double lambda_SIII_6312(double Z)
 *  \brief Wavelength of SIII 6312 line in microns.*/
double lambda_SIII_6312(void);

/*! \fn double lambda_NI_5755(double Z)
 *  \brief Wavelength of NI 5755 line in microns.*/
double lambda_NII_6548(void);

/*! \fn double lambda_NI_5755(double Z)
 *  \brief Wavelength of NI 5755 line in microns.*/
double lambda_NII_6583(void);

/*! \fn double lambda_SII_6716(double Z)
 *  \brief Wavelength of SII 6716 line in microns.*/
double lambda_SII_6716(void);

/*! \fn double lambda_SII_6730(double Z)
 *  \brief Wavelength of SII 6730 line in microns.*/
double lambda_SII_6730(void);

//double lambda_HeI_7065(void);

/*! \fn double lambda_ArIII_7136(double Z)
 *  \brief Wavelength of ArIII 7136 line in microns.*/
double lambda_ArIII_7136(void);

/*! \fn double lambda_OIII_5007(double Z)
 *  \brief Wavelength of OIII 5007 line in microns.*/
double lambda_OII_7320(void);

/*! \fn double lambda_OIII_5007(double Z)
 *  \brief Wavelength of OIII 5007 line in microns.*/
double lambda_OII_7331(void);

/*! \fn double lambda_ArIII_7751(double Z)
 *  \brief Wavelength of ArIII 7751 line in microns.*/
double lambda_ArIII_7751(void);

/*! \fn double lambda_SIII_9531(double Z)
 *  \brief Wavelength of SIII 9531 line in microns.*/
double lambda_SIII_6312(void);

/*! \fn double lambda_SIII_9069(double Z)
 *  \brief Wavelength of SIII 9069 line in microns.*/
double lambda_SIII_9069(void);

/*! \fn double lambda_SIII_9531(double Z)
 *  \brief Wavelength of SIII 9531 line in microns.*/
double lambda_SIII_9531(void);

/*! \fn double lambda_SIII_9531(double Z)
 *  \brief Wavelength of SIII 9531 line in microns.*/
double lambda_SII_10287(void);

/*! \fn double lambda_SIII_9531(double Z)
 *  \brief Wavelength of SIII 9531 line in microns.*/
double lambda_SII_10320(void);

/*! \fn double lambda_SIII_9531(double Z)
 *  \brief Wavelength of SIII 9531 line in microns.*/
double lambda_SII_10336(void);

/*! \fn double lambda_H(int n1, int n2)
 *  \brief Wavelength of a transition in the H atom between
 *  orbital states n1 and n2. The Rydbeg formula.*/
double lambda_H(int n1, int n2);

//HI recombination lines at ne=100, T=10^4, Case B, Table 4.4 of Osterbrock & Ferland

/*! \fn double lambda_H_alpha(void)
 *  \brief Wavelength of H-alpha */
double lambda_H_alpha(void);

/*! \fn double lambda_H_beta(void)
 *  \brief Wavelength of H-beta */
double lambda_H_beta(void);

/*! \fn double lambda_H_gamma(void)
 *  \brief Wavelength of H-gamma*/
double lambda_H_gamma(void);

/*! \fn double lambda_H_delta(void)
 *  \brief Wavelength of H-delta*/
double lambda_H_delta(void);

/*! \fn double lambda_H_10(void)
 *  \brief Wavelength of H-10*/
double lambda_H_10(void);

/*! \fn double lambda_H_15(void)
 *  \brief Wavelength of H-15*/
double lambda_H_15(void);

/*! \fn double lambda_H_20(void)
 *  \brief Wavelength of H-20*/
double lambda_H_20(void);


/*! \fn double lambda_H_Palpha(void)
 *  \brief Wavelength of Paschen-alpha line in microns. */
double lambda_H_Palpha(void);

/*! \fn double lambda_H_Pbeta(void)
 *  \brief Wavelength of Paschen-beta line in microns. */
double lambda_H_Pbeta(void);

/*! \fn double lambda_H_Pgamma(void)
 *  \brief Wavelength of Paschen-gamma line in microns. */
double lambda_H_Pgamma(void);

/*! \fn double lambda_H_Pdelta(void)
 *  \brief Wavelength of Paschen-delta line in microns. */
double lambda_H_Pdelta(void);

/*! \fn double lambda_H_P10(void)
 *  \brief Wavelength of Paschen-10 line in microns. */
double lambda_H_P10(void);

/*! \fn double lambda_H_P15(void)
 *  \brief Wavelength of Paschen-15 line in microns. */
double lambda_H_P15(void);

/*! \fn double lambda_H_P20(void)
 *  \brief Wavelength of Paschen-20 line in microns. */
double lambda_H_P20(void);

/*! \fn double lambda_H_Balpha(void)
 *  \brief Wavelength of H Brackett-alpha in microns. */
double lambda_H_Balpha(void);

/*! \fn double lambda_H_Bbeta(void)
 *  \brief Wavelength of H Brackett-beta in microns. */
double lambda_H_Bbeta(void);

/*! \fn double lambda_H_Bgamma(void)
 *  \brief Wavelength of H Brackett-gamma in microns. */
double lambda_H_Bgamma(void);

/*! \fn double lambda_H_Bdelta(void)
 *  \brief Wavelength of H Brackett-delta in microns. */
double lambda_H_Bdelta(void);

/*! \fn double lambda_H_B10(void)
 *  \brief Wavelength of H Brackett-10 in microns. */
double lambda_H_B10(void);

/*! \fn double lambda_H_B15(void)
 *  \brief Wavelength of H Brackett-15 in microns. */
double lambda_H_B15(void);

/*! \fn double lambda_H_B20(void)
 *  \brief Wavelength of H Brackett-20 in microns. */
double lambda_H_B20(void);

/*! \fn double F_H_alpha(void)
 *  \brief H-alpha line emission strength relative to H-beta */
double F_H_alpha(void);
/*! \fn double F_H_beta(void)
 *  \brief H-beta line emission strength relative to H-beta */
double F_H_beta(void);

/*! \fn double F_H_gamma(void)
 *  \brief H-gamma line emission strength relative to H-beta */
double F_H_gamma(void);

/*! \fn double F_H_delta(void)
 *  \brief H-delta line emission strength relative to H-beta */
double F_H_delta(void);

/*! \fn double F_H_10(void)
 *  \brief H-10 line emission strength relative to H-beta */
double F_H_10(void);

/*! \fn double F_H_15(void)
 *  \brief H-15 line emission strength relative to H-beta */
double F_H_15(void);

/*! \fn double F_H_20(void)
 *  \brief H-20 line emission strength relative to H-beta */
double F_H_20(void);

/*! \fn double F_H_Palpha(void)
 *  \brief Emission strength of Paschen-alpha relative to H-beta */
double F_H_Palpha(void);

/*! \fn double F_H_Pbeta(void)
 *  \brief Emission strength of Paschen-beta relative to H-beta */
double F_H_Pbeta(void);

/*! \fn double F_H_Pgamma(void)
 *  \brief Emission strength of Paschen-gamma relative to H-beta */
double F_H_Pgamma(void);

/*! \fn double F_H_Pdelta(void)
 *  \brief Emission strength of Paschen-delta relative to H-beta */
double F_H_Pdelta(void);

/*! \fn double F_H_P10(void)
 *  \brief Emission strength of Paschen-10 relative to H-beta */
double F_H_P10(void);

/*! \fn double F_H_P15(void)
 *  \brief Emission strength of Paschen-15 relative to H-beta */
double F_H_P15(void);

/*! \fn double F_H_P20(void)
 *  \brief Emission strength of Paschen-20 relative to H-beta */
double F_H_P20(void);

/*! \fn double F_H_Balpha(void)
 *  \brief Emission strength of Brackett-alpha relative to H-beta */
double F_H_Balpha(void);

/*! \fn double F_H_Bbeta(void)
 *  \brief Emission strength of Brackett-beta relative to H-beta */
double F_H_Bbeta(void);

/*! \fn double F_H_Bgamma(void)
 *  \brief Emission strength of Brackett-gamma relative to H-beta */
double F_H_Bgamma(void);

/*! \fn double F_H_Bdelta(void)
 *  \brief Emission strength of Brackett-delta relative to H-beta */
double F_H_Bdelta(void);

/*! \fn double F_H_B10(void)
 *  \brief Emission strength of Brackett-10 relative to H-beta */
double F_H_B10(void);

/*! \fn double F_H_B15(void)
 *  \brief Emission strength of Brackett-15 relative to H-beta */
double F_H_B15(void);

/*! \fn double F_H_B20(void)
 *  \brief Emission strength of Brackett-20 relative to H-beta */
double F_H_B20(void);

/*! \fn double lambda_He(int n1, int n2)
 *  \brief Wavelength of a transition in the He atom between
 *  orbital states n1 and n2. The Z-dependent Rydbeg formula.*/
double lambda_HeII(int n1, int n2);

//HeII recombination lines at ne=100, T=10^4, Case B, Table 4.5 of Osterbrock & Ferland

/*! \fn double lambda_HeII_32(void)
 *  \brief Wavelength of HeII 3->2 line, in microns. */
double lambda_HeII_32(void);

/*! \fn double lambda_HeII_42(void)
 *  \brief Wavelength of HeII 4->2 line, in microns. */
double lambda_HeII_42(void);

/*! \fn double lambda_HeII_52(void)
 *  \brief Wavelength of HeII 5->2 line, in microns. */
double lambda_HeII_52(void);

/*! \fn double lambda_HeII_62(void)
 *  \brief Wavelength of HeII 6->2 line, in microns. */
double lambda_HeII_62(void);

/*! \fn double lambda_HeII_102(void)
 *  \brief Wavelength of HeII 10->2 line, in microns. */
double lambda_HeII_102(void);

/*! \fn double lambda_HeII_152(void)
 *  \brief Wavelength of HeII 15->2 line, in microns. */
double lambda_HeII_152(void);

/*! \fn double lambda_HeII_202(void)
 *  \brief Wavelength of HeII 20->2 line, in microns. */
double lambda_HeII_202(void);

/*! \fn double F_HeII_norm(double n_Hepp, double n_p)
 *  \brief Normalization of the HeII 4686 (n=4 -> n=3) line strength
 *  relative to H-beta */
double F_HeII_norm(double n_Hepp, double n_p);

/*! \fn double F_HeII_4686(void)
 *  \brief Line strength relative to HeII 4686.*/
double F_HeII_4686(void);

/*! \fn double F_HeII_32(void)
 *  \brief Line strength of HeII 3->2, relative to HeII 4->3.*/
double F_HeII_32(void);

/*! \fn double F_HeII_42(void)
 *  \brief Line strength of HeII 4->2, relative to HeII 4->3.*/
double F_HeII_42(void);

/*! \fn double F_HeII_52(void)
 *  \brief Line strength of HeII 5->2, relative to HeII 4->3.*/
double F_HeII_52(void);

/*! \fn double F_HeII_62(void)
 *  \brief Line strength of HeII 6->2, relative to HeII 4->3.*/
double F_HeII_62(void);

/*! \fn double F_HeII_102(void)
 *  \brief Line strength of HeII 10->2, relative to HeII 4->3.*/
double F_HeII_102(void);

/*! \fn double F_HeII_152(void)
 *  \brief Line strength of HeII 15->2, relative to HeII 4->3.*/
double F_HeII_152(void);

/*! \fn double F_HeII_202(void)
 *  \brief Line strength of HeII 20->2, relative to HeII 4->3.*/
double F_HeII_202(void);

/*! \fn double lambda_HeII_43(void)
 *  \brief Wavelength of HeII 4->3 line, in microns --- 4686. */
double lambda_HeII_43(void);

/*! \fn double lambda_HeII_53(void)
 *  \brief Wavelength of HeII 5->3 line, in microns. */
double lambda_HeII_53(void);

/*! \fn double lambda_HeII_63(void)
 *  \brief Wavelength of HeII 6->3 line, in microns. */
double lambda_HeII_63(void);

/*! \fn double lambda_HeII_73(void)
 *  \brief Wavelength of HeII 7->3 line, in microns. */
double lambda_HeII_73(void);

/*! \fn double lambda_HeII_103(void)
 *  \brief Wavelength of HeII 10->3 line, in microns. */
double lambda_HeII_103(void);

/*! \fn double lambda_HeII_153(void)
 *  \brief Wavelength of HeII 15->3 line, in microns. */
double lambda_HeII_153(void);

/*! \fn double lambda_HeII_203(void)
 *  \brief Wavelength of HeII 20->3 line, in microns. */
double lambda_HeII_203(void);


/*! \fn double F_HeII_43(void)
 *  \brief Line strength relative to HeII 4686.*/
double F_HeII_43(void);

/*! \fn double F_HeII_53(void)
 *  \brief Line strength of HeII 5->3, relative to HeII 4->3.*/
double F_HeII_53(void);

/*! \fn double F_HeII_63(void)
 *  \brief Line strength of HeII 6->3, relative to HeII 4->3.*/
double F_HeII_63(void);

/*! \fn double F_HeII_73(void)
 *  \brief Line strength of HeII 7->3, relative to HeII 4->3.*/
double F_HeII_73(void);

/*! \fn double F_HeII_103(void)
 *  \brief Line strength of HeII 10->3, relative to HeII 4->3.*/
double F_HeII_103(void);

/*! \fn double F_HeII_103(void)
 *  \brief Line strength of HeII 10->3, relative to HeII 4->3.*/
double F_HeII_153(void);

/*! \fn double F_HeII_203(void)
 *  \brief Line strength of HeII 20->3, relative to HeII 4->3.*/
double F_HeII_203(void);

/*! \fn double F_HeII_54(void)
 *  \brief Line strength of HeII 5->4, relative to HeII 4->3.*/
double lambda_HeII_54(void);

/*! \fn double lambda_HeII_64(void)
 *  \brief Wavelength of HeII 6->4 line, in microns. */
double lambda_HeII_64(void);

/*! \fn double lambda_HeII_74(void)
 *  \brief Wavelength of HeII 7->4 line, in microns. */
double lambda_HeII_74(void);

/*! \fn double lambda_HeII_84(void)
 *  \brief Wavelength of HeII 8->4 line, in microns. */
double lambda_HeII_84(void);

/*! \fn double lambda_HeII_104(void)
 *  \brief Wavelength of HeII 10->4 line, in microns. */
double lambda_HeII_104(void);

/*! \fn double lambda_HeII_154(void)
 *  \brief Wavelength of HeII 15->4 line, in microns. */
double lambda_HeII_154(void);

/*! \fn double lambda_HeII_204(void)
 *  \brief Wavelength of HeII 20->4 line, in microns. */
double lambda_HeII_204(void);

/*! \fn double F_HeII_54(void)
 *  \brief Line strength of HeII 5->4, relative to HeII 4->3.*/
double F_HeII_54(void);

/*! \fn double F_HeII_64(void)
 *  \brief Line strength of HeII 6->4, relative to HeII 4->3.*/
double F_HeII_64(void);
/*! \fn double F_HeII_74(void)
 *  \brief Line strength of HeII 7->4, relative to HeII 4->3.*/
double F_HeII_74(void);
/*! \fn double F_HeII_84(void)
 *  \brief Line strength of HeII 8->4, relative to HeII 4->3.*/
double F_HeII_84(void);
/*! \fn double F_HeII_104(void)
 *  \brief Line strength of HeII 10->4, relative to HeII 4->3.*/
double F_HeII_104(void);
/*! \fn double F_HeII_154(void)
 *  \brief Line strength of HeII 15->4, relative to HeII 4->3.*/
double F_HeII_154(void);
/*! \fn double F_HeII_204(void)
 *  \brief Line strength of HeII 20->4, relative to HeII 4->3.*/
double F_HeII_204(void);


//HeI Case B lines, ne=100, T = 10^4, table 4.6 of Osterbrock & Ferland

/*! \fn double lambda_HeI_4471(void)
 *  \brief Wavelength of HeI 4471 in microns. */
double lambda_HeI_4471(void);

/*! \fn double lambda_HeI_5876(void)
 *  \brief Wavelength of HeI 5876 in microns. */
double lambda_HeI_5876(void);

/*! \fn double lambda_HeI_4026(void)
 *  \brief Wavelength of HeI 4026 in microns. */
double lambda_HeI_4026(void);

/*! \fn double lambda_HeI_7065(void)
 *  \brief Wavelength of HeI 7065 in microns. */
double lambda_HeI_7065(void);

/*! \fn double lambda_HeI_10830(void)
 *  \brief Wavelength of HeI 10830 in microns. */
double lambda_HeI_10830(void);

/*! \fn double lambda_HeI_3889(void)
 *  \brief Wavelength of HeI 3889 in microns. */
double lambda_HeI_3889(void);

/*! \fn double lambda_HeI_3187(void)
 *  \brief Wavelength of HeI 3187 in microns. */
double lambda_HeI_3187(void);

/*! \fn double lambda_HeI_6678(void)
 *  \brief Wavelength of HeI 6678 in microns. */
double lambda_HeI_6678(void);

/*! \fn double lambda_HeI_4922(void)
 *  \brief Wavelength of HeI 4922 in microns. */
double lambda_HeI_4922(void);

/*! \fn double lambda_HeI_5016(void)
 *  \brief Wavelength of HeI 5016 in microns. */
double lambda_HeI_5016(void);

/*! \fn double lambda_HeI_3965(void)
 *  \brief Wavelength of HeI 3965 in microns. */
double lambda_HeI_3965(void);

/*! \fn double F_HeI_norm(double n_Hep, double n_p)
 *  \brief Normalization of HeII 4471 line strength relative to H-beta. */
double F_HeI_norm(double n_Hep, double n_p);

/*! \fn double F_HeI_4471(void)
 *  \brief Normalization of HeI 4471 line strength relative to HeI 4471. */
double F_HeI_4471(void);

/*! \fn double F_HeI_5876(void)
 *  \brief Normalization of HeI 5876 line strength relative to HeI 4471. */
double F_HeI_5876(void);

/*! \fn double F_HeI_4026(void)
 *  \brief Normalization of HeI 4026 line strength relative to HeI 4471. */
double F_HeI_4026(void);

/*! \fn double F_HeI_7065(void)
 *  \brief Normalization of HeI 7065 line strength relative to HeI 4471. */
double F_HeI_7065(void);

/*! \fn double F_HeI_10830(void)
 *  \brief Normalization of HeI 10830 line strength relative to HeI 4471. */
double F_HeI_10830(void);

/*! \fn double F_HeI_3889(void)
 *  \brief Normalization of HeI 3889 line strength relative to HeI 4471. */
double F_HeI_3889(void);

/*! \fn double F_HeI_3187(void)
 *  \brief Normalization of HeI 3187 line strength relative to HeI 4471. */
double F_HeI_3187(void);

/*! \fn double F_HeI_6678(void)
 *  \brief Normalization of HeI 6678 line strength relative to HeI 4471. */
double F_HeI_6678(void);

/*! \fn double F_HeI_4922(void)
 *  \brief Normalization of HeI 4922 line strength relative to HeI 4471. */
double F_HeI_4922(void);

/*! \fn double F_HeI_5016(void)
 *  \brief Normalization of HeI 5016 line strength relative to HeI 4471. */
double F_HeI_5016(void);

/*! \fn double F_HeI_3965(void)
 *  \brief Normalization of HeI 3965 line strength relative to HeI 4471. */
double F_HeI_3965(void);




#endif //NEBULAR_LINE_EMISSION
