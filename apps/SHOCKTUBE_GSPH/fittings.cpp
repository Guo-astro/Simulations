/*
 * fittings.h
 *
 *  Created on: 2018/04/29
 *      Author: guo
 */
#include "fittings.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

double fitLambda(double T) {
	double F = 1. / (-3.09e-2 - 3.21e-7 * pow(T, 1.133));
	double lambda = pow(10, F) * exp(-8000. / T)
			+ 2e-26
					* (1e7 * sqrt(500. / pow(T, 1.4))
							* exp(-114800. / (T + 1000.))
							+ 9.e-3 * sqrt(12. * T) * exp(-92. / T));
	return lambda;
}
double fitGamma(double T) {
	return 1.e-26;
}
double fitTEMP(double numdens) {
	double x = log10(numdens);
	double fit_T_NUMDENS_a = 3.626416807912980E+00;
	double fit_T_NUMDENS_b = -1.224933568952209E+00;
	double fit_T_NUMDENS_c = -1.274578242438200E+00;
	double fit_T_NUMDENS_d = 5.219886021076703E-01;
	double fit_T_NUMDENS_e = 1.053208709509926E+00;
	double fit_T_NUMDENS_f = -3.879459648821958E-01;
	double fit_T_NUMDENS_g = -4.660615655301823E-01;
	double fit_T_NUMDENS_h = 2.356638469442452E-01;
	double fit_T_NUMDENS_i = 8.724038787337413E-02;
	double fit_T_NUMDENS_j = -7.718373457853663E-02;
	double fit_T_NUMDENS_k = 4.013294530466450E-03;
	double fit_T_NUMDENS_l = 1.047836002918579E-02;
	double fit_T_NUMDENS_m = -3.511934998554032E-03;
	double fit_T_NUMDENS_n = 3.528079479227492E-05;
	double fit_T_NUMDENS_o = 2.594946122742665E-04;
	double fit_T_NUMDENS_p = -8.283047406203526E-05;
	double fit_T_NUMDENS_q = 1.383358597907899E-05;
	double fit_T_NUMDENS_r = -1.423434735457324E-06;
	double fit_T_NUMDENS_s = 9.102911931652432E-08;
	double fit_T_NUMDENS_t = -3.340355558070870E-09;
	double fit_T_NUMDENS_u = 5.400936208142200E-11;
	double T_DENS = fit_T_NUMDENS_a + +fit_T_NUMDENS_b * pow(x, 1)
			+ fit_T_NUMDENS_c * pow(x, 2) + fit_T_NUMDENS_d * pow(x, 3)
			+ fit_T_NUMDENS_e * pow(x, 4) + fit_T_NUMDENS_f * pow(x, 5)
			+ fit_T_NUMDENS_g * pow(x, 6) + fit_T_NUMDENS_h * pow(x, 7)
			+ fit_T_NUMDENS_i * pow(x, 8) + fit_T_NUMDENS_j * pow(x, 9)
			+ fit_T_NUMDENS_k * pow(x, 10) + fit_T_NUMDENS_l * pow(x, 11)
			+ fit_T_NUMDENS_m * pow(x, 12) + fit_T_NUMDENS_n * pow(x, 13)
			+ fit_T_NUMDENS_o * pow(x, 14) + fit_T_NUMDENS_p * pow(x, 15)
			+ fit_T_NUMDENS_q * pow(x, 16) + fit_T_NUMDENS_r * pow(x, 17)
			+ fit_T_NUMDENS_s * pow(x, 18) + fit_T_NUMDENS_t * pow(x, 19)
			+ fit_T_NUMDENS_u * pow(x, 20);
	return pow(10.0, T_DENS);
}

double fitAB_E(double numdens) {
	double x = log10(numdens);
	double fit_AB_E_NUMDENS_a = -1.872011804950556E+00;
	double fit_AB_E_NUMDENS_b = -8.937317358549671E-01;
	double fit_AB_E_NUMDENS_c = -4.049467520228735E-01;
	double fit_AB_E_NUMDENS_d = 2.340570540550370E-01;
	double fit_AB_E_NUMDENS_e = 3.399024284192451E-01;
	double fit_AB_E_NUMDENS_f = -1.688149498189914E-01;
	double fit_AB_E_NUMDENS_g = -1.423033769865771E-01;
	double fit_AB_E_NUMDENS_h = 9.325026512578551E-02;
	double fit_AB_E_NUMDENS_i = 2.200234025037475E-02;
	double fit_AB_E_NUMDENS_j = -2.777434952732058E-02;
	double fit_AB_E_NUMDENS_k = 3.150027423127273E-03;
	double fit_AB_E_NUMDENS_l = 3.351038047856944E-03;
	double fit_AB_E_NUMDENS_m = -1.344002678105646E-03;
	double fit_AB_E_NUMDENS_n = 6.892793938024374E-05;
	double fit_AB_E_NUMDENS_o = 8.392156974890075E-05;
	double fit_AB_E_NUMDENS_p = -2.999536252357328E-05;
	double fit_AB_E_NUMDENS_q = 5.300810713536076E-06;
	double fit_AB_E_NUMDENS_r = -5.689318414940446E-07;
	double fit_AB_E_NUMDENS_s = 3.771620461037251E-08;
	double fit_AB_E_NUMDENS_t = -1.429771105375033E-09;
	double fit_AB_E_NUMDENS_u = 2.382784941187684E-11;
	double AB_E_DENS = fit_AB_E_NUMDENS_a + +fit_AB_E_NUMDENS_b * pow(x, 1)
			+ fit_AB_E_NUMDENS_c * pow(x, 2) + fit_AB_E_NUMDENS_d * pow(x, 3)
			+ fit_AB_E_NUMDENS_e * pow(x, 4) + fit_AB_E_NUMDENS_f * pow(x, 5)
			+ fit_AB_E_NUMDENS_g * pow(x, 6) + fit_AB_E_NUMDENS_h * pow(x, 7)
			+ fit_AB_E_NUMDENS_i * pow(x, 8) + fit_AB_E_NUMDENS_j * pow(x, 9)
			+ fit_AB_E_NUMDENS_k * pow(x, 10) + fit_AB_E_NUMDENS_l * pow(x, 11)
			+ fit_AB_E_NUMDENS_m * pow(x, 12) + fit_AB_E_NUMDENS_n * pow(x, 13)
			+ fit_AB_E_NUMDENS_o * pow(x, 14) + fit_AB_E_NUMDENS_p * pow(x, 15)
			+ fit_AB_E_NUMDENS_q * pow(x, 16) + fit_AB_E_NUMDENS_r * pow(x, 17)
			+ fit_AB_E_NUMDENS_s * pow(x, 18) + fit_AB_E_NUMDENS_t * pow(x, 19)
			+ fit_AB_E_NUMDENS_u * pow(x, 20);
	return pow(10.0, AB_E_DENS);
}
double fitAB_H2(double numdens) {
	double x = log10(numdens);
	double fit_AB_E_NUMDENS_a = -7.360098032317123E+00;
	double fit_AB_E_NUMDENS_b = -3.423714872708186E-01;
	double fit_AB_E_NUMDENS_c = -2.137800396209498E-01;
	double fit_AB_E_NUMDENS_d = 1.695543107865461E+00;
	double fit_AB_E_NUMDENS_e = 6.265122212733395E-01;
	double fit_AB_E_NUMDENS_f = -1.382649164788568E+00;
	double fit_AB_E_NUMDENS_g = -3.758096546109677E-01;
	double fit_AB_E_NUMDENS_h = 6.720367129380326E-01;
	double fit_AB_E_NUMDENS_i = 3.241466672648317E-02;
	double fit_AB_E_NUMDENS_j = -1.862316720489950E-01;
	double fit_AB_E_NUMDENS_k = 3.500919192325033E-02;
	double fit_AB_E_NUMDENS_l = 2.190695186294473E-02;
	double fit_AB_E_NUMDENS_m = -1.046778351894664E-02;
	double fit_AB_E_NUMDENS_n = 7.087659697485427E-04;
	double fit_AB_E_NUMDENS_o = 6.495692527397854E-04;
	double fit_AB_E_NUMDENS_p = -2.446944640902132E-04;
	double fit_AB_E_NUMDENS_q = 4.416837349827385E-05;
	double fit_AB_E_NUMDENS_r = -4.790959830696257E-06;
	double fit_AB_E_NUMDENS_s = 3.190719716023346E-07;
	double fit_AB_E_NUMDENS_t = -1.210306966686387E-08;
	double fit_AB_E_NUMDENS_u = 2.012517239806765E-10;
	double AB_H2_DENS = fit_AB_E_NUMDENS_a + +fit_AB_E_NUMDENS_b * pow(x, 1)
			+ fit_AB_E_NUMDENS_c * pow(x, 2) + fit_AB_E_NUMDENS_d * pow(x, 3)
			+ fit_AB_E_NUMDENS_e * pow(x, 4) + fit_AB_E_NUMDENS_f * pow(x, 5)
			+ fit_AB_E_NUMDENS_g * pow(x, 6) + fit_AB_E_NUMDENS_h * pow(x, 7)
			+ fit_AB_E_NUMDENS_i * pow(x, 8) + fit_AB_E_NUMDENS_j * pow(x, 9)
			+ fit_AB_E_NUMDENS_k * pow(x, 10) + fit_AB_E_NUMDENS_l * pow(x, 11)
			+ fit_AB_E_NUMDENS_m * pow(x, 12) + fit_AB_E_NUMDENS_n * pow(x, 13)
			+ fit_AB_E_NUMDENS_o * pow(x, 14) + fit_AB_E_NUMDENS_p * pow(x, 15)
			+ fit_AB_E_NUMDENS_q * pow(x, 16) + fit_AB_E_NUMDENS_r * pow(x, 17)
			+ fit_AB_E_NUMDENS_s * pow(x, 18) + fit_AB_E_NUMDENS_t * pow(x, 19)
			+ fit_AB_E_NUMDENS_u * pow(x, 20);

	return pow(10.0, AB_H2_DENS);
}

double fitAB_CO(double numdens) {
	double x = log10(numdens);
	double fit_AB_E_NUMDENS_a = -2.051274717031706E+01;
	double fit_AB_E_NUMDENS_b = 1.656155893322455E+00;
	double fit_AB_E_NUMDENS_c = -2.263458598301474E-01;
	double fit_AB_E_NUMDENS_d = 1.704379890823728E+00;
	double fit_AB_E_NUMDENS_e = 6.540640533501590E-01;
	double fit_AB_E_NUMDENS_f = -1.401277757381409E+00;
	double fit_AB_E_NUMDENS_g = -3.956320154899900E-01;
	double fit_AB_E_NUMDENS_h = 6.889568533047997E-01;
	double fit_AB_E_NUMDENS_i = 3.647545368628379E-02;
	double fit_AB_E_NUMDENS_j = -1.930445428728992E-01;
	double fit_AB_E_NUMDENS_k = 3.602797485651373E-02;
	double fit_AB_E_NUMDENS_l = 2.291844991727394E-02;
	double fit_AB_E_NUMDENS_m = -1.093717655323075E-02;
	double fit_AB_E_NUMDENS_n = 7.384884969774629E-04;
	double fit_AB_E_NUMDENS_o = 6.831606115135190E-04;
	double fit_AB_E_NUMDENS_p = -2.577595889622206E-04;
	double fit_AB_E_NUMDENS_q = 4.662860736266569E-05;
	double fit_AB_E_NUMDENS_r = -5.069249559221394E-06;
	double fit_AB_E_NUMDENS_s = 3.383564575974674E-07;
	double fit_AB_E_NUMDENS_t = -1.286218414053287E-08;
	double fit_AB_E_NUMDENS_u = 2.143159560810884E-10;
	double AB_H2_DENS = fit_AB_E_NUMDENS_a + +fit_AB_E_NUMDENS_b * pow(x, 1)
			+ fit_AB_E_NUMDENS_c * pow(x, 2) + fit_AB_E_NUMDENS_d * pow(x, 3)
			+ fit_AB_E_NUMDENS_e * pow(x, 4) + fit_AB_E_NUMDENS_f * pow(x, 5)
			+ fit_AB_E_NUMDENS_g * pow(x, 6) + fit_AB_E_NUMDENS_h * pow(x, 7)
			+ fit_AB_E_NUMDENS_i * pow(x, 8) + fit_AB_E_NUMDENS_j * pow(x, 9)
			+ fit_AB_E_NUMDENS_k * pow(x, 10) + fit_AB_E_NUMDENS_l * pow(x, 11)
			+ fit_AB_E_NUMDENS_m * pow(x, 12) + fit_AB_E_NUMDENS_n * pow(x, 13)
			+ fit_AB_E_NUMDENS_o * pow(x, 14) + fit_AB_E_NUMDENS_p * pow(x, 15)
			+ fit_AB_E_NUMDENS_q * pow(x, 16) + fit_AB_E_NUMDENS_r * pow(x, 17)
			+ fit_AB_E_NUMDENS_s * pow(x, 18) + fit_AB_E_NUMDENS_t * pow(x, 19)
			+ fit_AB_E_NUMDENS_u * pow(x, 20);

	return pow(10.0, AB_H2_DENS);
}

double fitAB_HI(double numdens) {
	double x = log10(numdens);
	double fit_AB_E_NUMDENS_a = -6.137627239621324E-03;
	double fit_AB_E_NUMDENS_b = 1.072591293581314E-02;
	double fit_AB_E_NUMDENS_c = -1.268358087584409E-03;
	double fit_AB_E_NUMDENS_d = 8.298703190565435E-04;
	double fit_AB_E_NUMDENS_e = -9.380209094011755E-03;
	double fit_AB_E_NUMDENS_f = 8.149209759993337E-06;
	double fit_AB_E_NUMDENS_g = 7.699114612721374E-03;
	double fit_AB_E_NUMDENS_h = -1.073911173255449E-03;
	double fit_AB_E_NUMDENS_i = -2.938966080876961E-03;
	double fit_AB_E_NUMDENS_j = 9.529549555644621E-04;
	double fit_AB_E_NUMDENS_k = 4.251472373259463E-04;
	double fit_AB_E_NUMDENS_l = -2.620579798064097E-04;
	double fit_AB_E_NUMDENS_m = 1.330515545506811E-05;
	double fit_AB_E_NUMDENS_n = 2.122109753312525E-05;
	double fit_AB_E_NUMDENS_o = -6.482216455042173E-06;
	double fit_AB_E_NUMDENS_p = 6.283216303757292E-07;
	double fit_AB_E_NUMDENS_q = 5.133851349512566E-08;
	double fit_AB_E_NUMDENS_r = -1.964694860379738E-08;
	double fit_AB_E_NUMDENS_s = 2.148351439547345E-09;
	double fit_AB_E_NUMDENS_t = -1.114467075592349E-10;
	double fit_AB_E_NUMDENS_u = 2.331356984330223E-12;
	double AB_H2_DENS = fit_AB_E_NUMDENS_a + +fit_AB_E_NUMDENS_b * pow(x, 1)
			+ fit_AB_E_NUMDENS_c * pow(x, 2) + fit_AB_E_NUMDENS_d * pow(x, 3)
			+ fit_AB_E_NUMDENS_e * pow(x, 4) + fit_AB_E_NUMDENS_f * pow(x, 5)
			+ fit_AB_E_NUMDENS_g * pow(x, 6) + fit_AB_E_NUMDENS_h * pow(x, 7)
			+ fit_AB_E_NUMDENS_i * pow(x, 8) + fit_AB_E_NUMDENS_j * pow(x, 9)
			+ fit_AB_E_NUMDENS_k * pow(x, 10) + fit_AB_E_NUMDENS_l * pow(x, 11)
			+ fit_AB_E_NUMDENS_m * pow(x, 12) + fit_AB_E_NUMDENS_n * pow(x, 13)
			+ fit_AB_E_NUMDENS_o * pow(x, 14) + fit_AB_E_NUMDENS_p * pow(x, 15)
			+ fit_AB_E_NUMDENS_q * pow(x, 16) + fit_AB_E_NUMDENS_r * pow(x, 17)
			+ fit_AB_E_NUMDENS_s * pow(x, 18) + fit_AB_E_NUMDENS_t * pow(x, 19)
			+ fit_AB_E_NUMDENS_u * pow(x, 20);
	return pow(10.0, AB_H2_DENS);
}



