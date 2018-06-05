#include "header.h"
void evolve_abundance(RealPtcl &hydro, const PS::F64 oneDynTimeStep);
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
PS::F64 getTimeStepGlobal(PS::ParticleSystem<RealPtcl>& sph_system) {
	PS::F64 dt = 1.0e+30; //set VERY LARGE VALUE
//	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
////		if(sph_system[i].eng < 0.0){
////			dt =0.8* sph_system[i].eng/sph_system[i].eng_dot;
////		}
//		dt = std::min(dt, sph_system[i].dt);
//	}

	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
		const RealPtcl& ith = sph_system[i];
		//		minsmth = fmin(minsmth, ith.smth);
		//		maxsnds = fmax(maxsnds, ith.snds);
		//		maxsgrav = fmax(maxsgrav, sqrt(ith.grav * ith.grav));
		//		maxbhgrav = fmax(maxbhgrav, sqrt(ith.extF * ith.extF));
		sph_system[i].bh_timescale = PARAM::ST
				* sqrt(ith.smth / sqrt(ith.extF * ith.extF)) / PARAM::yr;
		sph_system[i].sg_timescale = PARAM::ST
				* sqrt(ith.smth / sqrt(ith.grav * ith.grav)) / PARAM::yr;
		sph_system[i].hydro_timescale = PARAM::ST * ith.smth / (ith.snds)
				/ PARAM::yr;

		dt = fmin(0.05 * ith.smth / (ith.snds),
				0.02 * sqrt(ith.smth / sqrt(ith.grav * ith.grav)));
		dt = fmin(dt, 0.02 * sqrt(ith.smth / sqrt(ith.extF * ith.extF)));


	}
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
		const RealPtcl& ith = sph_system[i];
		sph_system[i].dt = dt;
		double etmp = sph_system[i].eng + .5 * dt * sph_system[i].eng_dot;
		//		std::cout << etmp << std::endl;
		double etmp1 = sph_system[i].eng + dt * sph_system[i].eng_dot;
		sph_system[i].eng_timescale = fabs(ith.eng / ith.eng_dot);
		if (etmp <= 0.0 || etmp1 <= 0.0) {
			sph_system[i].minus_eng = 1;
//			const PS::S32 ind[1] = { i };
//			sph_system.removeParticle(ind, 1);
			std::cout << sph_system[i].eng << "gbtime" << std::endl;
			dt = fmin(dt, 0.8 * fabs(ith.eng / ith.eng_dot));
		}
	}
	return PS::Comm::getMinValue(dt);
}

void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
		sph_system[i].vel_half = sph_system[i].vel
				+ 0.5 * dt
						* (sph_system[i].acc + sph_system[i].grav
								+ sph_system[i].extF);
		sph_system[i].eng_half = sph_system[i].eng
				+ 0.5 * dt * sph_system[i].eng_dot;
//		sph_system[i].BoverDens_half = sph_system[i].BoverDens + 0.5 * dt * sph_system[i].BoverDens_dot;
//		sph_system[i].MagneticB = sph_system[i].BoverDens_half * sph_system[i].dens;
//		sph_system[i].MagneticB.z = 0.0;

	}
}

void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt) {
	//time becomes t + dt;
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
		sph_system[i].pos += dt * sph_system[i].vel_half;

	}
}

void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
		sph_system[i].vel += dt
				* (sph_system[i].acc + sph_system[i].grav + sph_system[i].extF);
		sph_system[i].eng += dt * sph_system[i].eng_dot;
		if (sph_system[i].eng < 0.0) {
			sph_system[i].eng = 1e-4;
		}
//		sph_system[i].BoverDens += dt * sph_system[i].BoverDens_dot;
//		sph_system[i].MagneticB = sph_system[i].BoverDens * sph_system[i].dens;
//		sph_system[i].MagneticB.z = 0.0;

	}
}

void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt,
		const PS::F64 time) {
//	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//
//		sph_system[i].vel = sph_system[i].vel_half + 0.5 * dt * (sph_system[i].acc + sph_system[i].grav + sph_system[i].extF);
//		sph_system[i].eng = sph_system[i].eng_half + 0.5 * dt * sph_system[i].eng_dot;
//		if (sph_system[i].eng < 0.0) {
//			sph_system[i].eng = 1e-4;
//		}
//
////		sph_system[i].BoverDens = sph_system[i].BoverDens_half + 0.5 * dt * sph_system[i].BoverDens_dot;
////		sph_system[i].MagneticB = sph_system[i].BoverDens * sph_system[i].dens;
////		sph_system[i].MagneticB.z = 0.0;
//	}
#pragma omp parallel for
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		if (sph_system[i].dens * PARAM::SMassDens / (2.34 * PARAM::PROTONMASS)
				> 1.e10 || sqrt(sph_system[i].grav * sph_system[i].grav) > 50.0
				|| ((sph_system[i].pos.x < -20 || sph_system[i].pos.x > 20)
						&& time > 1.5)
				|| sph_system[i].cooling_timescale > 1.2e6 && time > 1.3) {
			const PS::S32 ind[1] = { i };
			sph_system.removeParticle(ind, 1);
		}
		sph_system[i].vel_half = sph_system[i].vel
				+ 0.5 * dt
						* (sph_system[i].acc + sph_system[i].grav
								+ sph_system[i].extF);
//
		sph_system[i].vel = sph_system[i].vel_half
				+ 0.5 * dt
						* (sph_system[i].acc + sph_system[i].grav
								+ sph_system[i].extF);
//		sph_system[i].eng = sph_system[i].eng + dt * sph_system[i].eng_dot;
//		sph_system[i].pos += dt * sph_system[i].vel;
		sph_system[i].pos += sph_system[i].vel * dt;
		sph_system[i].vel += dt
				* (sph_system[i].acc + sph_system[i].grav + sph_system[i].extF);
		sph_system[i].eng += dt * sph_system[i].eng_dot;
		//TODO Beaware that eng can already <0

		if (sph_system[i].eng < 0.0) {
			std::cout << sph_system[i].eng << "in integral.cpp " << std::endl;
			double NUMDENS_CGS = sph_system[i].dens * PARAM::SMassDens
					/ (2.34 * PARAM::PROTONMASS);
			double T = fitTEMP(NUMDENS_CGS);
			double abundance_e = fitAB_E(NUMDENS_CGS);
			double abundance_H2 = fitAB_H2(NUMDENS_CGS);
			double energy = T
					* (PARAM::KBOLTZ_CGS * (1.1 + abundance_e - abundance_H2))
					/ (2.34 * PARAM::PROTONMASS_CGS * (5. / 3. - 1.0));
			sph_system[i].eng = energy / PARAM::SEng_per_Mass;

		} else {
			evolve_abundance(sph_system[i], dt);
			if (sph_system[i].eng < 0.0) {
				std::cout << sph_system[i].eng << "in integral.cccccc "
						<< std::endl;

			}
		}

	}

//	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//			sph_system[i].pos += (sph_system[i].acc + sph_system[i].grav) * sph_system[i].dt * sph_system[i].dt * .5;
//	//		if(sqrt(sph_system[i].pos*sph_system[i].pos) >= PARAM::xi ){
//	//			sph_system[i].pos = 0.001;
//	//		}
//	//		if(sph_system[i].pos.x != sph_system[i].pos.x){
//	//			std::cout<<sph_system[i].acc<<std::endl;
//	//		}
//		}
}

PS::F64 RATE_XR_HII_HeII(PS::F64 T, PS::F64 abundance_e, PS::F64 abundance_HI,
		PS::F64 col_dens) {

	PS::F64 p1 = log10(col_dens / 1.e18);
	//	if (abundance_HI < 1e-10) {
	//		p1 = log10(1e-10 * col_dens / 1.e18);
	//	}
	PS::F64 p2 = log10(abundance_e);
	PS::F64 f_4 = 1.06 + 4.08e-2 * p2 + 6.51e-3 * pow(p2, 2.);
	PS::F64 f_5 = 1.90 + 0.678 * p2 + 0.113 * p2 * p2;
	PS::F64 rate = pow(10.,
			f_4 * (-15.6 - 1.1 * p1 + 9.13e-2 * p1 * p1)
					+ f_5 * 0.87 * exp(-pow((p1 - 0.41) / 0.84, 2.)));
	return rate;
}
// s-1///
PS::F64 RATE_CR_HII_HeII(PS::F64 abundance_e) {
	PS::F64 E = 35.;
	PS::F64 phi_H = (E / 13.6 - 1) * .3908
			* pow(1. - pow(abundance_e, 0.4092), 1.7592)
			* pow(E / 1.e3, abundance_e * 2.313);
	PS::F64 phi_He = (E / 24.6 - 1.) * 0.0554
			* pow(1. - pow(abundance_e, 0.4614), 1.6660)
			* pow(E / 1.e3, abundance_e * 2.313);
	PS::F64 xi = PARAM::CR_ION_RATE * (1. + phi_H + phi_He);
	return xi;
}
// cm3 s-1///
PS::F64 RATE_HH2_3H(PS::F64 T) {

	PS::F64 ch1 = 0.;
	ch1 = (T > 2000.) ? (3.4e-9 * exp(-43900. / T)) : 0;
	return ch1;
}
// cm3 s-1///
PS::F64 RATE_HHm_H2_e(PS::F64 T, PS::F64 abundance_HI, PS::F64 abundance_HII) {
	PS::F64 k8 = 0.;
	k8 =
			8.5e-16 * pow(T / 1000, 0.8)
					* (abundance_HI
							/ (abundance_HI
									+ 53.0 * pow(T / 1000, -.4) * abundance_HII));
	return k8;
}
// cm3 s-1///
PS::F64 RATE_H_col_HII(PS::F64 T) {
	return 5.8e-9 * pow(T / 1e4, .5) * exp(-15.8 / (T / 1e4));
}
// s-1///
PS::F64 RATE_H2_UV_2H(PS::F64 T, PS::F64 col_den_H2, PS::F64 mu,
		PS::F64 abundance_H2) {

	PS::F64 snds = sqrt(PARAM::KBOLTZ_cgs * T / (PARAM::PROTONMASS)) / 1e5;
	PS::F64 NH2 = abundance_H2 * col_den_H2;
	PS::F64 tau = 1.2e-14 * NH2 / snds;
	PS::F64 R_pump = 0.0;
	if (tau > 10) {
		PS::F64 b = 9.2e-3 / snds;
		PS::F64 v1 = 5.0e2 / snds;
		PS::F64 beta_tau = erfc(sqrt((pow(v1, -2) * tau * b / M_PI)))
				* (pow(log(tau / sqrt(M_PI)), -.5) / tau + sqrt(b / tau));
		R_pump = 3.4e-10 * PARAM::G0 * beta_tau;
	} else {
		int N = 25;
		PS::F64 beta_tau = 0;
		for (int n = 0; n < N; n++) {
			long long factorial = 1;
			for (int i = 1; i <= n; i++) {

				factorial = factorial * i;
			}
			beta_tau += pow(-tau, n)
					/ (factorial * sqrt(n + 1) * pow(M_PI, .5 * n));
		}
		R_pump = 3.4e-10 * PARAM::G0 * beta_tau;
	}
	return R_pump;
}
// s-1///
PS::F64 RATE_H2_CR_2H(PS::F64 abundance_e) {
	return 2.29 * RATE_CR_HII_HeII(abundance_e);
}
//cm3 s-1///
PS::F64 RATE_H2_e_2H(PS::F64 T) {
	PS::F64 qT =
			(T > 3000.) ?
					(241600.0 * pow(T, -0.912) * exp(-55800 / T) * 1e-9) : 0;
	return qT;

}
// cm3 s-1///
double RATE_HHGrain_H2(double T, double tdust) {
	double numerator = 3.e-17 * pow(T / 100.0, .5);
	double denominator = 1.0 + 0.04 * pow(T + tdust, .5) + 2e-3 * T
			+ 8e-6 * T * T;
	double rate = numerator / denominator;
	return rate;
}
// cm3 s-1///
double RATE_HIIe_HPhoton(double T) {
	double phi_T = 0;
	if (T <= 15780) {
		phi_T = .45 * log(1.578e5 / T);
	} else {
		phi_T = .4 * pow(1.578e5 / T, .5);
	}
	double ch13 = 2.06e-11 * pow(T, -.5) * phi_T;

	return ch13;
}

// erg s-1///
double Gamma_UV(double T, double NUMDENS_CGS, double abundance_H2,
		double abundance_H, double col_den_H2, double mu) {

	double n_cr_H2 = 1e6 * pow(T, -.5)
			/ (1.6 * abundance_H * exp(-pow(400.0 / T, 2))
					+ 1.4 * 2 * abundance_H2 * exp(-(12000. / (T + 1200.))));
	double R_pump = RATE_H2_UV_2H(T, col_den_H2, mu, abundance_H2);
	double gamma = (9.0 * R_pump * (2.2 * PARAM::EV_CGS)
			/ (1.0 + n_cr_H2 / (NUMDENS_CGS))) * abundance_H2;
	if (gamma < 1e-32) {
		gamma = 0.0;
	}
	return gamma;
}
// erg s-1///
double Gamma_H2_associative_forming(double T, double NUMDENS_CGS,
		double abundance_H2, double abundance_H, double abundance_Hp,
		double abundance_e) {
	double x2 = 2 * abundance_H2;
	double rates_chem = 0.0;
	double Rate = RATE_HHm_H2_e(T, abundance_H, abundance_Hp);
	double n_cr_H2 = 1e6 * pow(T, -.5)
			/ (1.6 * abundance_H * exp(-pow(400.0 / T, 2))
					+ 1.4 * x2 * exp(-12.0 / (T / 1000.0 + 1.2)));
	rates_chem = Rate * abundance_H * abundance_e * NUMDENS_CGS * NUMDENS_CGS
			* (3.53 / (1.0 + n_cr_H2 / NUMDENS_CGS)) * PARAM::EV_CGS;

	if (rates_chem < 1e-32) {
		rates_chem = 0.0;
	}
	return rates_chem;
}
// erg s-1///
double Gamma_H2_grain_forming(double T, double tdust, double NUMDENS_CGS,
		double abundance_H2, double abundance_HI, double col_den_H2) {
	double x2 = 2 * abundance_H2;
	double rates_chem = 0.0;
	double Rate = RATE_HHGrain_H2(T, tdust) * NUMDENS_CGS * abundance_HI
			* NUMDENS_CGS;
	double n_cr_H2 = 1.e6 * pow(T, -.5)
			/ (1.6 * abundance_HI * exp(-pow(400.0 / T, 2))
					+ 1.4 * x2 * exp(-12000. / (T + 1200.)));
	rates_chem = Rate * (0.2 + 4.2 / (1.0 + n_cr_H2 / NUMDENS_CGS))
			* PARAM::EV_CGS * abundance_HI;

	if (rates_chem < 1e-32) {
		rates_chem = 0.0;
	}
	return rates_chem;
}
// erg  s-1///
double Lambda_Mol(double T, double abundance_e, double NUMDENS_CGS) {
	if (NUMDENS_CGS > 100) {

		double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	}
	double alpha = 5.e-30;
	double beta = 0.6 + 0.41 * log10(NUMDENS_CGS);
	double delta = 2.3 - 0.18 * log10(NUMDENS_CGS);
	double lyman_alpha_cooling = alpha * pow(T / 10., beta)
			* pow(NUMDENS_CGS, delta);
	return lyman_alpha_cooling;
}
// erg  s-1///
double Lambda_Lya(double T, double abundance_e, double NUMDENS_CGS) {
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double lyman_alpha_cooling = NUMDENS_CGS_e * 7.5e-19 * exp(-118400.0 / T);
	return lyman_alpha_cooling;
}

double Lambda_Col_Dust(double T, double abundance_H, double abundance_e,
		double NUMDENS_CGS, double grain_size, double T_gr) {
	double rates_chem5 = 1.2e-31 * NUMDENS_CGS * sqrt(T / 1000)
			* sqrt(100 / grain_size) * (1 - .8 * exp(-75.0 / T)) * (T - T_gr);
	return rates_chem5 * NUMDENS_CGS;
}
// erg s-1///
double Lambda_Recomb_Grain_PAH(double T, double abundance_e,
		double NUMDENS_CGS) {
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double eta = pow(T, .94)
			* pow(PARAM::G0 * pow(T, .5) / NUMDENS_CGS_e, .735 / pow(T, .068));
	double cooling = 4.65e-30 * eta * NUMDENS_CGS_e;
	return cooling;
}
// erg s-1///
double Gamma_pe(double T, double NUMDENS_CGS, double abundance_e) {
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double var = PARAM::G0 * sqrt(T) / NUMDENS_CGS_e;
	double epsilon = 4.9e-2 / (1. + pow(var / 1925., .73))
			+ 3.7e-2 * pow(T / 1.e4, .7) / (1. + var / 5000.0);
	double Gamma_pe = 1. * 1.e-24 * epsilon * PARAM::G0;
	return Gamma_pe;

}
// erg s-1///
double Gamma_CR(double NUMDENS_CGS, double abundance_e) {
	double Xi_CR = 1.8e-17;

	double Gamma_CR = 0.;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double E = 35.;

	double f_1 = .9971 * (1. - pow(1. - pow(abundance_e, 0.2663), 1.3163));
	double f_2 = 2 * E - 13.6;
	double f_3 = 51.0375 + 22.7858 * (log10(abundance_e))
			+ 2.7152 * pow(log10(abundance_e), 2);
	double E_h = E * (f_1 + pow(f_2 / f_3 + 1. / (1. - f_1), -1.));

	Gamma_CR = RATE_CR_HII_HeII(abundance_e) * E_h * PARAM::EV_CGS;
	return Gamma_CR;

}
// erg s-1///
double Gamma_XR(double abundance_e, double col_dens) {

	double p1 = log10(col_dens / 1e18);
	double p2 = log10(abundance_e);
	double f6 = .990 - 2.74e-3 * p2 + 1.13e-3 * p2 * p2;
	double Gamma_X_Ray = pow(10.0,
			f6 * (-26.5 - 0.920 * p1 + 5.89e-2 * p1 * p1)
					+ f6 * 0.96 * exp(-pow((p1 - 0.38) / 0.87, 2.0)));
	return Gamma_X_Ray;
}
double Lambda_Meta_FeII(double T, double NUMDENS_CGS, double abundance_H,
		double abundance_H2, double abundance_e, double abundance_Hp,
		double abundance_FeII, double tau) {
	double NUMDENS_CGS_H = abundance_H * NUMDENS_CGS;
	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double NUMDENS_CGS_Hp = abundance_Hp * NUMDENS_CGS;
	double NUMDENS_CGS_FeII = abundance_FeII * NUMDENS_CGS;
	//Aij
	double A21 = 2.8e-4;
	double A31 = 1.0e-4;
	double A32 = 6.1e-3;
	//Eij
	double E21 = 2.7e3 * PARAM::KBOLTZ_CGS;
	double E32 = 8.1e2 * PARAM::KBOLTZ_CGS;
	double E31 = E21 + E32;
	//Hydrogen 1-0
	double C21_H = 1e-12;
	double C31_H = 1e-12;

	double C32_H = 1e-12;

	double C21_e = 0.0;
	double C31_e = 0.0;
	double C32_e = 0.0;
	//electron
	C21_e = T < 1.e4 ?
			2.2e-8 * pow(T / 10000, -.5) : 2.2e-8 * pow(T / 10000, -.5);
	C31_e = T < 1.e4 ?
			7.1e-9 * pow(T / 10000, -.5) : 7.1e-9 * pow(T / 10000, -.5);
	C32_e = T < 1.e4 ?
			3.8e-8 * pow(T / 10000, -.5) : 3.8e-8 * pow(T / 10000, -.5);

	double C21_H2 = 0.0;
	double C21_Hp = 0.0;
	double C31_H2 = 0.0;
	double C31_Hp = 0.0;
	double C32_H2 = 0.0;
	double C32_Hp = 0.0;
	double C21 = C21_H * NUMDENS_CGS_H + C21_H2 * NUMDENS_CGS_H2
			+ C21_e * NUMDENS_CGS_e + C21_Hp * NUMDENS_CGS_Hp;
	double C31 = C31_H * NUMDENS_CGS_H + C31_H2 * NUMDENS_CGS_H2
			+ C31_e * NUMDENS_CGS_e + C31_Hp * NUMDENS_CGS_Hp;
	double C32 = C32_H * NUMDENS_CGS_H + C32_H2 * NUMDENS_CGS_H2
			+ C32_e * NUMDENS_CGS_e + C32_Hp * NUMDENS_CGS_Hp;
	double cooling = 0.;

	double g1 = 10, g2 = 10, g3 = 8;

	double C12 = (g2 / g1) * C21 * exp(-E21 / (PARAM::KBOLTZ_CGS * T));
	double C13 = (g3 / g1) * C31 * exp(-E31 / (PARAM::KBOLTZ_CGS * T));
	double C23 = (g3 / g2) * C32 * exp(-E32 / (PARAM::KBOLTZ_CGS * T));
	double B21, B31, B32, B12, B13, B23 = 0.0;

	double x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	double _alpha, _beta, _gamma, _delta;
	B21 = A21 / (exp(x) - 1.);
	B12 = (g2 / g1) * B21;
	x = E31 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B31 = A31 / (exp(x) - 1.);
	B13 = (g3 / g1) * B31;
	x = E32 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B32 = A32 / (exp(x) - 1.);
	B23 = (g3 / g2) * B32;

	double Beta_esc = 1.0, N1, N2, N3;

//	Beta_esc =
//			tau < 7.0 ?
//					(1.0 - exp(-2.34 * tau) / (4.68 * tau)) :
//					1.0 / (4.0 * tau * sqrt(log(tau / sqrt(M_PI))));
	_alpha = A21 * Beta_esc + B21 * B21 + C21;
	_beta = A31 * Beta_esc + B31 * Beta_esc + C31 + A32 * Beta_esc
			+ B32 * Beta_esc + C32;
	_gamma = B21 * Beta_esc + C12;
	_delta = (A21 * Beta_esc + B21 * Beta_esc + C21) / (B21 * Beta_esc + C12)
			+ B23 * Beta_esc + C23;
	N2 = _gamma * _alpha * _beta
			/ (_alpha * _alpha * _beta + _alpha * _beta
					+ _delta * _gamma * _gamma);
	N1 = N2 * _alpha / _gamma;
	N3 = N2 * _delta * _gamma / (_alpha * _beta);
	double var1 = (A21 + B21) * E21 * N2;
	double var2 = B12 * E21 * N1;
	double var3 = (A31 + B31) * E31 * N3;
	double var4 = B13 * E31 * N1;
	double var5 = (A32 + B32) * E32 * N3;
	double var6 = B23 * E32 * N2;
	cooling = ((var1 - var2) * (g2 / g1) * exp(-E21 / T)
			+ (var3 - var4) * (g3 / g1) * exp(-E31 / T)
			+ (var5 - var6) * (g3 / g2) * exp(-E32 / T)) * abundance_FeII
			* NUMDENS_CGS;
//	cooling = NUMDENS_CGS * abundance_e * NUMDENS_CGS * abundance_FeII
//			* (exp(-2694 / T)) * 4.8e-18 / sqrt(T)
//			+ NUMDENS_CGS * abundance_e * NUMDENS_CGS * abundance_FeII
//					* (exp(-3496 / T)) * 7.8e-18 / sqrt(T);
	return cooling;
}
double Lambda_FeII(double T, double NUMDENS_CGS, double abundance_HI,
		double abundance_H2, double abundance_e, double abundance_HII,
		double abundance_FeII, double tau) {
	double NUMDENS_CGS_H = abundance_HI * NUMDENS_CGS;
	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double NUMDENS_CGS_Hp = abundance_HII * NUMDENS_CGS;
	double NUMDENS_CGS_FeII = abundance_FeII * NUMDENS_CGS;
	//Aij
	double A21 = 2.1e-3;
	double A31 = 1.5e-9;
	double A32 = 1.6e-3;
	//Eij
	double E21 = 5.6e2 * PARAM::KBOLTZ_CGS;
	double E32 = 4.1e2 * PARAM::KBOLTZ_CGS;
	double E31 = E21 + E32;
	//Hydrogen 1-0
	double C21_H = 9.5e-10;
	//electron
	double C21_e = 1.8e-6 * pow(T / 100, -.5);
	//2-0
	double C31_H = 5.7e-10;
	double C31_e = 1.8e-6 * pow(T / 100, -.5);
	//2-1
	double C32_H = 4.7e-10;
	double C32_e = 8.7e-7 * pow(T / 100, -.5);

	double C21_H2 = 0.0;
	double C21_Hp = 0.0;
	double C31_H2 = 0.0;
	double C31_Hp = 0.0;
	double C32_H2 = 0.0;
	double C32_Hp = 0.0;
	double C21 = C21_H * NUMDENS_CGS_H + C21_H2 * NUMDENS_CGS_H2
			+ C21_e * NUMDENS_CGS_e + C21_Hp * NUMDENS_CGS_Hp;
	double C31 = C31_H * NUMDENS_CGS_H + C31_H2 * NUMDENS_CGS_H2
			+ C31_e * NUMDENS_CGS_e + C31_Hp * NUMDENS_CGS_Hp;
	double C32 = C32_H * NUMDENS_CGS_H + C32_H2 * NUMDENS_CGS_H2
			+ C32_e * NUMDENS_CGS_e + C32_Hp * NUMDENS_CGS_Hp;
	double cooling = 0.;

	double g1 = 10, g2 = 8, g3 = 6;

	double C12 = (g2 / g1) * C21 * exp(-E21 / (PARAM::KBOLTZ_CGS * T));
	double C13 = (g3 / g1) * C31 * exp(-E31 / (PARAM::KBOLTZ_CGS * T));
	double C23 = (g3 / g2) * C32 * exp(-E32 / (PARAM::KBOLTZ_CGS * T));
	double B21, B31, B32, B12, B13, B23 = 0.0;

	double x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	double _alpha, _beta, _gamma, _delta;
	B21 = A21 / (exp(x) - 1.);
	B12 = (g2 / g1) * B21;
	x = E31 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B31 = A31 / (exp(x) - 1.);
	B13 = (g3 / g1) * B31;
	x = E32 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B32 = A32 / (exp(x) - 1.);
	B23 = (g3 / g2) * B32;

	double Beta_esc = 1.0, N1, N2, N3, _theta;
	//	tau = z;
	//	if (tau < 7.0) {
	//		Beta_esc = (1.0 - exp(-2.34 * tau)) / (4.68 * tau);
	//	} else {
	//		Beta_esc = 1.0 / (4.0 * tau * sqrt(log(tau / sqrt(M_PI))));
	//	}

	_alpha = A21 * Beta_esc + B21 * Beta_esc + C21;
	_beta = C12 + B12 * Beta_esc;
	_gamma = B13 * Beta_esc + C13;
	_theta = A31 * Beta_esc + B31 * Beta_esc + C31 + A32 * Beta_esc
			+ B32 * Beta_esc + C32;
	_delta = B23 * Beta_esc + C23;
	N2 = 1.0
			/ (_alpha / _beta + 1.0 + _delta / _theta
					+ _gamma * _alpha / (_beta * _theta));
	N1 = N2 * _alpha / _beta;
	N3 = N2 * (_delta / _theta + _gamma * _alpha / (_beta * _theta));
	double var1 = (A21 + B21) * E21 * N2;
	double var2 = B12 * E21 * N1;
	double var3 = (A31 + B31) * E31 * N3;
	double var4 = B13 * E31 * N1;
	double var5 = (A32 + B32) * E32 * N3;
	double var6 = B23 * E32 * N2;
	cooling = ((var1 - var2) * (g2 / g1) * exp(-E21 / T)
			+ (var3 - var4) * (g3 / g1) * exp(-E31 / T)
			+ (var5 - var6) * (g3 / g2) * exp(-E32 / T)) * abundance_FeII
			* NUMDENS_CGS;
//	cooling = NUMDENS_CGS * abundance_e * NUMDENS_CGS * abundance_FeII
//			* (exp(-554./ T) + 1.3 * exp(-961. / T)) * 1.1e-18 / sqrt(T)
//			+ NUMDENS_CGS * abundance_HI * NUMDENS_CGS * abundance_FeII
//					* (exp(-554. / T) + 1.3 * exp(-961. / T)) * 1.1e-22;
	return cooling;
}
double Lambda_Meta_OI(double T, double NUMDENS_CGS, double abundance_H,
		double abundance_H2, double abundance_e, double abundance_Hp,
		double abundance_OI, double tau) {

	double NUMDENS_CGS_H = abundance_H * NUMDENS_CGS;
	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double NUMDENS_CGS_Hp = abundance_Hp * NUMDENS_CGS;
	double NUMDENS_CGS_OI = abundance_OI * NUMDENS_CGS;
	//Aij
	double A21 = 6.7e-3;
	double A31 = 6.7e-2;
	double A32 = 1.3;
	//Eij
	double E21 = 2.3e4 * PARAM::KBOLTZ_CGS;
	double E32 = 2.6e4 * PARAM::KBOLTZ_CGS;
	double E31 = E21 + E32;
	//Hydrogen 1-0
	double C21_H = 1e-12;
	double C31_H = 1e-12;

	double C32_H = 1e-12;

	double C21_e = 0.0;
	double C31_e = 0.0;
	double C32_e = 0.0;
	//electron
	C21_e = T < 1.e4 ?
			5.1e-9 * pow(T / 10000., .57) : 5.1e-9 * pow(T / 10000., .17);
	C31_e = T < 1.e4 ?
			2.5e-9 * pow(T / 10000., .57) : 2.5e-9 * pow(T / 10000., .13);
	C32_e = T < 1.e4 ?
			5.2e-9 * pow(T / 10000., .57) : 5.2e-9 * pow(T / 10000., .15);

	double C21_H2 = 0.0;
	double C21_Hp = 0.0;
	double C31_H2 = 0.0;
	double C31_Hp = 0.0;
	double C32_H2 = 0.0;
	double C32_Hp = 0.0;
	double C21 = C21_H * NUMDENS_CGS_H + C21_H2 * NUMDENS_CGS_H2
			+ C21_e * NUMDENS_CGS + C21_Hp * NUMDENS_CGS_Hp;
	double C31 = C31_H * NUMDENS_CGS_H + C31_H2 * NUMDENS_CGS_H2
			+ C31_e * NUMDENS_CGS + C31_Hp * NUMDENS_CGS_Hp;
	double C32 = C32_H * NUMDENS_CGS_H + C32_H2 * NUMDENS_CGS_H2
			+ C32_e * NUMDENS_CGS + C32_Hp * NUMDENS_CGS_Hp;
	double cooling = 0.;

	double g1 = 9., g2 = 5., g3 = 1.;

	double C12 = (g2 / g1) * C21 * exp(-E21 / (PARAM::KBOLTZ_CGS * T));
	double C13 = (g3 / g1) * C31 * exp(-E31 / (PARAM::KBOLTZ_CGS * T));
	double C23 = (g3 / g2) * C32 * exp(-E32 / (PARAM::KBOLTZ_CGS * T));
	double B21, B31, B32, B12, B13, B23 = 0.0;

	double x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	double _alpha, _beta, _gamma, _delta, _theta;

	B21 = A21 / (exp(x) - 1.);

	B12 = (g2 / g1) * B21;

	x = E31 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B31 = 0.0;
	B31 = A31 / (exp(x) - 1.);

	B13 = (g3 / g1) * B31;
	B32 = 0.0;
	x = E32 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B32 = A32 / (exp(x) - 1.);
	B23 = (g3 / g2) * B32;
	double Beta_esc = 1.0, N1, N2, N3;
//	tau = z;
//	if (tau < 7.0) {
//		Beta_esc = (1.0 - exp(-2.34 * tau)) / (4.68 * tau);
//	} else {
//		Beta_esc = 1.0 / (4.0 * tau * sqrt(log(tau / sqrt(M_PI))));
//	}

	_alpha = A21 * Beta_esc + B21 * Beta_esc + C21;
	_beta = C12 + B12 * Beta_esc;
	_gamma = B13 * Beta_esc + C13;
	_theta = A31 * Beta_esc + B31 * Beta_esc + C31 + A32 * Beta_esc
			+ B32 * Beta_esc + C32;
	_delta = B23 * Beta_esc + C23;
	N2 = 1.0
			/ (_alpha / _beta + 1.0 + _delta / _theta
					+ _gamma * _alpha / (_beta * _theta));
	N1 = N2 * _alpha / _beta;
	N3 = N2 * (_delta / _theta + _gamma * _alpha / (_beta * _theta));
	double var1 = (A21 + B21) * E21 * N2;
	double var2 = B12 * E21 * N1;
	double var3 = (A31 + B31) * E31 * N3;
	double var4 = B13 * E31 * N1;
	double var5 = (A32 + B32) * E32 * N3;
	double var6 = B23 * E32 * N2;
	cooling = ((var1 - var2) * (g2 / g1) * exp(-E21 / T)
			+ (var3 - var4) * (g3 / g1) * exp(-E31 / T)
			+ (var5 - var6) * (g3 / g2) * exp(-E32 / T)) * abundance_OI
			* NUMDENS_CGS;
//	cooling = ((var1 - var2) * (g2 / g1) * exp(-E21 / T)) * abundance_OI
//			* NUMDENS_CGS;
////
//	double n_H = NUMDENS_CGS_H + 2. * NUMDENS_CGS_H2;
//	double f_CII = 2. * NUMDENS_CGS_H2 / n_H;
//	cooling = 1.8e-24 * NUMDENS_CGS * abundance_OI
//			* (NUMDENS_CGS_H + NUMDENS_CGS_H2) * exp(-22800 / T);
//	cooling = NUMDENS_CGS * abundance_OI * NUMDENS_CGS * abundance_e
//			* (9.4e-23 * sqrt(T) * exp(-22700 / T));
//	if (cooling < 1e31) {
//		cooling = 0.0;
//	}
	return cooling;
}

void three_level_population(double r01, double r02, double r12, double r10,
		double r20, double r21, double &n0, double &n1, double &n2) {

	double a1, a2, a3, b1, b2, b3;

	a1 = r01 + r02;
	a2 = -r10;
	a3 = -r20;
	b1 = r01;
	b2 = -(r10 + r12);
	b3 = r21;

	n2 = (r01 == 0.0 && r02 == 0.0) ?
			0.0 :
			(-a1 * (a1 * b2 - b1 * a2)
					/ ((a1 - a2) * (a1 * b3 - b1 * a3)
							- (a1 - a3) * (a1 * b2 - b1 * a2)));

	n1 = (r01 == 0.0 && r02 == 0.0) ?
			0.0 : ((a1 / (a1 - a2)) - ((a1 - a3) / (a1 - a2)) * n2);

	n0 = (r01 == 0.0 && r02 == 0.0) ? 1.0 : (1.0 - n1 - n2);
}
double Lambda_OI(double T, double NUMDENS_CGS, double abundance_H,
		double abundance_H2, double abundance_e, double abundance_Hp,
		double abundance_OI, double tau, double col_dens_H2) {
//	double NUMDENS_CGS_H = abundance_H * NUMDENS_CGS;
//	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
//	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
//	double NUMDENS_CGS_Hp = abundance_Hp * NUMDENS_CGS;
//	double NUMDENS_CGS_FeII = abundance_OI * NUMDENS_CGS;
//	//Aij
//	double A21 = 9.0e-5;
//	double A31 = 1.0e-10;
//	double A32 = 1.7e-5;
//	//Eij
//	double E21 = 2.3e2 * PARAM::KBOLTZ_CGS;
//	double E32 = 9.8e1 * PARAM::KBOLTZ_CGS;
//	double E31 = E21 + E32;
//	//Hydrogen 1-0
//	double T2 = T * 0.01;
//	double C21_H = 9.2e-11 * pow(T2, 0.67);
//	double C31_H = 4.3e-11 * pow(T2, 0.80);
//
//	double C32_H = 1.1e-10 * pow(T2, 0.44);
//
//	double C21_e = 1.4e-8;
//	double C31_e = 1.4e-8;
//	double C32_e = 5.0e-9;
//	//electron
//
//	double C21_H2 = 0.0;
//	double C21_Hp = 0.0;
//	double C31_H2 = 0.0;
//	double C31_Hp = 0.0;
//	double C32_H2 = 0.0;
//	double C32_Hp = 0.0;
//	double C21 = C21_H * NUMDENS_CGS_H + C21_H2 * NUMDENS_CGS_H2
//			+ C21_e * NUMDENS_CGS_e + C21_Hp * NUMDENS_CGS_Hp;
//	double C31 = C31_H * NUMDENS_CGS_H + C31_H2 * NUMDENS_CGS_H2
//			+ C31_e * NUMDENS_CGS_e + C31_Hp * NUMDENS_CGS_Hp;
//	double C32 = C32_H * NUMDENS_CGS_H + C32_H2 * NUMDENS_CGS_H2
//			+ C32_e * NUMDENS_CGS_e + C32_Hp * NUMDENS_CGS_Hp;
//	double cooling = 0.;
//
//	double g1 = 5., g2 = 3., g3 = 1.;
//
//	double C12 = (g2 / g1) * C21 * exp(-E21 / (PARAM::KBOLTZ_CGS * T));
//	double C13 = (g3 / g1) * C31 * exp(-E31 / (PARAM::KBOLTZ_CGS * T));
//	double C23 = (g3 / g2) * C32 * exp(-E32 / (PARAM::KBOLTZ_CGS * T));
//	double B21, B31, B32, B12, B13, B23 = 0.0;
//
//	double x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
//	double _alpha, _beta, _gamma, _delta;
//	B21 = A21 / (exp(x) - 1.);
//	B12 = (g2 / g1) * B21;
//	x = E31 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
//	B31 = A31 / (exp(x) - 1.);
//	B13 = (g3 / g1) * B31;
//	x = E32 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
//	B32 = A32 / (exp(x) - 1.);
//	B23 = (g3 / g2) * B32;
//
//	double Beta_esc = 1.0, N1, N2, N3;
//
//	//	Beta_esc =
//	//			tau < 7.0 ?
//	//					(1.0 - exp(-2.34 * tau) / (4.68 * tau)) :
//	//					1.0 / (4.0 * tau * sqrt(log(tau / sqrt(M_PI))));
//	_alpha = A21 * Beta_esc + B21 * B21 + C21;
//	_beta = A31 * Beta_esc + B31 * Beta_esc + C31 + A32 * Beta_esc
//			+ B32 * Beta_esc + C32;
//	_gamma = B21 * Beta_esc + C12;
//	_delta = (A21 * Beta_esc + B21 * Beta_esc + C21) / (B21 * Beta_esc + C12)
//			+ B23 * Beta_esc + C23;
//	N2 = _gamma * _alpha * _beta
//			/ (_alpha * _alpha * _beta + _alpha * _beta
//					+ _delta * _gamma * _gamma);
//	N1 = N2 * _alpha / _gamma;
//	N3 = N2 * _delta * _gamma / (_alpha * _beta);
//	double var1 = (A21 + B21) * E21 * N2;
//	double var2 = B12 * E21 * N1;
//	double var3 = (A31 + B31) * E31 * N3;
//	double var4 = B13 * E31 * N1;
//	double var5 = (A32 + B32) * E32 * N3;
//	double var6 = B23 * E32 * N2;
//	cooling = ((var1 - var2) * (g2 / g1) * exp(-E21 / T)
//			+ (var3 - var4) * (g3 / g1) * exp(-E31 / T)
//			+ (var5 - var6) * (g3 / g2) * exp(-E32 / T)) * abundance_OI
//			* NUMDENS_CGS;
	double NUMDENS_CGS_H = abundance_H * NUMDENS_CGS;
	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double NUMDENS_CGS_Hp = abundance_Hp * NUMDENS_CGS;
	double NUMDENS_CGS_OI = abundance_OI * NUMDENS_CGS;

	//	double cooling = 0;
	//	double C21_H = 0.;
	//	double C21_e = 0;
	//	double temp2 = T * 1e-2;
	//	double f = 0;
	//	double hh = 0;
	//
	//	double C12 = 0;
	//	double PARAM::CMB_T = 2.7;
	//	double x, g1, g2;
	//	double Beta_esc = 0.0;
	//
	//	g1 = 5;
	//	g2 = 3;
	//	if (tau < 7.0) {
	//		Beta_esc = 1.0 - exp(-2.34 * tau) / (4.68 * tau);
	//	} else {
	//		Beta_esc = 1.0 / (4.0 * tau * sqrt(log(tau / sqrt(M_PI))));
	//	}
	//	Beta_esc = 1.0;
	//	C21_H = 9.2e-11 * pow(temp2, 0.67);
	//	C21_e = 1.4e-8;
	//	C21 = C21_H * NUMDENS_CGS_H + C21_e * NUMDENS_CGS_e;
	//	C12 = C21 * (g2 / g1) * exp(-E21 / (PARAM::KBOLTZ_CGS * T));
	//o
	//	x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	//	B21 = A21 / (exp(x) - 1.);
	//	B12 = (g2 / g1) * B21;
	//	LHS = (C21 + A21 * Beta_esc + B21 * Beta_esc);
	//	RHS = (C12 + B12 * Beta_esc);
	//	N1 = LHS / (LHS + RHS);
	//	N2 = RHS / (LHS + RHS);
	//	var1 = (A21 + B21) * E21 * N2;
	//	var2 = B12 * E21 * N1;
	//	cooling = (var1 - var2) * (g2 / g1) * exp(-E21 / T) * abundance_OI * NUMDENS_CGS;
	double z = PARAM::col_dens_H2 / NUMDENS_CGS;
	//Aij
	//	double A21 = 8.86e-5;
	//	double A31 = 1.275e-10;
	//	double A32 = 1.772e-5;
	double A21 = 9e-5;
	double A31 = 1.e-10;
	double A32 = 1.7e-5;
	//Eij
	double E21 = 2.2771e2 * PARAM::KBOLTZ_CGS;
	double E32 = 9.886e1 * PARAM::KBOLTZ_CGS;
	double E31 = E21 + E32;

	//Hydrogen 1-0
	double C21_H = 9.2e-11 * pow(T / 100, .67);
	//electron
	double C21_e = 1.4e-8;
	//2-0
	double C31_H = 4.3e-11 * pow(T / 100, .80);
	double C31_e = 1.4e-8;
	//2-1
	double C32_H = 1.1e-10 * pow(T / 100, .44);
	double C32_e = 5.0e-9;

	double C21_H2 = 0.0;
	double C21_Hp = 0.0;
	double C31_H2 = 0.0;
	double C31_Hp = 0.0;
	double C32_H2 = 0.0;
	double C32_Hp = 0.0;
	double tinv = 1. / T;
	double tinth = pow(tinv, 1 / 3);
	double tinq = pow(tinv, 0.25);
	double tintq = pow(tinv, 0.75);
	double tsqrt = sqrt(T);
	double tisqt = 1. / tsqrt;
	double tlog = log10(T);
	double t4log = tlog - 4.0;
	if (T < 50) {
		C32_H = 2.62e-12 * pow(T, 0.74);
	}
	if (T < 5) {
		double tfix = 5;
		double tfintq = 1 / pow(tfix, 0.75);
		C21_H = (5e-11 / 3)
				* exp(
						4.581 - 156.118 * tfintq + 2679.979 * pow(tfintq, 2)
								- 78996.962 * pow(tfintq, 3)
								+ 1308323.468 * pow(tfintq, 4)
								- 13011761.861 * pow(tfintq, 5)
								+ 71010784.971 * pow(tfintq, 6)
								- 162826621.855 * pow(tfintq, 7))
				* exp(E21 / (PARAM::KBOLTZ_CGS * tfix));
		C31_H = 5e-11
				* exp(
						3.297 - 168.382 * tfintq + 1844.099 * pow(tfintq, 2)
								- 68362.889 * pow(tfintq, 3)
								+ 1376864.737 * pow(tfintq, 4)
								- 17964610.169 * pow(tfintq, 5)
								+ 134374927.808 * pow(tfintq, 6)
								- 430107587.886 * pow(tfintq, 7))
				* exp(E31 / (PARAM::KBOLTZ_CGS * tfix));
	} else if (T < 1e3) {
		C21_H = (5e-11 / 3)
				* exp(
						4.581 - 156.118 * tintq + 2679.979 * pow(tintq, 2)
								- 78996.962 * pow(tintq, 3)
								+ 1308323.468 * pow(tintq, 4)
								- 13011761.861 * pow(tintq, 5)
								+ 71010784.971 * pow(tintq, 6)
								- 162826621.855 * pow(tintq, 7))
				* exp(E21 / (PARAM::KBOLTZ_CGS * T));
		C31_H = 5e-11
				* exp(
						3.297 - 168.382 * tintq + 1844.099 * pow(tintq, 2)
								- 68362.889 * pow(tintq, 3)
								+ 1376864.737 * pow(tintq, 4)
								- 17964610.169 * pow(tintq, 5)
								+ 134374927.808 * pow(tintq, 6)
								- 430107587.886 * pow(tintq, 7))
				* exp(E31 / (PARAM::KBOLTZ_CGS * T));
		C32_H = 3e-11
				* exp(
						3.437 + 17.443 * tisqt - 618.761 * pow(tisqt, 2)
								+ 3757.156 * pow(tisqt, 3)
								- 12736.468 * pow(tisqt, 4)
								+ 22785.266 * pow(tisqt, 5)
								- 22759.228 * pow(tisqt, 6)
								+ 12668.261 * pow(tisqt, 7))
				* exp(E21 / (PARAM::KBOLTZ_CGS * T));
	} else {
		C21_H = 6.81e-11 * pow(T, 0.376);
		C31_H = 6.34e-11 * pow(T, 0.36);
		C32_H = 3.61e-10 * pow(T, 0.158);

	}
	double tfix;
	if (T > 1e4) {
		tfix = 1e4;
	} else {
		tfix = T;
	}
	double opratio = 2.4;
	double fortho = opratio / (1. + opratio);
	double fpara = 1. - fortho;
	C21_H2 = fortho * 2.70e-11 * pow(tfix, 0.362)
			+ fpara * 3.46e-11 * pow(tfix, 0.316);
	C31_H2 = fortho * 5.49e-11 * pow(tfix, 0.317)
			+ fpara * 7.07e-11 * pow(tfix, 0.268);
	C32_H2 = fortho * 2.74e-14 * pow(tfix, 1.060)
			+ fpara * 3.33e-15 * pow(tfix, 1.360);
	C21_e = 5.12e-10 * pow(T, (-0.075));
	C31_e = 4.863e-10 * pow(T, (-0.026));
	C32_e = 1.082e-14 * pow(T, (0.926));

	if (T < 193) {
		C21_Hp = 6.38e-11 * pow(T, 0.40);
	} else if (T < 3686) {
		C21_Hp = 7.75e-12 * pow(T, 0.80);
	} else {
		C21_Hp = 2.65e-10 * pow(T, 0.37);
	}

	if (T < 511) {
		C31_Hp = 6.10e-13 * pow(T, 1.10);
	} else if (T < 7510) {
		C31_Hp = 2.12e-12 * pow(T, 0.90);
	} else {
		C31_Hp = 4.49e-10 * pow(T, 0.30);
	}
	if (T < 2090) {
		C32_Hp = 2.029e-11 * pow(T, 0.56);
	} else {
		C32_Hp = 3.434e-10 * pow(T, 0.19);
	}

	//	double C21 = C21_H * NUMDENS_CGS_H + C21_H2 * NUMDENS_CGS_H2 + C21_e * NUMDENS_CGS_e + C21_Hp * NUMDENS_CGS_Hp;
	//	double C31 = C31_H * NUMDENS_CGS_H + C31_H2 * NUMDENS_CGS_H2 + C31_e * NUMDENS_CGS_e + C31_Hp * NUMDENS_CGS_Hp;
	//	double C32 = C32_H * NUMDENS_CGS_H + C32_H2 * NUMDENS_CGS_H2 + C32_e * NUMDENS_CGS_e + C32_Hp * NUMDENS_CGS_Hp;

	C21_e = 1.4e-8;
	C31_e = 1.4e-8;
	C32_e = 5.0e-9;
	C21_H = 9.2e-11 * pow(T / 100, 0.67);
	C31_H = 4.3e-11 * pow(T / 100, 0.8);
	C32_H = 1.1e-10 * pow(T / 100, 0.44);
	double C21 = C21_H * NUMDENS_CGS_H + C21_e * NUMDENS_CGS_e
			+ C21_Hp * NUMDENS_CGS_Hp;
	double C31 = C31_H * NUMDENS_CGS_H + C31_e * NUMDENS_CGS_e
			+ C31_Hp * NUMDENS_CGS_Hp;
	double C32 = C32_H * NUMDENS_CGS_H + C32_e * NUMDENS_CGS_e
			+ C32_Hp * NUMDENS_CGS_Hp;
	double cooling = 0.;

	double g1 = 5, g2 = 3, g3 = 1;

	double C12 = (g2 / g1) * C21 * exp(-E21 / (PARAM::KBOLTZ_CGS * T));
	double C13 = (g3 / g1) * C31 * exp(-E31 / (PARAM::KBOLTZ_CGS * T));
	double C23 = (g3 / g2) * C32 * exp(-E32 / (PARAM::KBOLTZ_CGS * T));
	double B21, B31, B32, B12, B13, B23 = 0.0;

	double x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	double _alpha, _beta, _gamma, _delta, _theta;
	if (x < 5.0) {
		B21 = 0.0;
	} else {
		B21 = A21 / (exp(x) - 1.);
	}
	B12 = (g2 / g1) * B21;

	x = E31 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	if (x < 5.0) {
		B31 = 0.0;
	} else {
		B31 = A31 / (exp(x) - 1.);
	}

	B13 = (g3 / g1) * B31;
	if (x < 5.0) {
		B32 = 0.0;
	} else {
		x = E32 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	}
	B32 = A32 / (exp(x) - 1.);
	B23 = (g3 / g2) * B32;

	double Beta_esc = 1.0, N1, N2, N3;
	double R12 = C12 + B12;
	double R13 = C13 + B13;
	double R23 = C23 + B23;
	double R21 = C21 + A21 + B21;
	double R31 = C31 + A31 + B31;
	double R32 = C32 + A32 + B32;
	three_level_population(R12, R13, R23, R21, R31, R32, N1, N2, N3);
	if (abundance_OI < 1e-9) {
		cooling = 0.0;
	} else {
		//		double oxlam1 = (A21 + B21) * E21 * N2 + ((A31 + B31) * E31 + (A32 + B32) * E32) * N3;
		double oxlam1 = (A21 + B21) * E21 * N2;
		double oxlam2 = (B12 * E21 + B13 * E31) * N1 + B23 * E32 * N2;
		cooling = NUMDENS_CGS * (oxlam1 - oxlam2) * abundance_OI;
	}
	return cooling;
}

double Lambda_CII(double T, double NUMDENS_CGS, double abundance_H,
		double abundance_H2, double abundance_e, double abundance_CII,
		double abundance_C, double tau_init, double col_dens_H2, double fpara,
		double fortho) {
	double NUMDENS_CGS_H = abundance_H * NUMDENS_CGS;
	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double NUMDENS_CGS_Cp = abundance_CII * NUMDENS_CGS;
	double z = PARAM::col_dens_H2 / NUMDENS_CGS;
	double A21 = 2.291e-6;
	double E21 = 91.25 * PARAM::KBOLTZ_CGS;
	double C21 = 0.;
	double B21 = 0;
	double B12 = 0;
	double Q21 = 0.0;
	double _alpha = 0;
	double _beta = 0;
	double var1 = 0;
	double var2 = 0;
	double N1 = 0;
	double N2 = 0;

	double cooling = 0;
	double C21_H = 0.;
	double C21_e = 0;
	double temp2 = T * 1e-2;
	double hh = 0;

	double C12 = 0;
	double x, g1, g2;
	double Beta_esc = 1.0;
//	tau = z;
	g1 = 1;
	g2 = 2;
	double lambda21 = 157.7e-4;
	double C21_H2_Ortho = 0.0;
	double C21_H2_Para = 0.0, C21_H2 = 0.0;
//	if (tau_init < 7.0) {
//		Beta_esc = (1.0 - exp(-2.34 * tau_init)) / (4.68 * tau_init);
//	} else {
//		Beta_esc = 1.0 / (4.0 * tau_init * sqrt(log(tau_init / sqrt(M_PI))));
//	}
//	Beta_esc = 1 / (1 + .5 * tau_init);
	C21_H = T < 2.e3 ?
			8.e-10 * pow(temp2, 0.07) : 3.11361e-10 * pow(temp2, 0.385);
	C21_H2_Ortho =
			(T < 250.) ?
					(4.7e-10 + 4.6e-13 * T) * fortho :
					(5.85e-10 * pow(T, 0.07)) * fortho;
	C21_H2_Para =
			(T < 250.) ?
					(2.5e-10 * pow(T, 0.12)) * fpara :
					(4.85e-10 * pow(T, 0.07)) * fpara;
	C21_H2 = C21_H2_Ortho + C21_H2_Para;

	C21_e = (T < 2.e3) ?
			3.86e-7 * pow(temp2, -.5) : 2.426206e-7 / pow(temp2, 0.345);

	C21 = C21_H * NUMDENS_CGS_H + C21_e * NUMDENS_CGS_e
			+ C21_H2 * NUMDENS_CGS_H2;

	C12 = C21 * (g2 / g1) * exp(-E21 / (PARAM::KBOLTZ_CGS * T));
	x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B21 = A21 / (exp(x) - 1.);
	B12 = (g2 / g1) * B21;
	_alpha = (C21 + A21 * Beta_esc + B21 * Beta_esc);
	_beta = (C12 + B12 * Beta_esc);
	N1 = _alpha / (_alpha + _beta);
	N2 = _beta / (_alpha + _beta);
////////////////////////////////////////////////////////////
//
//		double tau = abundance_Cp * (g2 / g1) * (A21 * pow(lambda21, 3) / (4 * pow(2 * M_PI, 1.5) * PARAM::snds * 1e5))
//				* PARAM::col_dens_H2 * N1 * (1.0 - N2 * g1 / (N1 * g2));
//		double Beta_esc_old = Beta_esc;
////	Beta_esc = 1 / (1 + .5 * tau);
////	double iter = 0;
////	while (fabs(Beta_esc - Beta_esc_old) > pow(10, -5)) {
////		Beta_esc_old = Beta_esc;
////		double beta = pop_CII(T, NUMDENS_CGS, abundance_H, abundance_H2, abundance_e, abundance_Cp, abundance_C, tau, N1,
////				N2, B21, B12);
////		tau = abundance_Cp * (g2 / g1) * (A21 * pow(lambda21, 3) / (4 * pow(2 * M_PI, 1.5) * PARAM::snds * 1e5))
////				* PARAM::col_dens_H2 * N1 * (1.0 - N2 * g1 / (N1 * g2));
////		Beta_esc = 1 / (1 + .5 * tau);
////
////		iter = iter + 1;
////	}
/////////////////////////////////////////////////////////////
	var1 = (A21 + B21) * E21 * N2;
	var2 = B12 * E21 * N1;
	cooling = (var1 - var2) * abundance_CII * NUMDENS_CGS;
//		cooling /= 10;
//
//		double n_h = abundance_H * NUMDENS_CGS + abundance_H2 * NUMDENS_CGS;
//		double f = 2.0 * abundance_H2 * NUMDENS_CGS / n_h;
//		double n_cp = abundance_Cp * NUMDENS_CGS;
//		double N_d = 6.6e20;
//		double _tau = col_dens_H2 / N_d;
//		double eps = 1.0 / (1.0 + _tau * pow(2 * M_PI * log(2.13 + _tau * _tau), 1.5));
//		double n_crit = 8.7 * sqrt(T / 100) + 3e2 * pow(T / 100, -0.07);
//		cooling = n_h * (1.0 - .5 * f) * n_cp * C12 * E21
//				/ (1.0 + (1.0 + (g2 / g1) * exp(-E21 / ( PARAM::KBOLTZ_CGS*T))) * (NUMDENS_CGS / (n_crit * eps)));
//	} else {

//	C21_H = 8.e-10 * pow(temp2, 0.07);
//	C21_e = 2.86e-7 * pow(temp2, -.5);
//	C21 = C21_H + C21_e / NUMDENS_CGS / NUMDENS_CGS_e;
//	C12 = C21 * (g2 / g1) * exp(-E21 / (PARAM::KBOLTZ_CGS * T));
//	double n_h = abundance_H * NUMDENS_CGS + abundance_H2 * NUMDENS_CGS;
//	double f = 2.0 * abundance_H2 * NUMDENS_CGS / n_h;i
//	double n_cp = abundance_Cp * NUMDENS_CGS;
//	double N_d = 6.5e20;
//	double _tau = col_dens_H2 / N_d;
//	double eps = 1.0 / (1.0 + _tau * sqrt(2 * M_PI * log(2.13 + _tau * _tau)));
//	double n_crit = 8.7 * sqrt(T / 100) + 3e2 * pow(T / 100, -0.07);
//	cooling = n_h * (1.0 - .5 * f) * n_cp * C12 * E21
//			/ (1.0 + (1.0 + (g2 / g1) * exp(-E21 / (PARAM::KBOLTZ_CGS * T))) * (NUMDENS_CGS / (n_crit * eps)));
//		double n_h = abundance_H * NUMDENS_CGS + abundance_H2 * NUMDENS_CGS;
//				double f = 2.0 * abundance_H2 * NUMDENS_CGS / n_h;
//				double n_cp = abundance_Cp * NUMDENS_CGS;
//				double N_d = 6.6e20;
//				double _tau = col_dens_H2 / N_d;
//				double n_crit = 8.7 * sqrt(T / 100) + 3e2 * pow(T / 100, -0.07);
//		cooling = 2.0*2.2e-23*NUMDENS_CGS * n_cp * exp(-92/T)
//						/ (1.0 + (1.0 + (g2 / g1) * exp(-E21 / ( PARAM::KBOLTZ_CGS*T))) * (NUMDENS_CGS / (n_crit )));
//		cooling = 7.9e-20*pow(T,-.5)*n_h * (1.0 - .5 * f) * n_cp * exp(-92/T);
//		double OmegaT = 1.8 + .484 * T / 10000 + 4.01 * pow(T / 10000, 2) - 3.39 * pow(T / 10000, 3);
//		cooling = 2.54e-14 * abundance_Cp
//				* ((2.1e-7 / sqrt(T / 100)) * OmegaT * abundance_e * NUMDENS_CGS + 8.86e-10 * NUMDENS_CGS * abundance_H)
//				* exp(-92 / T);
//	}
	return cooling;

}

double Lambda_Meta_CII(double T, double NUMDENS_CGS, double abundance_H,
		double abundance_H2, double abundance_e, double abundance_Cp,
		double abundance_C, double tau) {
	double NUMDENS_CGS_H = abundance_H * NUMDENS_CGS;
	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double NUMDENS_CGS_Cp = abundance_Cp * NUMDENS_CGS;

	double A21 = 3.6;
	double E21 = 6.2e4 * PARAM::KBOLTZ_CGS;
	double C21 = 0.;
	double B21 = 0;
	double B12 = 0;
	double Q21 = 0.0;
	double _alpha = 0;
	double _beta = 0;
	double var1 = 0;
	double var2 = 0;
	double N1 = 0;
	double N2 = 0;

	double cooling = 0;
	double C21_H = 0.;
	double C21_e = 0;
	double temp2 = T * 1e-2;
	double f = 0;
	double hh = 0;

	double C12 = 0;
	double x, g1, g2;
	double Beta_esc = 1.0;

	g1 = 6;
	g2 = 12;
//	if (tau < 7.0) {
//		Beta_esc = 1.0 - exp(-2.34 * tau) / (4.68 * tau);
//	} else {
//		Beta_esc = 1.0 / (4.0 * tau * sqrt(log(tau / sqrt(M_PI))));
//	}
	C21_H = 1e-12;
	C21_e = 2.3e-8 * pow(T / 1e4, -.5);
	C21 = C21_H * NUMDENS_CGS_H + C21_e * NUMDENS_CGS_e;
	C12 = C21 * (g2 / g1) * exp(-E21 / (PARAM::KBOLTZ_CGS * T));

	x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B21 = A21 / (exp(x) - 1.);
	B12 = (g2 / g1) * B21;
	_alpha = (C21 + A21 * Beta_esc + B21 * Beta_esc);
	_beta = (C12 + B12 * Beta_esc);
	N1 = _alpha / (_alpha + _beta);
	N2 = _beta / (_alpha + _beta);
	var1 = (A21 + B21) * E21 * N2;
	var2 = B12 * E21 * N1;

	cooling = Beta_esc * (var1 - var2) * (g2 / g1) * exp(-E21 / T)
			* abundance_Cp * NUMDENS_CGS;
//	cooling = NUMDENS_CGS * abundance_Cp * NUMDENS_CGS * abundance_e
//			* (3.e-17 * exp(-61900 / T) / sqrt(T));
	return cooling;

}
double Lambda_SiII(double T, double NUMDENS_CGS, double abundance_H,
		double abundance_H2, double abundance_e, double abundance_Sip,
		double tau) {

	double NUMDENS_CGS_H = abundance_H * NUMDENS_CGS;
	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double NUMDENS_CGS_Sip = abundance_Sip * NUMDENS_CGS;

	double A21 = 2.1e-4;
	double E21 = 410 * PARAM::KBOLTZ_CGS;
	double C21 = 0.;
	double B21 = 0;
	double B12 = 0;
	double Q21 = 0.0;
	double _alpha = 0;
	double _beta = 0;
	double N1 = 0;
	double N2 = 0;

	double cooling = 0;
	double C21_H = 0.;
	double C21_e = 0;
	double temp2 = T * 1e-2;
	double f = 0;
	double hh = 0;
	double var1, var2;
	double C12 = 0;
	double x, g1, g2;
	double Beta_esc = 1.0;

	g1 = 1;
	g2 = 2;
//	if (tau < 7.0) {
//		Beta_esc = 1.0 - exp(-2.34 * tau) / (4.68 * tau);
//	} else {
//		Beta_esc = 1.0 / (4.0 * tau * sqrt(log(tau / sqrt(M_PI))));
//	}
	C21_H = 8.e-10 * pow(temp2, 0.07);
	C21_e = 1.7e-6 * pow(temp2, -.5);
	C21 = C21_H * NUMDENS_CGS_H + C21_e * NUMDENS_CGS_e;
	C12 = C21 * (g2 / g1) * exp(-E21 / (PARAM::KBOLTZ_CGS * T));

	x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B21 = A21 / (exp(x) - 1.);
	B12 = (g2 / g1) * B21;
	_alpha = (C21 + A21 * Beta_esc + B21 * Beta_esc);
	_beta = (C12 + B12 * Beta_esc);
	N1 = _alpha / (_alpha + _beta);
	N2 = _beta / (_alpha + _beta);
	var1 = (A21 + B21) * E21 * N2;
	var2 = B12 * E21 * N1;

	cooling = Beta_esc * (var1 - var2) * (g2 / g1) * exp(-E21 / T)
			* abundance_Sip * NUMDENS_CGS;
//	cooling = NUMDENS_CGS * abundance_e * NUMDENS_CGS * abundance_Sip
//			* (exp(-413 / T)) * 1.9e-18 / sqrt(T)
//			+ NUMDENS_CGS * abundance_H * NUMDENS_CGS * abundance_Sip
//					* (exp(-413 / T)) * 7.4e-23;
	return cooling;

}

double Lambda_Meta_SiII(double T, double NUMDENS_CGS, double abundance_HI,
		double abundance_H2, double abundance_e, double abundance_SiII,
		double tau) {

	double NUMDENS_CGS_H = abundance_HI * NUMDENS_CGS;
	double NUMDENS_CGS_H2 = abundance_H2 * NUMDENS_CGS;
	double NUMDENS_CGS_e = abundance_e * NUMDENS_CGS;
	double NUMDENS_CGS_Sip = abundance_SiII * NUMDENS_CGS;

	double A21 = 6.4e3;
	double E21 = 6.2e4 * PARAM::KBOLTZ_CGS;
	double C21 = 0.;
	double B21 = 0;
	double B12 = 0;
	double Q21 = 0.0;
	double LHS = 0;
	double RHS = 0;
	double var1 = 0;
	double var2 = 0;
	double N1 = 0;
	double N2 = 0;

	double cooling = 0.;
	double C21_H = 0.;
	double C21_e = 0;
	double temp2 = T * 1e-2;
	double f = 0;
	double hh = 0;

	double C12 = 0;
	double x, g1, g2;
	double Beta_esc = 1.0;

	g1 = 6;
	g2 = 12;
//	if (tau < 7.0) {
//		Beta_esc = 1.0 - exp(-2.34 * tau) / (4.68 * tau);
//	} else {
//		Beta_esc = 1.0 / (4.0 * tau * sqrt(log(tau / sqrt(M_PI))));
//	}
	Beta_esc = 1.0;
	C21_H = 1e-12;
	C21_e = 6.5e-8 * pow(T / 1e4, -.5);
	C21 = C21_H * NUMDENS_CGS_H + C21_e * NUMDENS_CGS_e;
	C12 = C21 * (g2 / g1) * exp(-E21 / (PARAM::KBOLTZ_CGS * T));

	x = E21 / (PARAM::KBOLTZ_CGS * PARAM::CMB_T);
	B21 = A21 / (exp(x) - 1.);
	B12 = (g2 / g1) * B21;
	LHS = (C21 + A21 * Beta_esc + B21 * Beta_esc);
	RHS = (C12 + B12 * Beta_esc);
	N1 = LHS / (LHS + RHS);
	N2 = RHS / (LHS + RHS);
	var1 = (A21 + B21) * E21 * N2;
	var2 = B12 * E21 * N1;
	cooling = (var1 - var2) * (g2 / g1) * exp(-E21 / T) * abundance_SiII
			* NUMDENS_CGS;
//	cooling = NUMDENS_CGS * abundance_e * NUMDENS_CGS * abundance_SiII
//			* (exp(-63600 / T)) * 3.0e-17 / sqrt(T);
	return cooling;
}
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
double Lambda_H2_rovib(double T, double NUMDENS_CGS, double abundance_H2,
		double abundance_H) {
	double NUMDENS_CGS_H = NUMDENS_CGS * abundance_H;
	double NUMDENS_CGS_H2 = NUMDENS_CGS * abundance_H2;
	double E_10, E_20, Lambda_H2 = 0.0;
	if (abundance_H > 0.0 && abundance_H2 > 0.0 && T > 100) {
		double L_r_H, L_r_H2, L_v_H, L_v_H2, L_r_H2_0, L_r_H_0, L_v_H_0,
				L_v_H2_0, n_cr_H_rot_n_H, n_cr_H_vib_n_H, n_cr_H2_rot_n_H2,
				n_cr_H2_vib_n_H2, L_vr_H, L_vr_H2;
		double temp = T;
		double gamma_10_H, gamma_20_H, gamma_10_H2, gamma_20_H2;
		double E_10, E_20;

		L_r_H = (1. / NUMDENS_CGS_H)
				* (9.5e-22 * exp(-pow(.13 / (temp / 1000.), 3))
						* pow(temp / 1000., 3.76)
						/ (1. + .12 * pow(temp / 1000., 2.1))
						+ 3.e-24 * exp(-.51 / (temp / 1000.)));
		L_r_H2 = (1. / NUMDENS_CGS_H2)
				* (9.5e-22 * exp(-pow(.13 / (temp / 1000.), 3))
						* pow(temp / 1000., 3.76)
						/ (1. + .12 * pow(temp / 1000., 2.1))
						+ 3.e-24 * exp(-.51 / (temp / 1000.)));
		L_v_H = (1. / NUMDENS_CGS_H)
				* (6.7e-19 * exp(-5.86 / (temp / 1000.))
						+ 1.6e-18 * exp(-11.7 / (temp / 1000.)));
		L_v_H2 = (1. / NUMDENS_CGS_H2)
				* (6.7e-19 * exp(-5.86 / (temp / 1000.))
						+ 1.6e-18 * exp(-11.7 / (temp / 1000.)));
		L_r_H_0 = pow(10,
				-103.0 + 97.59 * log10(temp) - 48.05 * pow(log10(temp), 2)
						+ 10.8 * pow(log10(temp), 3)
						- .9032 * pow(log10(temp), 4));
		L_r_H2_0 = L_r_H_0;
		gamma_10_H = 1.e-12 * pow(temp, .5) * exp(-1000. / temp);
		gamma_20_H = 1.6e-12 * pow(temp, .5) * exp(-pow(400. / temp, 2));
		gamma_10_H2 = 1.4e-12 * pow(temp, .5) * exp(-12000. / (temp + 1200.));
		gamma_20_H2 = 0.;
		E_10 = 5860.;
		E_20 = 5860. * 2.;
		L_v_H_0 = gamma_10_H * exp(-E_10 / temp) * E_10 * PARAM::KBOLTZ_CGS
				+ gamma_20_H * exp(-E_20 / temp) * E_20 * PARAM::KBOLTZ_CGS;
		L_v_H2_0 = gamma_10_H2 * exp(-E_10 / temp) * E_10 * PARAM::KBOLTZ_CGS
				+ gamma_20_H2 * exp(-E_20 / temp) * E_20 * PARAM::KBOLTZ_CGS;
		n_cr_H_rot_n_H = L_r_H / L_r_H_0;
		n_cr_H_vib_n_H = L_v_H / L_v_H_0;
		n_cr_H2_rot_n_H2 = L_r_H2 / L_r_H2_0;
		n_cr_H2_vib_n_H2 = L_v_H2 / L_v_H2_0;
		L_vr_H = L_r_H / (1. + n_cr_H_rot_n_H) + L_v_H / (1. + n_cr_H_vib_n_H);
		L_vr_H2 = L_r_H2 / (1. + n_cr_H2_rot_n_H2)
				+ L_v_H2 / (1. + n_cr_H2_vib_n_H2);
		Lambda_H2 = abundance_H2
				* (abundance_H * L_vr_H + 2.0 * abundance_H2 * L_vr_H2);

		Lambda_H2 *= NUMDENS_CGS * NUMDENS_CGS;

	}
	return Lambda_H2;
}
// erg s-1///
double Lambda_CO_rot(double T, double NUMDENS_CGS, double abundance_CO) {
	double A_0 = 9.7e-8;
	double E_0 = 2.76 * PARAM::KBOLTZ_CGS;
	double n_cr = 3.3e6 * pow(T / 1000., 3. / 4.);
	double Lambda_CO_rot = 4. * pow(PARAM::KBOLTZ_CGS * T, 2) * A_0
			/ (E_0
					* (1. + n_cr / NUMDENS_CGS
							+ 1.5 * pow(n_cr / NUMDENS_CGS, .5)));
	Lambda_CO_rot *= abundance_CO;
	if (Lambda_CO_rot < 1e-30) {
		Lambda_CO_rot = 0.0;
	}
	return Lambda_CO_rot;
}
// erg s-1///
double Lambda_CO_vib(double T, double abundance_H, double abundance_H2,
		double abundance_CO) {
	double DeltaE_10 = 3080. * PARAM::KBOLTZ_CGS;
	double gamma_01_H = 3.e-12 * pow(T, .5) * exp(-pow(2000. / T, 3.43))
			* exp(-3080. / T) * abundance_H;
	double gamma_01_H2 = 4.3e-14 * T * exp(-pow(3.14e5 / T, 0.333))
			* exp(-3080. / T) * abundance_H2;
	double Lambda_CO_vib = gamma_01_H * DeltaE_10 + gamma_01_H2 * DeltaE_10;
	Lambda_CO_vib *= abundance_CO;
	return Lambda_CO_vib;
}

double recomb_HpeGrain_HGrain_rate(double T, double numdens, double abundance_e,
		double dust_to_gas_ratio) {

//	double shielding_column = numdens * shielding_length;
//	double AV = AV_conversion_factor * dust_to_gas_ratio * shielding_column;

//	double G_dust = isrf_chi * exp(-2.5 * AV);
	double fshield = 1.;
	double f_dust = 1.;
	double G_dust = 1.;
	double phi = 0.0;
	double k_hp_recomb_dust;
	if (abundance_e == 0.) {
		phi = 1e20;
	} else {
		phi = G_dust * pow(T, .5) / (numdens * abundance_e);
	}
	if (phi <= 1.e-6) {
		k_hp_recomb_dust = 1.225e-13 * dust_to_gas_ratio;
	} else {
		double temp_dep_1 = 5.087e2 * pow(T, 1.586e-2);

		double temp_dep_2 = -0.4723 - 1.102e-5 * log10(T);
		double hgrvar1 = 8.074e-6 * pow(phi, 1.378);
		double hgrvar2 = (1. + temp_dep_1 * pow(phi, temp_dep_2));
		k_hp_recomb_dust = 1.225e-13 * dust_to_gas_ratio
				/ (1. + hgrvar1 * hgrvar2);
	}
	return k_hp_recomb_dust;

}

PS::F64 rate_Gamma(PS::F64 T, PS::F64 NUMDENS_CGS, PS::F64 abundance_HI,
		PS::F64 abundance_e, PS::F64 abundance_H2, PS::F64 abundance_HII,
		PS::F64 abundance_CII, PS::F64 dust_to_gas_ratio, PS::F64 T_gr,
		PS::F64 col_den_H2, PS::F64 GAMMA, PS::F64 mu) {
	double G1 = Gamma_UV(T, NUMDENS_CGS, abundance_H2, abundance_HI,
			PARAM::col_dens_H2, mu) * NUMDENS_CGS;
	double G2 = Gamma_H2_grain_forming(T, T_gr, NUMDENS_CGS, abundance_H2,
			abundance_HI, col_den_H2);

	double G3 = Gamma_XR(abundance_e, PARAM::col_dens_H2) * NUMDENS_CGS;
	double G4 = Gamma_pe(T, NUMDENS_CGS, abundance_e) * NUMDENS_CGS;
	double G5 = Gamma_CR(NUMDENS_CGS, abundance_e) * NUMDENS_CGS;
	double G6 = Gamma_H2_associative_forming(T, NUMDENS_CGS, abundance_H2,
			abundance_HI, abundance_HII, abundance_e);
	PS::F64 tot = (G1 < PARAM::OVERFLOW_LIMIT ? 0.0 : G1)
			+ (G2 < PARAM::OVERFLOW_LIMIT ? 0.0 : G2)
			+ (G3 < PARAM::OVERFLOW_LIMIT ? 0.0 : G3)
			+ (G4 < PARAM::OVERFLOW_LIMIT ? 0.0 : G4)
			+ (G5 < PARAM::OVERFLOW_LIMIT ? 0.0 : G5)
			+ (G6 < PARAM::OVERFLOW_LIMIT ? 0.0 : G6);
//	if (tot != tot) {
//		std::cout << "Tmp" << T << " UV: " << ((G1)) << " " << (G2) << " H2:" << (G2 + G1) << " XR:" << G3 << " PE:" << G4 << " CR:" << G5 << std::endl;
//	}
	return tot;
}
PS::F64 rate_Lambda(PS::F64 T, PS::F64 NUMDENS_CGS, PS::F64 abundance_HI,
		PS::F64 abundance_e, PS::F64 abundance_H2, PS::F64 abundance_HII,
		PS::F64 abundance_CII, PS::F64 abundance_FeII, PS::F64 abundance_SiII,
		PS::F64 abundance_CI, PS::F64 abundance_OI, PS::F64 abundance_CO,
		PS::F64 dust_to_gas_ratio, PS::F64 T_gr, PS::F64 grain_size) {
	double L1 = Lambda_Lya(T, abundance_e, NUMDENS_CGS) * NUMDENS_CGS;
	double L2 = Lambda_CII(T, NUMDENS_CGS, abundance_HI, abundance_H2,
			abundance_e, abundance_CII, abundance_CI, PARAM::tau,
			PARAM::col_dens_H2, PARAM::fpara, PARAM::fortho);
	double L3 = Lambda_OI(T, NUMDENS_CGS, abundance_HI, abundance_H2,
			abundance_e, abundance_HII, abundance_OI, PARAM::tau,
			PARAM::col_dens_H2);
	double L4 = Lambda_Recomb_Grain_PAH(T, abundance_e, NUMDENS_CGS)
			* NUMDENS_CGS;
	double L5 = Lambda_Col_Dust(T, abundance_HI, abundance_e, NUMDENS_CGS,
			grain_size, T_gr);
	double L6 = Lambda_CO_rot(T, NUMDENS_CGS, abundance_CO) * NUMDENS_CGS;
	double L7 = Lambda_H2_rovib(T, NUMDENS_CGS, abundance_H2, abundance_HI);
	double L8 = Lambda_CO_vib(T, abundance_HI, abundance_H2, abundance_CO)
			* NUMDENS_CGS;

	double L9 = Lambda_SiII(T, NUMDENS_CGS, abundance_HI, abundance_H2,
			abundance_e, abundance_SiII, PARAM::tau);
	double L10 = Lambda_FeII(T, NUMDENS_CGS, abundance_HI, abundance_H2,
			abundance_e, abundance_HII, abundance_FeII, PARAM::tau);
	double L11 = Lambda_Meta_CII(T, NUMDENS_CGS, abundance_HI, abundance_H2,
			abundance_e, abundance_CII, abundance_CI, PARAM::tau);
	double L12 = Lambda_Meta_OI(T, NUMDENS_CGS, abundance_HI, abundance_H2,
			abundance_e, abundance_HII, abundance_OI, PARAM::tau);
	double L13 = Lambda_Meta_SiII(T, NUMDENS_CGS, abundance_HI, abundance_H2,
			abundance_e, abundance_SiII, PARAM::tau);
	double L14 = Lambda_Meta_FeII(T, NUMDENS_CGS, abundance_HI, abundance_H2,
			abundance_e, abundance_HII, abundance_FeII, PARAM::tau);
	double L15 = Lambda_Mol(T, abundance_e, NUMDENS_CGS * abundance_H2);
	PS::F64 tot = 0.0;
#ifdef KI2000
	tot = (L1 <PARAM:: OVERFLOW_LIMIT ? 0.0 : L1)
	+ (L4 <PARAM:: OVERFLOW_LIMIT ? 0.0 : L4)
	+ (L7 <PARAM:: OVERFLOW_LIMIT ? 0.0 : L7)
	+ (L2 <PARAM:: OVERFLOW_LIMIT ? 0.0 : L2)
	+ (L3 <PARAM:: OVERFLOW_LIMIT ? 0.0 : L3)
	+ (L5 <PARAM:: OVERFLOW_LIMIT ? 0.0 : L5)
	+ (L6 <PARAM:: OVERFLOW_LIMIT ? 0.0 : L6)
	+ (L8 <PARAM:: OVERFLOW_LIMIT ? 0.0 : L8);

#else
	tot = (L1 < PARAM::OVERFLOW_LIMIT ? 0.0 : L1)
			+ (L4 < PARAM::OVERFLOW_LIMIT ? 0.0 : L4)
			+ (L7 < PARAM::OVERFLOW_LIMIT ? 0.0 : L7)
			+ (L2 < PARAM::OVERFLOW_LIMIT ? 0.0 : L2)
			+ (L3 < PARAM::OVERFLOW_LIMIT ? 0.0 : L3)
			+ (L5 < PARAM::OVERFLOW_LIMIT ? 0.0 : L5)
			+ (L6 < PARAM::OVERFLOW_LIMIT ? 0.0 : L6)
			+ (L8 < PARAM::OVERFLOW_LIMIT ? 0.0 : L8);
	if (tot != tot) {
		std::cout << " Lya:" << ((L1) / NUMDENS_CGS) << " CII:"
				<< ((L2) / NUMDENS_CGS) << " CO:" << ((L6 + L8) / NUMDENS_CGS)
				<< " OI:" << (L3 / NUMDENS_CGS) << " Rec:" << (L4 / NUMDENS_CGS)
				<< " GR:" << (L5 / NUMDENS_CGS) << "H2 rovib"
				<< (L7 / NUMDENS_CGS) << "Tot " << tot << std::endl;
	}

#endif
	if (T < 5.0) {
		//	tot = 0.0;
	}
	return tot;
}

//cm-3 s-1///
PS::F64 _C_CO(PS::F64 NUMDENS_CGS, PS::F64 abundance_OI, PS::F64 abundance_CII,
		PS::F64 abundance_H2) {
	PS::F64 k0 = 5.e-16;
	PS::F64 k1 = 5.e-10;
	PS::F64 gamma_chx = k1 * PARAM::G0;
	PS::F64 G_co = 1.e-10 * PARAM::G0;
	PS::F64 beta = k1 * (abundance_OI)
			/ (k1 * (abundance_OI)
					+ PARAM::G0 * gamma_chx / (abundance_H2 * NUMDENS_CGS));

	PS::F64 ydot = k0 * abundance_CII * beta * NUMDENS_CGS * NUMDENS_CGS;
	return ydot;

}
// s-1///
PS::F64 _D_CO() {
	PS::F64 G_co = 1.e-10 * PARAM::G0;
	PS::F64 ydot = G_co;
	return ydot;

}
void evolve_abundance(RealPtcl &hydro, const PS::F64 oneDynTimeStep) {
	PS::F64 abundance_e = hydro.abundances[0];
	PS::F64 abundance_HI = hydro.abundances[1];
	PS::F64 abundance_HeI = hydro.abundances[2];
	PS::F64 abundance_OI = hydro.abundances[3];
	PS::F64 abundance_CI = hydro.abundances[4];
	PS::F64 abundance_H2 = hydro.abundances[5];
	PS::F64 abundance_CO = hydro.abundances[6];
	PS::F64 abundance_HII = hydro.abundances[7];
	PS::F64 abundance_CII = hydro.abundances[8];
	PS::F64 abundance_FeII = hydro.abundances[9];
	PS::F64 abundance_SiII = hydro.abundances[10];
	PS::F64 abundance_e_old = abundance_e;
	PS::F64 abundance_HI_old = abundance_HI;
	PS::F64 abundance_HeI_old = abundance_HeI;
	PS::F64 abundance_OI_old = abundance_OI;
	PS::F64 abundance_CI_old = abundance_CI;
	PS::F64 abundance_H2_old = abundance_H2;
	PS::F64 abundance_CO_old = abundance_CO;
	PS::F64 abundance_HII_old = abundance_HII;
	PS::F64 abundance_CII_old = abundance_CII;
	PS::F64 abundance_FeII_old = abundance_FeII;
	PS::F64 abundance_SiII_old = abundance_SiII;
	double mu = abundance_HI + abundance_HII + 2.0 * abundance_H2
			+ 4.0 * abundance_HeI;
	mu /= (abundance_HI + abundance_HII + abundance_H2 + abundance_HeI
			+ abundance_e);
	PS::F64 NUMDENS_CGS = hydro.dens * PARAM::SMassDens
			/ (mu * PARAM::PROTONMASS);
	PS::F64 MASSDENS_CGS = hydro.dens * PARAM::SMassDens;

	PS::F64 tdust = PARAM::Grain_T;
	PS::F64 dust_to_gas_ratio = PARAM::dust_to_gas_ratio;

	//	PS::F64 temprature_old = T;

	//	PS::F64 energy_old = hydro.eng * hydro.mass * PARAM::GSPH_SEng;
	//	energy_old/=(hydro.mass*PARAM::SM/(PARAM::PROTONMASS*PARAM::MU ));
	//	PS::F64 energy_init = hydro.eng * hydro.mass * PARAM::GSPH_SEng;
	//	energy_init/=(hydro.mass*PARAM::SM/(PARAM::PROTONMASS*PARAM::MU ));
	PS::F64 energy_old = hydro.eng * PARAM::SEng_per_Mass;

	PS::F64 energy_init = hydro.eng * PARAM::SEng_per_Mass;

	PS::F64 abundance_e_init = abundance_e;
	PS::F64 abundance_HI_init = abundance_HI;
	PS::F64 abundance_HeI_init = abundance_HeI;
	PS::F64 abundance_OI_init = abundance_OI;
	PS::F64 abundance_CI_init = abundance_CI;
	PS::F64 abundance_H2_init = abundance_H2;
	PS::F64 abundance_CO_init = abundance_CO;
	PS::F64 abundance_HII_init = abundance_HII;
	PS::F64 abundance_CII_init = abundance_CII;
	PS::F64 abundance_FeII_init = abundance_FeII;
	PS::F64 abundance_SiII_init = abundance_SiII;
	PS::F64 GAMMA_GAS = 5. / 3.;
	PS::F64 temprature_init = mu * PARAM::PROTONMASS_CGS * hydro.eng
			* (GAMMA_GAS - 1.0) * PARAM::SEng_per_Mass
			/ (PARAM::KBOLTZ_CGS * (1.1 + abundance_e - abundance_H2));

	PS::F64 ydot_H2 = 0.0;
	PS::F64 ydot_Hp = 0.0;
	PS::F64 ydot_CO = 0.0;

	//	PS::F64 Gamma = rate_Gamma(T, NUMDENS_CGS, abundance_HI, abundance_e, abundance_H2, abundance_HII, abundance_CII, dust_to_gas_ratio, PARAM::Grain_T, PARAM::col_dens_H2, GAMMA_GAS, mu);
	//	PS::F64 Lambda = rate_Lambda(T, NUMDENS_CGS, abundance_HI, abundance_e, abundance_H2, abundance_HII, abundance_CII, abundance_FeII, abundance_SiII, abundance_CI, abundance_OI, abundance_CO, dust_to_gas_ratio, PARAM::Grain_T, PARAM::grain_size);
	//	PS::F64 energy_dot = hydro.eng_dot * hydro.dens * PARAM::SM / (PARAM::SL * PARAM::ST * PARAM::ST * PARAM::ST);
	PS::F64 Gamma, Lambda;
	PS::F64 ylam = 0.0;

	PS::F64 energy = 0.0;
	PS::U64 main_try_steps = 1e8, sub_try_steps = 1000, chem_try_steps = 100;
	PS::F64 t_start = 0.;

	PS::F64 fac = 0.01;
	PS::F64 time_passed = 0.0, tleft = 0.0;
	PS::F64 dt = oneDynTimeStep * PARAM::ST;
	PS::F64 dt_substep = dt;
	PS::F64 abs_tol_H = 1e-7, rel_tol_H = fac, rel_diff_H;
	PS::F64 abs_tol_HII = 1e-7, rel_tol_HII = fac, rel_diff_HII;
	PS::F64 abs_tol_H2 = 1e-7, rel_tol_H2 = fac, rel_diff_H2;
	PS::F64 abs_tol_CII = 1e-7, rel_tol_CII = fac, rel_diff_CII;
	PS::F64 abs_tol_CO = 1e-7, rel_tol_CO = fac, rel_diff_CO;
	PS::F64 abs_tol_e = 1e-7, rel_tol_e = fac, rel_diff_e;
	PS::F64 abs_tol_eng, rel_tol_eng = fac, rel_eng_diff;
	double loop = 1000;
	double T = temprature_init;

	Gamma = rate_Gamma(T, NUMDENS_CGS, abundance_HI, abundance_e, abundance_H2,
			abundance_HII, abundance_CII, dust_to_gas_ratio, PARAM::Grain_T,
			PARAM::col_dens_H2, GAMMA_GAS, mu);
	Lambda = rate_Lambda(T, NUMDENS_CGS, abundance_HI, abundance_e,
			abundance_H2, abundance_HII, abundance_CII, abundance_FeII,
			abundance_SiII, abundance_CI, abundance_OI, abundance_CO,
			dust_to_gas_ratio, 8, 100);
	energy = T * (PARAM::KBOLTZ_CGS * (1.1 + abundance_e - abundance_H2))
			/ (mu * PARAM::PROTONMASS_CGS * (GAMMA_GAS - 1.0));
//	std::cout << (energy / fabs(Lambda)) * MASSDENS_CGS << std::endl;

	for (PS::U64 step = 0; step < main_try_steps; step++) {

//			std::cout << dt << " " << time_passed << " " << tleft << "out loop" << std::endl;
		PS::F32 iter_success = 0;
		tleft = oneDynTimeStep * PARAM::ST - time_passed;
		dt = (tleft < dt) ? tleft : dt;

		int reset1 = 1, loopcnt = 0;
		while (reset1) {
			loopcnt++;
			abundance_e = abundance_e_init;
			abundance_HI = abundance_HI_init;
			abundance_OI = abundance_OI_init;
			abundance_CI = abundance_CI_init;
			abundance_H2 = abundance_H2_init;
			abundance_CO = abundance_CO_init;
			abundance_HII = abundance_HII_init;
			abundance_CII = abundance_CII_init;
			abundance_FeII = abundance_FeII_init;
			abundance_SiII = abundance_SiII_init;
			T = temprature_init;
			energy = energy_init;
			//
			double K1 = RATE_CR_HII_HeII(abundance_e);
			double K2 = RATE_XR_HII_HeII(T, abundance_e, abundance_HI,
					PARAM::col_dens_H2);
			double K3 = RATE_H_col_HII(T);
			double K4 = RATE_HIIe_HPhoton(T);
			double K5 = RATE_HHGrain_H2(T, tdust);
			double K6 = RATE_HHm_H2_e(T, abundance_HI, abundance_HII);
			double K7 = RATE_H2_UV_2H(T, PARAM::col_dens_H2, mu, abundance_H2);
			double K8 = RATE_HH2_3H(T);
			double K9 = RATE_H2_CR_2H(abundance_e);

			double C_HII = K1 * abundance_HI * NUMDENS_CGS
					+ K2 * abundance_HI * NUMDENS_CGS
					+ K3 * NUMDENS_CGS * abundance_e * NUMDENS_CGS
							* abundance_HI;
			double D_HII = K4 * abundance_e * NUMDENS_CGS;
			double C_HI = K4 * abundance_e * abundance_HII * NUMDENS_CGS
					* NUMDENS_CGS + K7 * abundance_H2 * NUMDENS_CGS
					+ K8 * abundance_H2 * abundance_HI * NUMDENS_CGS
					+ K9 * abundance_H2 * NUMDENS_CGS;
			double D_HI = K1 + K2 + K3 * abundance_e * NUMDENS_CGS
					+ K5 * abundance_HI * NUMDENS_CGS;
			double C_H2 = K5 * abundance_HI * abundance_HI * NUMDENS_CGS
					* NUMDENS_CGS
					+ K6 * abundance_e * abundance_HI * NUMDENS_CGS
							* abundance_HI * NUMDENS_CGS;
			double D_H2 = K7 + K8 * NUMDENS_CGS * abundance_HI + K9;
			double C_CO = _C_CO(NUMDENS_CGS, abundance_OI, abundance_CII,
					abundance_H2);
			double D_CO = _D_CO();

			abundance_HII = C_HII / (NUMDENS_CGS) / D_HII
					+ (abundance_HII_init - C_HII / (NUMDENS_CGS) / D_HII)
							* exp(-D_HII * .5 * dt);
			abundance_HI = C_HI / (NUMDENS_CGS) / D_HI
					+ (abundance_HI_init - C_HI / (NUMDENS_CGS) / D_HI)
							* exp(-D_HI * .5 * dt);
			abundance_H2 = C_H2 / (NUMDENS_CGS) / D_H2
					+ (abundance_H2_init - C_H2 / (NUMDENS_CGS) / D_H2)
							* exp(-D_H2 * .5 * dt);
			abundance_CO = C_CO / (NUMDENS_CGS) / D_CO
					+ (abundance_CO_init - C_CO / (NUMDENS_CGS) / D_CO)
							* exp(-D_CO * .5 * dt);

			abundance_HII = abundance_HII
					/ (abundance_HII + abundance_HI + 2. * abundance_H2);
			abundance_HI = abundance_HI
					/ (abundance_HII + abundance_HI + 2. * abundance_H2);
			abundance_H2 = abundance_H2
					/ (abundance_HII + abundance_HI + 2. * abundance_H2);
			abundance_CII = fmax(abundance_CI - abundance_CO, 0.0);
			abundance_e = fmin(abundance_HII + abundance_CII + abundance_SiII,
					1.0);
			GAMMA_GAS = 5. / 3.;
			Gamma = rate_Gamma(T, NUMDENS_CGS, abundance_HI, abundance_e,
					abundance_H2, abundance_HII, abundance_CII,
					dust_to_gas_ratio, PARAM::Grain_T, PARAM::col_dens_H2,
					GAMMA_GAS, mu);
			Lambda = rate_Lambda(T, NUMDENS_CGS, abundance_HI, abundance_e,
					abundance_H2, abundance_HII, abundance_CII, abundance_FeII,
					abundance_SiII, abundance_CI, abundance_OI, abundance_CO,
					dust_to_gas_ratio, 8, 100);

			rel_eng_diff = (fabs(Gamma - Lambda) * dt)
					/ (energy * MASSDENS_CGS);

			if (fabs(rel_eng_diff) + 1e-5 > 0.2) {
//				if (loopcnt > 1000) {
//
//					std::cout << loopcnt << " " << rel_eng_diff << " "
//							<< fabs(rel_eng_diff - rel_tol_eng) / rel_tol_eng
//							<< std::endl;
//
//				}
				reset1 = 1;
				dt = fac * (energy * MASSDENS_CGS) / fabs(Gamma - Lambda);

			} else {
				reset1 = 0;
			}
		}
		energy = (energy_init + (Gamma - Lambda) * dt / MASSDENS_CGS);
//		if (energy < 0.0) {
//			std::cout << loopcnt << " " << rel_eng_diff << " "
//					<< fabs(rel_eng_diff - rel_tol_eng) / rel_tol_eng
//					<< std::endl;
//
//		}
		T = mu * PARAM::PROTONMASS_CGS * energy * (GAMMA_GAS - 1.0)
				* PARAM::SEng_per_Mass
				/ (PARAM::SEng_per_Mass * PARAM::KBOLTZ_CGS
						* (1.1 + abundance_e - abundance_H2));

		time_passed += dt;

//		if (time_passed >= oneDynTimeStep * PARAM::ST) {
////				std::cout << dt << " " << time_passed << " " << tleft << "out loop" << hydro.id << std::endl;
//
//			hydro.Lambda = Lambda * 1e26;
//
//			hydro.Gamma = Gamma * 1e26;
//			hydro.cooling_timescale = ((energy / fabs(Lambda)) * MASSDENS_CGS)
//					/ PARAM::yr;
//			hydro.eng = energy / PARAM::SEng_per_Mass;
//			if (hydro.eng < 0.0) {
//				std::cout << hydro.eng << "in integral.cpp " << std::endl;
//
//			}
//			hydro.abundances[0] = abundance_e;
//			hydro.abundances[1] = abundance_HI;
//			hydro.abundances[2] = abundance_HeI;
//			hydro.abundances[3] = abundance_OI;
//			hydro.abundances[4] = abundance_CI;
//			hydro.abundances[5] = abundance_H2;
//			hydro.abundances[6] = abundance_CO;
//			hydro.abundances[7] = abundance_HII;
//			hydro.abundances[8] = abundance_CII;
//			hydro.abundances[9] = abundance_FeII;
//			hydro.abundances[10] = abundance_SiII;
//			break;
//		}
		if (fabs(T - temprature_init) < 1.e-14) {
			//	std::cout << time_passed << " " << oneDynTimeStep * PARAM::ST
			//			<< " " << energy << std::endl;
			hydro.Lambda = Lambda * 1e26;
			hydro.Gamma = Gamma * 1e26;
			hydro.cooling_timescale = ((energy / fabs(Lambda)) * MASSDENS_CGS)
					/ PARAM::yr;
			hydro.eng = energy / PARAM::SEng_per_Mass;
			hydro.abundances[0] = abundance_e;
			hydro.abundances[1] = abundance_HI;
			hydro.abundances[2] = abundance_HeI;
			hydro.abundances[3] = abundance_OI;
			hydro.abundances[4] = abundance_CI;
			hydro.abundances[5] = abundance_H2;
			hydro.abundances[6] = abundance_CO;
			hydro.abundances[7] = abundance_HII;
			hydro.abundances[8] = abundance_CII;
			hydro.abundances[9] = abundance_FeII;
			hydro.abundances[10] = abundance_SiII;
			break;

		}
		abundance_e_init = abundance_e;
		abundance_HI_init = abundance_HI;
		abundance_OI_init = abundance_OI;
		abundance_CI_init = abundance_CI;
		abundance_H2_init = abundance_H2;
		abundance_CO_init = abundance_CO;
		abundance_HII_init = abundance_HII;
		abundance_CII_init = abundance_CII;
		abundance_FeII_init = abundance_FeII;
		abundance_SiII_init = abundance_SiII;
		temprature_init = T;
		energy_init = energy;
		dt *= .5;

	}

}

