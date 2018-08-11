#include "header.h"
void Initialize(PS::ParticleSystem<RealPtcl>& sph_system) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		sph_system[i].smth = PARAM::SMTH
				* pow(sph_system[i].mass / sph_system[i].dens,
						1.0 / (PS::F64) (PARAM::Dim));

		sph_system[i].setPressure();
	}
}
double f1(double t, double x, double v) {
	return v;
}

double f2(double t, double x, double v) {
	return -pow(x, 1.5) - 2.0 / t * v;
}

double getPhi(double x) {

	double a = 1.000744397;

	double b = 1.50738116e-2;

	double c = -2.334693e-1;

	double d = 7.8109735e-2;

	double e = -1.021152772e-2;

	double f = 4.833905487e-4;

	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3)
			+ e * pow(x, 4) + f * pow(x, 5);

	return phi;

}

double getPhi_dash(double x) {

	double a = 8.841062613e-3;

	double b = -3.86087802e-1;

	double c = 6.92279238e-2;

	double d = 6.159126e-3;

	double e = -2.44979111e-3;

	double f = 1.6204056e-4;

	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3)
			+ e * pow(x, 4) + f * pow(x, 5);

	return phi;

}

//double getPhi(double x) {
//	//		double a = 9.999995857999977E-01;
//	//		double b = -3.215355146247138E-04;
//	//		double c = -1.667864474453594E-01;
//	//		double d = 5.252845006191754E-04;
//	//		double e = 1.155630978128144E-02;
//	//		double f = 1.050795977407310E-03;
//	//		double g = -1.389661348847463E-03;
//	//		double h = 2.648078982351008E-04;
//	//		double i = -1.692465157938468E-05;
////		double i = 5.174e-16;
////		double h = - 4.433e-13;
////		double g = + 1.574e-10;
////		double f =- 2.992e-08;
////		double e = 3.286e-06;
////		double d =- 0.000209;
////		double c = 0.007358;
////		double b = - 0.1254;
////		double a =+ 0.7485;
////		double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7) + i * pow(x, 8);
//	double phi = pow(1.0 + pow(x / 2.88, 2.0), -1.47);
//	return phi;
//
//}
//double getPhi_dash(double x) {
//	//		double a = -3.257098683322130E-04;
//	//		double b = -3.335390631705378E-01;
//	//		double c = 1.494801635423298E-03;
//	//		double d = 4.629812975807299E-02;
//	//		double e = 5.242083986884144E-03;
//	//		double f = -8.358773124387748E-03;
//	//		double g = 1.867615193618750E-03;
//	//		double h = -1.388085203319444E-04;
//	//		double i = 2.998565502661576E-07;
////		double i = 1.193e-16;
////		double h = -9.659e-14;
////		double g = +3.172e-11;
////		double f = -5.37e-09;
////		double e = +4.87e-07;
////		double d = -2.107e-05;
////		double c = +0.0001497;
////		double b = +0.01751;
////		double a = -0.4394;
////		double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7) + i * pow(x, 8);
//
//	double a1 = 1.81178953;
//	double b1 = 3.2662366;
//	double c1 = 1.11317173;
//	double phi = a1 * pow(1 + pow(x / b1, 2), -c1) * 2 * x / (b1 * b1);
//
//	return phi;
//
//}
//
//double getMass(double x) {
//	//		double a = -3.257098683322130E-04;
//	//		double b = -3.335390631705378E-01;
//	//		double c = 1.494801635423298E-03;
//	//		double d = 4.629812975807299E-02;
//	//		double e = 5.242083986884144E-03;
//	//		double f = -8.358773124387748E-03;
//	//		double g = 1.867615193618750E-03;
//	//		double h = -1.388085203319444E-04;
//	//		double i = 2.998565502661576E-07;
//	double i = -5.353e-05;
//	double h = 0.002491;
//	double g = -0.04842;
//	double f = +0.5035;
//	double e = -2.924;
//	double d = +8.333;
//	double c = -2.768;
//	double b = +0.7079;
//	double a = -0.03603;
//	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7) + i * pow(x, 8);
//
//	return phi;
//
//}
double fRand(double fMin, double fMax) {
	double f = (double) rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
void SetupIC_Polytrope_Dens_Relaxation(PS::ParticleSystem<RealPtcl>& sph_system,
		PS::F64* end_time) {
	std::vector<RealPtcl> ptcl;

	const PS::F64 Radi = 1.5;
	*end_time = 10000;
	const PS::F64 dx = 1 / 8;
	PS::S32 id = 0;
	int ptcl_num = 200000;
//	int shell_num = 8000;

	double M_tot = -4.0 * M_PI * getPhi_dash(Radi) * Radi * Radi;
//	double M_tot = getMass(Radi);
//integrate [0,3] 4pi*x^2*(1+(x/2.88)^2)^-1.47
	std::cout << "Mtot: " << M_tot << std::endl;

//	double shell_width = PARAM::poly_r / shell_num;
//	double mass = PARAM::poly_M / ptcl_num;
	double mass = M_tot / ptcl_num;
	for (int i = 0; i < ptcl_num; i++) {
		double costheta = fRand(-1, 1);

		double Radi_ = Radi * fRand(0, 1);
		double phi = fRand(0, 2.0 * M_PI);
		double x = Radi_ * sqrt(1. - costheta * costheta) * cos(phi);
		double y = Radi_ * sqrt(1. - costheta * costheta) * sin(phi);
		double z = Radi_ * costheta;
		double r = sqrt(x * x + y * y + z * z);
//
		if (r >= Radi - 0.5 && r < 0.0001) {
			continue;
		}
		RealPtcl ith;
		ith.pos.x = x;
		ith.pos.y = y;
		ith.pos.z = z;
		ith.dens = pow(getPhi(r), 15.45);
//		std::cout << "Ptcls num: " << r << std::endl;

		ith.mass = mass;
		ith.pres = ith.dens;
		ith.eng = 2.5;
		ith.id = id++;
		ptcl.push_back(ith);
	}
//	for (int i = 0; i < shell_num - 1; i++) {
//
//		const PS::F64 _r = shell_width * i;
//		const PS::F64 _r_next = shell_width * (i + 1);
//
//		double current_phi = getPhi(_r);
////		double current_phi_dash = getPhi_dash(_r);
////
////		double current_phi_next = getPhi(_r_next);
////		double current_phi_dash_next = getPhi_dash(_r_next);
////
////		double inner_mass = -4.0 * M_PI * _r * _r * current_phi_dash;
////		double inner_mass_next = -4.0 * M_PI * _r_next * _r_next * current_phi_dash_next;
//
//		double inner_mass = getMass(_r);
//		double inner_mass_next =  getMass(_r_next);
//		int ptcl_num_inshell = (inner_mass_next - inner_mass) / mass;
//
//		for (int i = 0; i < ptcl_num_inshell; i++) {
//			double costheta = fRand(-1, 1);
//			double phi = fRand(0, 2.0 * M_PI);
//			double x = _r_next * sqrt(1. - costheta * costheta) * cos(phi);
//			double y = _r_next * sqrt(1. - costheta * costheta) * sin(phi);
//			double z = _r_next * costheta;
//			if (sqrt(x * x + y * y + z * z) <= _r || sqrt(x * x + y * y + z * z) >= Radi-0.2) {
//				continue;
//			}
//			RealPtcl ith;
//			ith.pos.x = x;
//			ith.pos.y = y;
//			ith.pos.z = z;
//			ith.dens = current_phi;
////			std::cout << "Ptcls num: " << current_phi << std::endl;
//
//			ith.mass = mass;
//			ith.pres = ith.dens;
//			ith.eng = 2.5;
//			ith.id = id++;
//			ptcl.push_back(ith);
//		}
//
//	}

	for (PS::U32 i = 0; i < ptcl.size(); ++i) {
		ptcl[i].smth = pow(ptcl[i].mass / ptcl[i].dens, 1 / 3);
	}

	if (PS::Comm::getRank() == 0) {
		const PS::S32 numPtclLocal = ptcl.size();
		sph_system.setNumberOfParticleLocal(numPtclLocal);
		for (PS::U32 i = 0; i < ptcl.size(); ++i) {
			sph_system[i] = ptcl[i];
		}

		std::cout << "Ptcls num: " << numPtclLocal << std::endl;
	} else {
		sph_system.setNumberOfParticleLocal(0);
	}
}

void SetupIC_BonnerEbert_Force_Relaxation(PS::ParticleSystem<RealPtcl>& sph_system,
		PS::F64 *end_time) {
/////////
//place ptcls
/////////
	*end_time = 10000000;

	FileHeader header;
#ifdef RESTART
	sph_system.readParticleAscii<FileHeader>(
			"CO_0.2_0.4_mass/ascii/HVCC_imp_7.00/HVCC_0313.txt", header);

#else
	sph_system.readParticleAscii("result/iso8000.dat", header);
	PS::U32 ptcl_size = sph_system.getNumberOfParticleGlobal();
	std::cout << "setup...ptcl num " << std::endl;

	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		double mu = 2.34;

		double mu_new = 2.34;

		const double T = 30;

		const double n_H = 10000;

		const double x_e = 1e-5;

		const double n_H2 = 4000;

		const double x_2 = 2. * n_H2 / n_H;

		double NUMDENS_CGS = 0;

		for (int j = 0; j < 100; j++) {

			mu = mu_new;


			NUMDENS_CGS = sph_system[i].dens * PARAM::SMassDens
					/ (mu * PARAM::PROTONMASS);

			sph_system[i].abundances[0] = fitAB_E(NUMDENS_CGS);

			sph_system[i].abundances[1] = fitAB_HI(NUMDENS_CGS);

			sph_system[i].abundances[2] = PARAM::i_abundance_HeI;

			sph_system[i].abundances[3] = PARAM::i_abundance_OI;

			sph_system[i].abundances[4] = PARAM::i_abundance_CI;

			sph_system[i].abundances[5] = fitAB_H2(NUMDENS_CGS);

			sph_system[i].abundances[6] = fitAB_CO(NUMDENS_CGS);

			sph_system[i].abundances[7] = fmax(
					1. - 2. * sph_system[i].abundances[5]
							- sph_system[i].abundances[1], 0);

			sph_system[i].abundances[8] = fmax(
					sph_system[i].abundances[0] - sph_system[i].abundances[7]
							- PARAM::i_abundance_Si, 0);

			sph_system[i].abundances[9] = PARAM::i_abundance_Fe;

			sph_system[i].abundances[10] = PARAM::i_abundance_Si;

			mu_new = sph_system[i].abundances[1] + sph_system[i].abundances[7]
					+ 2.0 * sph_system[i].abundances[5]
					+ 4.0 * sph_system[i].abundances[2];

			mu_new /= (sph_system[i].abundances[1] + sph_system[i].abundances[7]
					+ sph_system[i].abundances[5] + sph_system[i].abundances[2]
					+ sph_system[i].abundances[0]);

			if (fabs(mu - mu_new) < 1e-3) {

				break;

			}

		}
		if(fabs(sph_system[i].dens -1.0)<2e-2){
		std::cout << "mu is: "
					<< mu << std::endl;

		std::cout << "we want PARAM::SMassDens to be : "
				<< mu * n_H * PARAM::PROTONMASS_CGS << std::endl;


		std::cout << " Therefore the SM: "
				<< PARAM::SL * PARAM::SL * PARAM::SL * mu * n_H
						* PARAM::PROTONMASS_CGS / PARAM::M_SUN_cgs << std::endl;


		sph_system[i].NUMDENS = NUMDENS_CGS;

		sph_system[i].pres = (1.1 + sph_system[i].abundances[0]
				- sph_system[i].abundances[5]) * NUMDENS_CGS * T
				* PARAM::KBOLTZ_CGS / PARAM::SPres;

		sph_system[i].eng = sph_system[i].pres
				/ ((PARAM::GAMMA - 1.0) * sph_system[i].dens);

		double tem = mu * PARAM::PROTONMASS_CGS * sph_system[i].pres
				* PARAM::SEng_per_Mass
				/ (sph_system[i].dens * PARAM::KBOLTZ_cgs
						* (1.1 + PARAM::i_abundance_e - PARAM::i_abundance_H2));
		sph_system[i].temp = tem;
		sph_system[i].mu = mu;

		sph_system[i].smth = pow(sph_system[i].mass / sph_system[i].dens,
				1.0 / (PS::F64) (PARAM::Dim));
		std::cout << "dens "
				<<NUMDENS_CGS<< "K "
				<< (PARAM::GAMMA -1)*sph_system[i].eng <<" Temp "<<sph_system[i].temp<< std::endl;
		}

	}

#endif

}
