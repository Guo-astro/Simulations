#include "header.h"
void Initialize(PS::ParticleSystem<RealPtcl>& sph_system) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0 / (PS::F64) (PARAM::Dim));

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

	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5);

 	return phi;



}

double getPhi_dash(double x) {

	double a = 8.841062613e-3;

	double b = -3.86087802e-1;

	double c = 6.92279238e-2;

	double d = 6.159126e-3;

	double e = -2.44979111e-3;

	double f = 1.6204056e-4;

	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5);

	return phi;



}

double getMass(double x) {
	//		double a = -3.257098683322130E-04;
	//		double b = -3.335390631705378E-01;
	//		double c = 1.494801635423298E-03;
	//		double d = 4.629812975807299E-02;
	//		double e = 5.242083986884144E-03;
	//		double f = -8.358773124387748E-03;
	//		double g = 1.867615193618750E-03;
	//		double h = -1.388085203319444E-04;
	//		double i = 2.998565502661576E-07;
	double i = -5.353e-05;
	double h = 0.002491;
	double g = -0.04842;
	double f = +0.5035;
	double e = -2.924;
	double d = +8.333;
	double c = -2.768;
	double b = +0.7079;
	double a = -0.03603;
	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7) + i * pow(x, 8);

	return phi;

}
double fRand(double fMin, double fMax) {
	double f = (double) rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
void SetupIC_Polytrope_Dens_Relaxation(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time) {
	std::vector<RealPtcl> ptcl;

	const PS::F64 Radi = 6.5;
	*end_time = 10000;
	const PS::F64 dx = 1 / 8;
	PS::S32 id = 0;
	int ptcl_num = 10000;
//	int shell_num = 8000;

//	double M_tot = -4.0 * M_PI * getPhi_dash(Radi) * Radi * Radi;
	double M_tot = getMass(Radi);
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
		if (r >= Radi - 1.0) {
			continue;
		}
		RealPtcl ith;
		ith.pos.x = x;
		ith.pos.y = y;
		ith.pos.z = z;
		ith.dens = getPhi(r);
//		std::cout << "Ptcls num: " << r << std::endl;

		ith.mass = mass;
		ith.pres = .5 * pow(ith.dens, 2);
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

void SetupIC_BonnerEbert_Force_Relaxation(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time) {
/////////
//place ptcls
/////////
	*end_time = 10000000;

	FileHeader header;
#ifdef RESTART
	sph_system.readParticleAscii<FileHeader>(
			"CO_0.2_0.4_mass/ascii/HVCC_imp_7.00/HVCC_0313.txt", header);

#else
	sph_system.readParticleAscii("result/20180304_Polytrope_7651.dat", header);
	PS::U32 ptcl_size = sph_system.getNumberOfParticleGlobal();
	std::cout << "setup...ptcl num " << std::endl;
//Reset Energy
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		sph_system[i].eng = sph_system[i].pres / ((PARAM::GAMMA - 1) * sph_system[i].dens);
		sph_system[i].smth = pow(sph_system[i].mass / sph_system[i].dens, 1.0 / (PS::F64) (PARAM::Dim));
		std::cout << "dens " << pow(sph_system[i].mass / sph_system[i].dens, 1. / 3.) << " " << sph_system[i].smth << std::endl;

	}

#endif

}
