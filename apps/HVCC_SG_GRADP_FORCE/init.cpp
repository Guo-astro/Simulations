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
	double a = 9.999995857999977E-01;
	double b = -3.215355146247138E-04;
	double c = -1.667864474453594E-01;
	double d = 5.252845006191754E-04;
	double e = 1.155630978128144E-02;
	double f = 1.050795977407310E-03;
	double g = -1.389661348847463E-03;
	double h = 2.648078982351008E-04;
	double i = -1.692465157938468E-05;
	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7) + i * pow(x, 8);
	return phi;

}
double getPhi_dash(double x) {
	double a = -3.257098683322130E-04;
	double b = -3.335390631705378E-01;
	double c = 1.494801635423298E-03;
	double d = 4.629812975807299E-02;
	double e = 5.242083986884144E-03;
	double f = -8.358773124387748E-03;
	double g = 1.867615193618750E-03;
	double h = -1.388085203319444E-04;
	double i = 2.998565502661576E-07;
	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7) + i * pow(x, 8);
	return phi;

}
double fRand(double fMin, double fMax) {
	double f = (double) rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
void SetupIC_Polytrope_Dens_Relaxation(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time) {
	std::vector<RealPtcl> ptcl;

	const PS::F64 Radi = PARAM::xi - 0.005;
	*end_time = 10000;
	const PS::F64 dx = 1 / 8;
	PS::S32 id = 0;
	int ptcl_num = 10000;
	int shell_num = 5000;
	double shell_width = PARAM::poly_r / shell_num;
	double mass = PARAM::poly_M / ptcl_num;
	for (int i = 0; i < shell_num - 1; i++) {

		const PS::F64 _r = shell_width * i;
		const PS::F64 _r_next = shell_width * (i + 1);

		double current_phi = getPhi(_r);
		double current_phi_dash = getPhi_dash(_r);

		double current_phi_next = getPhi(_r_next);
		double current_phi_dash_next = getPhi_dash(_r_next);

		double inner_mass = -4.0 * M_PI * _r * _r * current_phi_dash;
		double inner_mass_next = -4.0 * M_PI * _r_next * _r_next * current_phi_dash_next;
		int ptcl_num_inshell = (inner_mass_next - inner_mass) / mass;

		for (int i = 0; i < ptcl_num_inshell; i++) {
			double costheta = fRand(-1, 1);
			double phi = fRand(0, 2.0 * M_PI);
			double x = _r_next * sqrt(1. - costheta * costheta) * cos(phi);
			double y = _r_next * sqrt(1. - costheta * costheta) * sin(phi);
			double z = _r_next * costheta;
			if (sqrt(x * x + y * y + z * z) <= _r) {
				continue;
			}
			RealPtcl ith;
			ith.pos.x = x;
			ith.pos.y = y;
			ith.pos.z = z;
			ith.dens = pow(current_phi, 1.5);
			ith.mass = mass;
			ith.pres = 0.5 * pow(ith.dens, 2.0);
			ith.eng = 2.5;
			ith.id = id++;
			ptcl.push_back(ith);
		}

	}

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

void SetupIC_Restart_Polytrope_Dens_Relaxation(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time) {
	*end_time = 1000000 * 0.74 * 1e6 * PARAM::yr / PARAM::ST;
	FileHeader header;
#ifdef RESTART
	sph_system.readParticleAscii<FileHeader>(
			"CO_0.2_0.4_mass/ascii/HVCC_imp_7.00/HVCC_0313.txt", header);

#else
	sph_system.readParticleAscii("result/0800.dat", header);
	PS::U32 ptcl_size = sph_system.getNumberOfParticleGlobal();
	std::cout << "setup...ptcl num " << ptcl_size << std::endl;
	//Reset Energy
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

//		sph_system[i].eng = sph_system[i].pres / ((PARAM::GAMMA - 1) * sph_system[i].dens);
//		sph_system[i].smth = pow(sph_system[i].mass / sph_system[i].dens, 1.0 / (PS::F64) (PARAM::Dim));
//		std::cout << "dens " << pow(sph_system[i].mass / sph_system[i].dens, 1. / 3.) << " " << sph_system[i].smth << std::endl;

	}

#endif
}
