#include "header.h"

void SetupIC_BLASTWAVE(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time, boundary *box) {
	// Place SPH particles
	*end_time = 0.77;
	std::vector<RealPtcl> ptcl;

	const PS::F64 dx = 8e-3;
	box->x = 1.0;
	box->y = 11 * .5 * sqrt(3) * dx;
	box->z = dx;
	PS::S32 id = 0;

	int rowL = -1;
	for (PS::F64 y = -box->y + .5 * sqrt(3) * dx; y <= box->y + .5 * sqrt(3) * dx; y += .5 * sqrt(3) * dx) {
		rowL += 1;
		for (PS::F64 x = -box->x + .5 * dx; x < .5 * dx; x += dx) {

			for (PS::F64 z = 0; z < box->z; z += dx) {

				RealPtcl ith;
				if (rowL % 2 == 1) {
					ith.pos.x = x + .5 * dx;
				} else {
					ith.pos.x = x;
				}
				ith.pos.y = y;
				ith.pos.z = z;

				ith.dens = 1.0;
				ith.mass = ith.dens * dx * dx * dx;
				ith.pres = 3000;
				ith.eng = ith.pres / ((PARAM::GAMMA - 1.0) * ith.dens);
				ith.id = id++;
				ith.smth = PARAM::SMTH * pow(ith.mass / ith.dens, 1.0 / (PS::F64) (PARAM::Dim));
				ptcl.push_back(ith);

			}
		}

	}
	int rowR = -1;
	for (PS::F64 y = -box->y + .5 * sqrt(3) * dx; y <= box->y + .5 * sqrt(3) * dx; y += .5 * sqrt(3) * dx) {
		rowR += 1;
		for (PS::F64 x = .5 * dx; x <= box->x + .5 * dx; x += dx) {

			for (PS::F64 z = 0; z < box->z; z += dx) {
				RealPtcl ith;
				if (rowR % 2 == 1) {
					ith.pos.x = x + .5 * dx;
				} else {
					ith.pos.x = x;
				}
				ith.pos.y = y;
				ith.pos.z = z;

				ith.dens = 1.0;
				ith.mass = dx * dx * dx;
				ith.pres = 1.0;
				ith.dens = 1.0;
				ith.mass = ith.dens * dx * dx * dx;
				ith.pres = 1.e-7;
				ith.eng = ith.pres / ((PARAM::GAMMA - 1.0) * ith.dens);
				ith.id = id++;
				ith.smth = PARAM::SMTH * pow(ith.mass / ith.dens, 1.0 / (PS::F64) (PARAM::Dim));
				ptcl.push_back(ith);

			}
		}

	}
	for (PS::U32 i = 0; i < ptcl.size(); ++i) {
		ptcl[i].mass = 4. * box->x * box->y * box->z / (PS::F64) (ptcl.size());
	}

	if (PS::Comm::getRank() == 0) {
		const PS::S32 numPtclLocal = ptcl.size();
		sph_system.setNumberOfParticleLocal(numPtclLocal);
		for (PS::U32 i = 0; i < ptcl.size(); ++i) {
			sph_system[i] = ptcl[i];
		}
	} else {
		sph_system.setNumberOfParticleLocal(0);
	}
}

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0 / (PS::F64) (PARAM::Dim));

		sph_system[i].setPressure();
	}
}
void SetupIC_BonnerEbert_Force_Relaxation(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time) {
	/////////
	//place ptcls
	/////////
	*end_time = 1000000;

	FileHeader header;
#ifdef RESTART
	sph_system.readParticleAscii<FileHeader>(
			"CO_0.2_0.4_mass/ascii/HVCC_imp_7.00/HVCC_0313.txt", header);

#else
//	sph_system.readParticleAscii("result/poly7000.dat", header);
	sph_system.readParticleAscii("result/1743.dat", header);
	PS::U32 ptcl_size = sph_system.getNumberOfParticleGlobal();
	std::cout << "setup...ptcl num " << ptcl_size << std::endl;
	//Reset Energy
//	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//
//		sph_system[i].eng = sph_system[i].pres / ((PARAM::GAMMA - 1) * sph_system[i].dens);
//		sph_system[i].smth = pow(sph_system[i].mass / sph_system[i].dens, 1.0 / (PS::F64) (PARAM::Dim));
//		std::cout << "dens " << sph_system[i].dens << " " << sph_system[i].eng << std::endl;
//
//	}

#endif

}

