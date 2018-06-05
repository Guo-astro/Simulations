#include "header.h"

void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time, boundary *box) {
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	/////////
	const PS::F64 Mass = 2.0;
	const PS::F64 Grav = 1.0;
	*end_time = 0.77;

	const PS::F64 dx = 4e-3;
	box->x = 1.0;
	box->y = 7.0 * .5 * sqrt(3) * dx;
	box->z = dx;
	PS::S32 id = 0;
	int rowL = -1;
	for (PS::F64 y = -box->y; y <= box->y; y += .5 * sqrt(3) * dx) {
		rowL += 1;
		for (PS::F64 x = -box->x; x < 0.0; x += dx) {

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
				ith.mass = dx * dx * dx;
				ith.pres = 20.0;
				ith.vel.x = 10.0;
				ith.vel.y = 0;
				ith.vel.z = 0;

				ith.MagneticB.x = 5. / sqrt(4. * 3.141592);
				ith.MagneticB.y = 5. / sqrt(4. * 3.141592);
				ith.MagneticB.z = 0.0;
				ith.BoverDens = ith.MagneticB / ith.dens;
				ith.eng = .5 * ith.vel * ith.vel + ith.pres / ((5. / 3. - 1.0) * ith.dens) + .5 * ith.MagneticB * ith.MagneticB / ith.dens;
				ith.id = id++;
				ptcl.push_back(ith);

			}
		}

	}
	int rowR = -1;
	for (PS::F64 y = -box->y; y <= box->y; y += .5 * sqrt(3) * dx) {
		rowR += 1;
		for (PS::F64 x = 0.0; x <= box->x; x += dx) {

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
				ith.vel.x = -10.0;
				ith.vel.y = 0;
				ith.vel.z = 0;

				ith.MagneticB.x = 5. / sqrt(4. * 3.141592);
				ith.MagneticB.y = 5. / sqrt(4. * 3.141592);
				ith.MagneticB.z = 0.0;
				ith.BoverDens = ith.MagneticB / ith.dens;
				ith.eng = .5 * ith.vel * ith.vel + ith.pres / ((5. / 3. - 1.0) * ith.dens) + .5 * ith.MagneticB * ith.MagneticB / ith.dens;
				ith.id = id++;
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
	//Fin.
	std::cout << "setup..." << ptcl.size() << std::endl;
}

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0 / (PS::F64) (PARAM::Dim));

		sph_system[i].setPressure();
	}
}

void SetupIC_HVCC_Chem_file(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time, boundary *box, PS::F64 imp) {
	/////////
	//place ptcls
	/////////
	//place ptcls
	/////////
	FileHeader header;
#ifdef RESTART
	sph_system.readParticleBinary<FileHeader>(
			"../../POLYTROPE_TESTS/INIT_CONDS/HVCC_chem_3450.dat",

			header);
#else
//	sph_system.readParticleBinary < FileHeader > ("../../POLYTROPE_TESTS/INIT_CONDS/FORCE1e5.dat",header);
	sph_system.readParticleBinary<FileHeader>("IC/FORCE_RELAXED_1e5.dat", header);

#endif

	/////////
	*end_time = 1.5 * 1e6 * PARAM::yr / PARAM::ST;
	//Fin.
	std::cout << "setup..." << std::endl;

}

