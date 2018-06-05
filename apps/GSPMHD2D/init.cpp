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

	const PS::F64 dx =4.e-3;
	box->x = 1.0;
	box->y = 0.1;
	PS::S32 id = 0;

	for (PS::F64 x = -box->x; x < 0.0; x += dx) {
		for (PS::F64 y = -box->y ; y <= box->y; y += dx) {
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;

				ith.dens = 1.0;
				ith.mass = dx * dx * dx;
				ith.pres = 20.0;
				ith.vel.x = 10.0;
				ith.vel.y = 0;
				ith.MagneticB.x = 5. / sqrt(4. * 3.141592);
				ith.MagneticB.y = 5. / sqrt(4. * 3.141592);
				ith.BoverDens = ith.MagneticB/ith.dens;
				ith.eng = .5 * ith.vel * ith.vel + ith.pres / ((5. / 3. - 1.0) * ith.dens) + .5 * ith.MagneticB * ith.MagneticB / ith.dens;
				ith.id = id++;
				ptcl.push_back(ith);

		}

	}
	for (PS::F64 x = 0.0; x <= box->x; x += dx) {
		for (PS::F64 y = -box->y; y <= box->y; y += dx) {
				RealPtcl ith;
				ith.pos.x = x;
				ith.pos.y = y;

				ith.dens = 1.0;
				ith.mass = dx * dx * dx;
				ith.pres = 1.0;
				ith.vel.x = -10.0;
				ith.vel.y = 0;

				ith.MagneticB.x = 5. / sqrt(4. * 3.141592);
				ith.MagneticB.y = 5. / sqrt(4. * 3.141592);
				ith.BoverDens = ith.MagneticB/ith.dens;
				ith.eng = .5 * ith.vel * ith.vel + ith.pres / ((5. / 3. - 1.0) * ith.dens) + .5 * ith.MagneticB * ith.MagneticB / ith.dens;
				ith.id = id++;
				ptcl.push_back(ith);

		}

	}
	for (PS::U32 i = 0; i < ptcl.size(); ++i) {
		ptcl[i].mass = 4.*box->x*box->y / (PS::F64) (ptcl.size());
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
