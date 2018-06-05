#include "header.h"

void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time,
		boundary *box) {
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	/////////
	const PS::F64 Mass = 2.0;
	const PS::F64 Grav = 1.0;
	*end_time = 0.77;

	const PS::F64 dx = 4e-3;
	box->x = 1.4;
	box->y = 7.0 * sqrt(3) * dx;
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
				ith.eng = .5 * ith.vel * ith.vel
						+ ith.pres / ((5. / 3. - 1.0) * ith.dens)
						+ .5 * ith.MagneticB * ith.MagneticB / ith.dens;
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
				ith.eng = .5 * ith.vel * ith.vel
						+ ith.pres / ((5. / 3. - 1.0) * ith.dens)
						+ .5 * ith.MagneticB * ith.MagneticB / ith.dens;
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

void SetupICSphericalBlastWaves(PS::ParticleSystem<RealPtcl>& sph_system,
		PS::F64* end_time, boundary *box) {
	/////////
	//place ptcls
	/////////
	std::vector<RealPtcl> ptcl;
	/////////
	const PS::F64 Mass = 1.0;
	const PS::F64 Grav = 1.0;
	*end_time = 0.77;

	const PS::F64 dx = 1e-2;
	box->x = 0.5;
	box->y = 0.5;
	box->z = dx;

	PS::S32 id = 0;

	int rowL = -1;
	for (PS::F64 y = -box->y; y <= box->y; y += .5 * sqrt(3) * dx) {
		rowL += 1;
		for (PS::F64 x = -box->x; x < box->x; x += dx) {

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
				ith.pres = 10.0;
				ith.vel.y = 0;
				ith.vel.z = 0;

				ith.MagneticB.x = sqrt(2. * 3.141592);
				ith.MagneticB.y = sqrt(2. * 3.141592);
				ith.MagneticB.z = 0.0;
				ith.BoverDens = ith.MagneticB / ith.dens;
				ith.eng = .5 * ith.vel * ith.vel
						+ ith.pres / ((5. / 3. - 1.0) * ith.dens)
						+ .5 * ith.MagneticB * ith.MagneticB / ith.dens;

				ith.id = id++;
				ptcl.push_back(ith);

			}
		}
	}

	for (PS::U32 i = 0; i < ptcl.size(); ++i) {
		RealPtcl ith = ptcl[i];
		double x, y, z;
		x = ith.pos.x;
		y = ith.pos.y;
		z = ith.pos.z;
		const PS::F64 r = sqrt(x * x + y * y + z * z);

		if (r < 0.1) {

			ptcl[i].pres = 10.0;
			ptcl[i].dens = 1.0;
			ptcl[i].MagneticB.z = 0.0;
			ptcl[i].eng = .5 * ptcl[i].vel * ptcl[i].vel
					+ ptcl[i].pres / ((5. / 3. - 1.0) * ith.dens)
					+ .5 * ptcl[i].MagneticB * ptcl[i].MagneticB / ptcl[i].dens;

		} else {

			ptcl[i].pres = 0.1;
			ptcl[i].dens = 1.0;
			ptcl[i].MagneticB.z = 0.0;
			ptcl[i].eng = .5 * ptcl[i].vel * ptcl[i].vel
					+ ptcl[i].pres / ((5. / 3. - 1.0) * ith.dens)
					+ .5 * ptcl[i].MagneticB * ptcl[i].MagneticB / ptcl[i].dens;
		}

//		std::cout << "setup..." << ptcl[i].eng-.5 * ptcl[i].vel * ptcl[i].vel- .5 * ptcl[i].MagneticB * ptcl[i].MagneticB / ptcl[i].dens<< std::endl;


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

		sph_system[i].smth = PARAM::SMTH
				* pow(sph_system[i].mass / sph_system[i].dens,
						1.0 / (PS::F64) (PARAM::Dim));

		sph_system[i].setPressure();
	}
}

//#include "header.h"
//
//void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time, boundary *box) {
//	/////////
//	//place ptcls
//	/////////
//	std::vector<RealPtcl> ptcl;
//	/////////
//	const PS::F64 Mass = 2.0;
//	const PS::F64 Grav = 1.0;
//	*end_time = 0.77;
//
//	const PS::F64 dx =4e-3;
//	box->x = 1.0;
//	box->y = 0.1;
//	box->z = dx;
//	PS::S32 id = 0;
//
//	for (PS::F64 x = -box->x; x < 0.0; x += dx) {
//		for (PS::F64 y = -box->y ; y <= box->y; y += dx) {
//			for (PS::F64 z = 0; z < box->z; z += dx) {
//				RealPtcl ith;
//				ith.pos.x = x;
//				ith.pos.y = y;
//				ith.pos.z = z;
//
//				ith.dens = 1.0;
//				ith.mass = dx * dx * dx;
//				ith.pres = 0.95;
//				ith.vel.x = 1.2;
//				ith.vel.y = 0.01;
//				ith.vel.z = 0.5;
//
//				ith.MagneticB.x =2. / sqrt(4. * 3.141592);
//				ith.MagneticB.y = 3.6 / sqrt(4. * 3.141592);
//				ith.MagneticB.z = 2. / sqrt(4. * 3.141592);
//				ith.BoverDens = ith.MagneticB/ith.dens;
//				ith.eng = .5 * ith.vel * ith.vel + ith.pres / ((5. / 3. - 1.0) * ith.dens) + .5 * ith.MagneticB * ith.MagneticB / ith.dens;
//				ith.id = id++;
//				ptcl.push_back(ith);
//
//			}
//		}
//
//	}
//	for (PS::F64 x = 0.0; x <= box->x; x += dx) {
//		for (PS::F64 y = -box->y; y <= box->y; y += dx) {
//			for (PS::F64 z = 0; z < box->z; z += dx) {
//				RealPtcl ith;
//				ith.pos.x = x;
//				ith.pos.y = y;
//				ith.pos.z = z;
//
//				ith.dens = 1.0;
//				ith.mass = dx * dx * dx;
//				ith.pres = 1.0;
//				ith.vel.x = 0;
//				ith.vel.y = 0;
//				ith.vel.z = 0;
//
//				ith.MagneticB.x = 2. / sqrt(4. * 3.141592);
//				ith.MagneticB.y = 4. / sqrt(4. * 3.141592);
//				ith.MagneticB.z = 2. / sqrt(4. * 3.141592);
//				ith.BoverDens = ith.MagneticB/ith.dens;
//				ith.eng = .5 * ith.vel * ith.vel + ith.pres / ((5. / 3. - 1.0) * ith.dens) + .5 * ith.MagneticB * ith.MagneticB / ith.dens;
//				ith.id = id++;
//				ptcl.push_back(ith);
//
//			}
//		}
//
//	}
//	for (PS::U32 i = 0; i < ptcl.size(); ++i) {
//		ptcl[i].mass = 4.*box->x*box->y*box->z / (PS::F64) (ptcl.size());
//	}
//
//	if (PS::Comm::getRank() == 0) {
//		const PS::S32 numPtclLocal = ptcl.size();
//		sph_system.setNumberOfParticleLocal(numPtclLocal);
//		for (PS::U32 i = 0; i < ptcl.size(); ++i) {
//			sph_system[i] = ptcl[i];
//		}
//	} else {
//		sph_system.setNumberOfParticleLocal(0);
//	}
//	//Fin.
//	std::cout << "setup..." << ptcl.size() << std::endl;
//}
//
//void Initialize(PS::ParticleSystem<RealPtcl>& sph_system) {
//	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
//
//		sph_system[i].smth = PARAM::SMTH * pow(sph_system[i].mass / sph_system[i].dens, 1.0 / (PS::F64) (PARAM::Dim));
//
//		sph_system[i].setPressure();
//	}
//}

