#include "header.h"

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system) {
	PS::F64 dt = 1.0e+30; //set VERY LARGE VALUE
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
		dt = std::min(dt, sph_system[i].dt);
	}
	return PS::Comm::getMinValue(dt);
}

void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
		sph_system[i].vel_half = sph_system[i].vel + 0.5 * dt * (sph_system[i].acc + sph_system[i].grav + sph_system[i].extF);
		sph_system[i].eng_half = sph_system[i].eng + 0.5 * dt * sph_system[i].eng_dot;
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
		sph_system[i].vel += dt * (sph_system[i].acc + sph_system[i].grav + sph_system[i].extF);
		sph_system[i].eng += dt * sph_system[i].eng_dot;
		if (sph_system[i].eng < 0.0) {
			sph_system[i].eng = 1e-4;
		}
//		sph_system[i].BoverDens += dt * sph_system[i].BoverDens_dot;
//		sph_system[i].MagneticB = sph_system[i].BoverDens * sph_system[i].dens;
//		sph_system[i].MagneticB.z = 0.0;

	}
}

void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		sph_system[i].vel = sph_system[i].vel_half + 0.5 * dt * (sph_system[i].acc + sph_system[i].grav + sph_system[i].extF);
		sph_system[i].eng = sph_system[i].eng_half + 0.5 * dt * sph_system[i].eng_dot;
		if (sph_system[i].eng < 0.0) {
			sph_system[i].eng = 1e-4;
		}
//		sph_system[i].BoverDens = sph_system[i].BoverDens_half + 0.5 * dt * sph_system[i].BoverDens_dot;
//		sph_system[i].MagneticB = sph_system[i].BoverDens * sph_system[i].dens;
//		sph_system[i].MagneticB.z = 0.0;
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
