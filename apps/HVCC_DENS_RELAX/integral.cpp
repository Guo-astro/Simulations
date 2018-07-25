#include "header.h"

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system) {
	double dt = 1e+30, maxsgrav = -1e+30, maxbhgrav = -1e+30, minsmth = 1e+30, maxsnds = -1e+30;
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		dt = std::min(dt, sph_system[i].dt);
	}

	return PS::Comm::getMinValue(dt);
}

void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt) {
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
		sph_system[i].pos += (sph_system[i].acc + sph_system[i].grav) * sph_system[i].dt * sph_system[i].dt * .5;
//		if(sqrt(sph_system[i].pos*sph_system[i].pos) >= PARAM::xi ){
//			sph_system[i].pos = 0.001;
//		}
//		if(sph_system[i].pos.x != sph_system[i].pos.x){
//			std::cout<<sph_system[i].acc<<std::endl;
//		}
	}
}

