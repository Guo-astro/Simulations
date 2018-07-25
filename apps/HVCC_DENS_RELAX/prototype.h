#pragma once
struct boundary{
   PS::F64 x, y, z;
};
void DisplayInfo(void);
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64*, boundary *box);
void SetupIC_Polytrope_Dens_Relaxation(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time);
void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system);

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system);

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system);
