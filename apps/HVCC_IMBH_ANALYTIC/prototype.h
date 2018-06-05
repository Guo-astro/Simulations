#pragma once
struct boundary{
   PS::F64 x, y, z;
};
void DisplayInfo(void);
void SetupIC_BLASTWAVE(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time, boundary *box);
void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt, const PS::F64 time);
void SetupIC_HVCC_file(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time);
void SetupIC_HVCC_RESTART_file(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time);
void Initialize(PS::ParticleSystem<RealPtcl>& sph_system);

PS::F64 getTimeStepGlobal( PS::ParticleSystem<RealPtcl>& sph_system);

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system);
