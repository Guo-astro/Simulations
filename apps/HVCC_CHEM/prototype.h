#pragma once
struct boundary{
   PS::F64 x, y, z;
};
void DisplayInfo(void);
void SetupIC(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64*, boundary *box);
void SetupIC_HVCC_Chem_file(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64 *end_time, boundary *box, PS::F64 imp);
void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system);

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system);

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system);
