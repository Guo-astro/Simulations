#pragma once
struct boundary{
   PS::F64 x, y, z;
};
void DisplayInfo(void);
void SetupICBlastWave(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64*, boundary *box);
void SetupICConvergentTest(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time,
		boundary *box);
void SetupICSlowShock(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time,
		boundary *box);
void SetupICOrszagTang(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time,
		boundary *box);
void SetupICBrioWu(PS::ParticleSystem<RealPtcl>& sph_system, PS::F64* end_time,
		boundary *box);
void SetupICCurrentSheet(PS::ParticleSystem<RealPtcl>& sph_system,
		PS::F64* end_time, boundary *box);
void SetupICFieldLoop(PS::ParticleSystem<RealPtcl>& sph_system,
		PS::F64* end_time, boundary *box);
void SetupICSphericalBlastWaves(PS::ParticleSystem<RealPtcl>& sph_system,
		PS::F64* end_time, boundary *box) ;
void InitialKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FullDrift(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void Predict(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);
void FinalKick(PS::ParticleSystem<RealPtcl>& sph_system, const PS::F64 dt);

void Initialize(PS::ParticleSystem<RealPtcl>& sph_system);

PS::F64 getTimeStepGlobal(const PS::ParticleSystem<RealPtcl>& sph_system);

void CheckConservativeVariables(const PS::ParticleSystem<RealPtcl>& sph_system);
