#pragma once

class CalcDensity {
	kernel_t kernel;
public:
	void operator ()(const EPI::Dens* const ep_i, const PS::S32 Nip,
			const EPJ::Dens* const ep_j, const PS::S32 Njp,
			RESULT::Dens* const dens) {
		for (PS::S32 i = 0; i < Nip; ++i) {
			dens[i].clear();
			const EPI::Dens& ith = ep_i[i];
			for (PS::S32 j = 0; j < Njp; ++j) {
				const EPJ::Dens& jth = ep_j[j];

				const PS::F64vec dr = ith.pos - jth.pos;
				dens[i].dens += jth.mass * kernel.W(dr, ith.smth);

			}
			dens[i].smth = PARAM::SMTH
					* pow(ith.mass / dens[i].dens, 1.0 / (PS::F64)(PARAM::Dim));
		}
	}
};

class CalcDerivative {
	kernel_t kernel;
public:
	void operator ()(const EPI::Drvt* ep_i, const PS::S32 Nip,
			const EPJ::Drvt* ep_j, const PS::S32 Njp,
			RESULT::Drvt* const drvt) {
		for (PS::S32 i = 0; i < Nip; ++i) {
			drvt[i].clear();
			const EPI::Drvt& ith = ep_i[i];

			for (PS::S32 j = 0; j < Njp; ++j) {
				const EPJ::Drvt& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64vec dv = ith.vel - jth.vel;
				if (dr == 0.0) {
					continue;
				}
				const PS::F64vec ndir = dr / sqrt(dr * dr);
				const PS::F64vec MagneticBPerpj = jth.MagneticB
						- (jth.MagneticB * ndir) * ndir;
				const PS::F64vec MagneticBPerpi = ith.MagneticB
						- (ith.MagneticB * ndir) * ndir;
				const PS::F64 Bperpj2 = MagneticBPerpj * MagneticBPerpj;
				const PS::F64 Bperpi2 = MagneticBPerpi * MagneticBPerpi;
				const PS::F64 PTj = jth.pres
						+ .5 * MagneticBPerpj * MagneticBPerpj;
				const PS::F64 PTi = ith.pres
						+ .5 * MagneticBPerpi * MagneticBPerpi;
				const PS::F64vec vperpj = jth.vel - (jth.vel * ndir) * ndir;
				const PS::F64vec vperpi = ith.vel - (ith.vel * ndir) * ndir;
				drvt[i].grad_dens += jth.mass * kernel.gradW(dr, ith.smth);

				drvt[i].gradP += jth.mass * (jth.pres - ith.pres)
						* kernel.gradW(dr, ith.smth) / jth.dens;

				drvt[i].gradVpara += jth.mass * ((jth.vel - ith.vel) * ndir)
						* kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradvperp_x += jth.mass * (vperpj.x - vperpi.x)
						* kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradvperp_y += jth.mass * (vperpj.y - vperpi.y)
						* kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradvperp_z += jth.mass * (vperpj.z - vperpi.z)
						* kernel.gradW(dr, ith.smth) / jth.dens;

				drvt[i].gradBperp_x += jth.mass
						* (MagneticBPerpj.x - MagneticBPerpi.x)
						* kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradBperp_y += jth.mass
						* (MagneticBPerpj.y - MagneticBPerpi.y)
						* kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradBperp_z += jth.mass
						* (MagneticBPerpj.z - MagneticBPerpi.z)
						* kernel.gradW(dr, ith.smth) / jth.dens;

				drvt[i].gradBperp2 += jth.mass * (Bperpj2 - Bperpi2)
						* kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradPT += jth.mass * (PTj - PTi)
						* kernel.gradW(dr, ith.smth) / jth.dens;

				drvt[i].ngrad_dens = ndir * drvt[i].grad_dens;
				drvt[i].ngradP = ndir * drvt[i].gradP;
				drvt[i].ngradPT = ndir * drvt[i].gradPT;
				drvt[i].ngradVpara = ndir * drvt[i].gradVpara;
				drvt[i].ngradBperp2 = ndir * drvt[i].gradBperp2;
				drvt[i].ngradvperp_x = ndir * drvt[i].gradvperp_x;
				drvt[i].ngradvperp_y = ndir * drvt[i].gradvperp_y;
				drvt[i].ngradvperp_z = ndir * drvt[i].gradvperp_z;
				drvt[i].ngradBperp_x = ndir * drvt[i].gradBperp_x;
				drvt[i].ngradBperp_y = ndir * drvt[i].gradBperp_y;
				drvt[i].ngradBperp_z = ndir * drvt[i].gradBperp_z;
				drvt[i].nrot_v = ndir * drvt[i].rot_v;
				drvt[i].ngradV = ndir * drvt[i].gradV;
				drvt[i].ngraddens = ndir * drvt[i].graddens;
				drvt[i].ngradpres = ndir * drvt[i].gradpres;
				drvt[i].ngradvel_x = ndir * drvt[i].gradvel_x;
				drvt[i].ngradvel_y = ndir * drvt[i].gradvel_y;
				drvt[i].ngradvel_z = ndir * drvt[i].gradvel_z;
				drvt[i].ngradvel = ndir * drvt[i].gradvel;

			}

		}
	}
};

class CalcHydroForce {
	const kernel_t kernel;
public:
	PS::S32 sgn(PS::F64 v) {
		return (v > 0) - (v < 0);
	}
	void calc_Vij2_and_ss(const EPI::Hydro ep_i, const EPJ::Hydro ep_j,
			PS::F64&Vij2_hi, PS::F64&Vij2_hj, PS::F64& ssP, const PS::F64 delta,
			const PS::F64vec eij) {

		const PS::F64 Vi = 1.0 / ep_i.dens;
		const PS::F64 Vj = 1.0 / ep_j.dens;
		const PS::F64 hi2 = ep_i.smth * ep_i.smth;
		const PS::F64 hj2 = ep_j.smth * ep_j.smth;
		const PS::F64 Pij = 0.5 * (ep_i.pres + ep_j.pres);
		PS::F64 Aij, Bij, Cij, Dij;
		Cij = (Vi - Vj) / delta;
		Dij = 0.5 * (Vi + Vj);
		Vij2_hi = 0.25 * hi2 * Cij * Cij + Dij * Dij;
		Vij2_hj = 0.25 * hj2 * Cij * Cij + Dij * Dij;
		ssP = (hi2 * Cij * Dij / (2 * Vij2_hi));
		//		const PS::F64 rhoij = 0.5 * (ep_i.dens + ep_j.dens);
		//		const PS::F64 Vij2crit = 10.0 * (.0 / (rhoij * rhoij));
		//		const PS::F64 dVi = ep_i.gradV * eij;
		//		const PS::F64 dVj = ep_j.gradV * eij;
		//		Aij = -2 * (Vi - Vj) * pow(delta, -3) + (dVi + dVj) * pow(delta, -2);
		//		Bij = 0.5 * (dVi - dVj) / delta;
		//		Cij = 1.5 * (Vi - Vj) / delta - 0.25 * (dVi + dVj);
		//		Dij = 0.5 * (Vi + Vj) - 0.125 * (dVi - dVj) * delta;
		//		Vij2_hi = ((0.234375 * hi2 * Aij * Aij + 0.1875 * (2 * Aij * Cij + Bij * Bij)) * hi2 + 0.25 * (2 * Bij * Dij + Cij * Cij)) * hi2 + Dij * Dij;
		//		Vij2_hj = ((0.234375 * hj2 * Aij * Aij + 0.1875 * (2 * Aij * Cij + Bij * Bij)) * hj2 + 0.25 * (2 * Bij * Dij + Cij * Cij)) * hj2 + Dij * Dij;
		//		ssP = 0.5 * (((0.46875 * hi2 * Aij * Bij + 0.375 * (Aij * Dij + Bij * Cij)) * hi2 + 0.5 * Cij * Dij) * hi2 / Vij2_hi + ((0.46875 * hj2 * Aij * Bij + 0.375 * (Aij * Dij + Bij * Cij)) * hj2 + 0.5 * Cij * Dij) * hj2 / Vij2_hj);
	}
	void operator ()(const EPI::Hydro* const ep_i, const PS::S32 Nip,
			const EPJ::Hydro* const ep_j, const PS::S32 Njp,
			RESULT::Hydro* const hydro) {
		for (PS::S32 i = 0; i < Nip; ++i) {
			hydro[i].clear();
			PS::F64 Gamma = PARAM::GAMMA;
			const EPI::Hydro& ith = ep_i[i];

			for (PS::S32 j = 0; j < Njp; ++j) {
				const EPJ::Hydro& jth = ep_j[j];

				const PS::F64vec dr = ith.pos - jth.pos;
				if (dr == 0.0) {
					continue;
				}
				const PS::F64 dr_norm = sqrt(dr * dr);
				const PS::F64vec ndir = dr / dr_norm;
				const PS::F64vec dv = ith.vel - jth.vel;

				const PS::F64 s_ij = dr_norm;
				PS::F64 VIJ_I, VIJ_J, sstar;
//				calc_Vij2_and_ss(ith, jth, VIJ_I, VIJ_J, sstar, dr_norm, ndir);
				const PS::F64 hbar_ij = .5 * (ith.smth + jth.smth);
//				const PS::F64 F_ij = (kernel.gradW_scoord(s_ij, sqrt(2.0) * ith.smth)  *VIJ_I + kernel.gradW_scoord(s_ij, sqrt(2.0) * jth.smth)  *VIJ_J);
				const PS::F64 F_ij = (kernel.gradW_scoord(s_ij, hbar_ij)
						/ (ith.dens * ith.dens)
						+ kernel.gradW_scoord(s_ij, hbar_ij)
								/ (jth.dens * jth.dens));

				PS::F64 Pt_RP = 0.0;
				PS::F64 VparaRP = 0.0;
				PS::F64vec BperpMoC = 0.0;
				PS::F64vec VperpMoC = 0.0;
				PS::F64 P_Bpara = 0.0;
				PS::F64 BparaStar = 0.0;

				calcRPMoC(ith.dt, Pt_RP, VparaRP, BperpMoC, VperpMoC, P_Bpara,
						BparaStar, Gamma, ndir, ith, jth);
				const PS::F64vec fac_acc_correc = BparaStar * ith.MagneticB;

				const PS::F64vec fac_acc = -(Pt_RP - P_Bpara) * ndir
						+ BparaStar * BperpMoC - fac_acc_correc;

				hydro[i].acc += jth.mass * F_ij * fac_acc;

			}
			const PS::F64vec vhalf = ith.vel * .5 * hydro[i].acc * ith.dt;
			PS::F64 divB = 0.0;
			for (PS::S32 j = 0; j < Njp; ++j) {
				const EPJ::Hydro& jth = ep_j[j];

				const PS::F64vec dr = ith.pos - jth.pos;
				if (dr == 0.0) {
					continue;
				}
				const PS::F64 dr_norm = sqrt(dr * dr);
				const PS::F64vec ndir = dr / dr_norm;
				const PS::F64vec dv = ith.vel - jth.vel;

				const PS::F64 s_ij = dr_norm;
				PS::F64 VIJ_I, VIJ_J, sstar;
				//				calc_Vij2_and_ss(ith, jth, VIJ_I, VIJ_J, sstar, dr_norm, ndir);
				const PS::F64 hbar_ij = .5 * (ith.smth + jth.smth);
				//				const PS::F64 F_ij = (kernel.gradW_scoord(s_ij, sqrt(2.0) * ith.smth)  *VIJ_I + kernel.gradW_scoord(s_ij, sqrt(2.0) * jth.smth)  *VIJ_J);
				const PS::F64 F_ij = (kernel.gradW_scoord(s_ij, hbar_ij)
						/ (ith.dens * ith.dens)
						+ kernel.gradW_scoord(s_ij, hbar_ij)
								/ (jth.dens * jth.dens));

				PS::F64 Pt_RP = 0.0;
				PS::F64 VparaRP = 0.0;
				PS::F64vec BperpMoC = 0.0;
				PS::F64vec VperpMoC = 0.0;
				PS::F64 P_Bpara = 0.0;
				PS::F64 BparaStar = 0.0;

				calcRPMoC(ith.dt, Pt_RP, VparaRP, BperpMoC, VperpMoC, P_Bpara,
						BparaStar, Gamma, ndir, ith, jth);
				const PS::F64vec fac_acc_correc = BparaStar * ith.MagneticB;
				const PS::F64 fac_eng_dot_correc = BparaStar * vhalf
						* ith.MagneticB;

				const PS::F64 fac_eng_dot = -(Pt_RP - P_Bpara) * VparaRP
						+ BparaStar * (BperpMoC * VperpMoC)
						- fac_eng_dot_correc;

				hydro[i].eng_dot += jth.mass * F_ij * fac_eng_dot;

				const PS::F64vec fac_BoverDens_dot = BparaStar
						* (VparaRP * ndir + VperpMoC - vhalf);
				hydro[i].BoverDens_dot += jth.mass * F_ij * fac_BoverDens_dot;
				hydro[i].div_B += jth.mass * F_ij * BparaStar;

			}

			const PS::F64 C_fast = sqrt(
					(Gamma * ith.pres + ith.MagneticB * ith.MagneticB)
							/ ith.dens);
			PS::F64 dt = PARAM::C_CFL * ith.smth / C_fast;
			hydro[i].dt = 0.6 * dt;

		}

	}
	double trimin(double x, double y, double z) {
		return x < y ? (x < z ? x : z) : (y < z ? y : z);
	}
	//i -> R; j->L
	void calcRPMoC(const PS::F64 dt, PS::F64 &Pt_RP, PS::F64 &V_RP,
			PS::F64vec &BperpMoC, PS::F64vec &VperpMoC, PS::F64 &P_Bpara,
			PS::F64 &BparaStar, const PS::F64 gamma, const PS::F64vec ndir,
			const EPI::Hydro ep_R, const EPJ::Hydro ep_L) {
		const PS::F64 deltaS = (ep_R.pos - ep_L.pos) * ndir;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////Original Variables Start/////////////////////////////////////////////////////////////////////////////
		/**Left j**/
		const PS::F64 orig_vparaL = ep_L.vel * ndir;
		const PS::F64vec orig_MagneticBL = ep_L.MagneticB;
		const PS::F64vec orig_BperpL = ep_L.MagneticB
				- (ep_L.MagneticB * ndir) * ndir;
		const PS::F64 orig_BparaL = ep_L.MagneticB * ndir;
		const PS::F64 orig_PresL = ep_L.pres;
		const PS::F64 orig_DensL = ep_L.dens;
		/**Right i**/
		const PS::F64 orig_vparaR = ep_R.vel * ndir;
		const PS::F64vec orig_MagneticBR = ep_R.MagneticB;
		const PS::F64vec orig_BperpR = ep_R.MagneticB
				- (ep_R.MagneticB * ndir) * ndir;
		const PS::F64 orig_BparaR = ep_R.MagneticB * ndir;
		const PS::F64 orig_PresR = ep_R.pres;
		const PS::F64 orig_DensR = ep_R.dens;
//		if (orig_BparaL != orig_BparaR) {
//			std::cout << orig_vparaL << " | " << orig_MagneticBL << " | " << orig_BperpL << " | " << orig_BparaL << " | " << orig_PresL << " | " << orig_DensL << std::endl;
//			std::cout << orig_vparaR << " | " << orig_MagneticBR << " | " << orig_BperpR << " | " << orig_BparaR << " | " << orig_PresR << " | " << orig_DensR << std::endl;
//		}
		//////////////////////////Original Variables Finish/////////////////////////////////////////////////////////////////////////////
		//////////////////////////Original Derivative Variables Start/////////////////////////////////////////////////////////////////////////////
		/**Left j**/
		const PS::F64 drv_densL = ep_L.ngrad_dens;
		const PS::F64 drv_velparaL = ep_L.ngradVpara;
		const PS::F64 drv_BperpLx = ep_L.ngradBperp_x;
		const PS::F64 drv_BperpLy = ep_L.ngradBperp_y;
		const PS::F64 drv_BperpLz = ep_L.ngradBperp_z;
		const PS::F64 drv_PresL = ep_L.ngradP;

		/**Right i**/
		const PS::F64 drv_densR = ep_R.ngrad_dens;
		const PS::F64 drv_velparaR = ep_R.ngradVpara;
		const PS::F64 drv_BperpRx = ep_R.ngradBperp_x;
		const PS::F64 drv_BperpRy = ep_R.ngradBperp_y;
		const PS::F64 drv_BperpRz = ep_R.ngradBperp_z;
		const PS::F64 drv_PresR = ep_R.ngradP;

		//////////////////////////Original Derivative Variables Finish/////////////////////////////////////////////////////////////////////////////
		//////////////////////////Original Derived Variables Start/////////////////////////////////////////////////////////////////////////////
		/**Left j**/
		const PS::F64 orig_BperpL2 = orig_BperpL * orig_BperpL;
		const PS::F64 orig_BparaL2 = orig_BparaL * orig_BparaL;
		const PS::F64 orig_idensL = 1. / orig_DensL;
		const PS::F64 orig_densLsqrt = sqrt(orig_DensL);
		const PS::F64 orig_idensLsqrt = 1. / orig_densLsqrt;
		/**Right i**/
		const PS::F64 orig_BperpR2 = orig_BperpR * orig_BperpR;
		const PS::F64 orig_BparaR2 = orig_BparaR * orig_BparaR;
		const PS::F64 orig_idensR = 1. / orig_DensR;
		const PS::F64 orig_densRsqrt = sqrt(orig_DensR);
		const PS::F64 orig_idensRsqrt = 1. / orig_densRsqrt;

		P_Bpara = .25 * (orig_BparaR2 + orig_BparaL2);
		//////////////////////////Original Derived Variables Finish/////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////Monotone Variables Start////////////////////////////////////////////////////////////////////////////
		const PS::F64 orig_PTL = orig_PresL + .5 * orig_BperpL2;
		const PS::F64 orig_PTR = orig_PresR + .5 * orig_BperpR2;

		const PS::F64 DiffVparaR = orig_vparaR - orig_vparaL;
		const PS::F64 DiffPTR = orig_PTR - orig_PTL;
		const PS::F64 DiffPresR = orig_PresR - orig_PresL;
		const PS::F64 DiffDensR = orig_DensR - orig_DensL;
		const PS::F64 DiffBperpRx = orig_BperpR.x - orig_BperpL.x;
		const PS::F64 DiffBperpRy = orig_BperpR.y - orig_BperpL.y;
		const PS::F64 DiffBperpRz = orig_BperpR.z - orig_BperpL.z;

		const PS::F64 DiffVparaL = orig_vparaL - orig_vparaR;
		const PS::F64 DiffPTL = orig_PTL - orig_PTR;
		const PS::F64 DiffPresL = orig_PresL - orig_PresR;
		const PS::F64 DiffDensL = orig_DensL - orig_DensR;
		const PS::F64 DiffBperpLx = orig_BperpL.x - orig_BperpR.x;
		const PS::F64 DiffBperpLy = orig_BperpL.y - orig_BperpR.y;
		const PS::F64 DiffBperpLz = orig_BperpL.z - orig_BperpR.z;
		/**Left j**/
		const PS::F64 DeltaVparaL = drv_velparaL * deltaS;
		const PS::F64 DeltaDensL = drv_densL * deltaS;
		const PS::F64 DeltaPresL = drv_PresL * deltaS;

		const PS::F64 DeltaBperpLx = drv_BperpLx * deltaS;
		const PS::F64 DeltaBperpLy = drv_BperpLy * deltaS;
		const PS::F64 DeltaBperpLz = drv_BperpLz * deltaS;
		const PS::F64 DeltaVpara_dashL = 2.0 * DeltaVparaL - DiffVparaL;
		const PS::F64 DeltaDens_dashL = 2.0 * DeltaDensL - DiffDensL;
		const PS::F64 DiffPres_dashL = 2.0 * DeltaPresL - DiffPresL;

		const PS::F64 DiffBperp_dashLx = 2.0 * DeltaBperpLx - DiffBperpLx;
		const PS::F64 DiffBperp_dashLy = 2.0 * DeltaBperpLy - DiffBperpLy;
		const PS::F64 DiffBperp_dashLz = 2.0 * DeltaBperpLz - DiffBperpLz;
//		/**Right i**/
		const PS::F64 DeltaVparaR = drv_velparaR * deltaS;
		const PS::F64 DeltaDensR = drv_densR * deltaS;
		const PS::F64 DeltaBperpRx = drv_BperpRx * deltaS;
		const PS::F64 DeltaBperpRy = drv_BperpRy * deltaS;
		const PS::F64 DeltaBperpRz = drv_BperpRz * deltaS;
		const PS::F64 DeltaPresR = drv_PresR * deltaS;

		const PS::F64 DeltaVpara_dashR = 2.0 * DeltaVparaR - DiffVparaR;
		const PS::F64 DeltaDens_dashR = 2.0 * DeltaDensR - DiffDensR;
		const PS::F64 DiffPres_dashR = 2.0 * DeltaPresR - DiffPresR;

		const PS::F64 DiffBperp_dashRx = 2.0 * DeltaBperpRx - DiffBperpRx;
		const PS::F64 DiffBperp_dashRy = 2.0 * DeltaBperpRy - DiffBperpRy;
		const PS::F64 DiffBperp_dashRz = 2.0 * DeltaBperpRz - DiffBperpRz;

//		//////////////////////////Monotone Variables Finish////////////////////////////////////////////////////////////////////////////
//		//////////////////////////Monotone Constraints Start////////////////////////////////////////////////////////////////////////////
		const PS::F64 DeltaVparaL_mono =
				(sgn(DiffVparaL) == sgn(DeltaVparaL)
						&& sgn(DeltaVparaL) == sgn(DeltaVpara_dashL)) ?
						sgn(DeltaVparaL)
								* trimin(2.0 * fabs(DiffVparaL),
										fabs(DeltaVparaL),
										2.0 * fabs(DeltaVpara_dashL)) :
						0.0;
		const PS::F64 DeltaDensL_mono =
				(sgn(DiffDensL) == sgn(DeltaDensL)
						&& sgn(DeltaDensL) == sgn(DeltaDens_dashL)) ?
						sgn(DeltaDensL)
								* trimin(2.0 * fabs(DiffDensL),
										fabs(DeltaDensL),
										2.0 * fabs(DeltaDens_dashL)) :
						0.0;
		const PS::F64 DeltaPresL_mono =
				(sgn(DiffPresL) == sgn(DeltaPresL)
						&& sgn(DeltaPresL) == sgn(DiffPres_dashL)) ?
						sgn(DeltaPresL)
								* trimin(2.0 * fabs(DiffPresL),
										fabs(DeltaPresL),
										2.0 * fabs(DiffPres_dashL)) :
						0.0;

		const PS::F64 DeltaBperpLx_mono =
				(sgn(DiffBperpLx) == sgn(DeltaBperpLx)
						&& sgn(DeltaBperpLx) == sgn(DiffBperp_dashLx)) ?
						sgn(DeltaBperpLx)
								* trimin(2.0 * fabs(DiffBperpLx),
										fabs(DeltaBperpLx),
										2.0 * fabs(DiffBperp_dashLx)) :
						0.0;
		const PS::F64 DeltaBperpLy_mono =
				(sgn(DiffBperpLy) == sgn(DeltaBperpLy)
						&& sgn(DeltaBperpLy) == sgn(DiffBperp_dashLy)) ?
						sgn(DeltaBperpLy)
								* trimin(2.0 * fabs(DiffBperpLy),
										fabs(DeltaBperpLy),
										2.0 * fabs(DiffBperp_dashLy)) :
						0.0;
		const PS::F64 DeltaBperpLz_mono =
				(sgn(DiffBperpLz) == sgn(DeltaBperpLz)
						&& sgn(DeltaBperpLz) == sgn(DiffBperp_dashLz)) ?
						sgn(DeltaBperpLz)
								* trimin(2.0 * fabs(DiffBperpLz),
										fabs(DeltaBperpLz),
										2.0 * fabs(DiffBperp_dashLz)) :
						0.0;
		const PS::F64 DeltaBperpRx_mono =
				(sgn(DiffBperpRx) == sgn(DeltaBperpRx)
						&& sgn(DeltaBperpRx) == sgn(DiffBperp_dashRx)) ?
						sgn(DeltaBperpRx)
								* trimin(2.0 * fabs(DiffBperpRx),
										fabs(DeltaBperpRx),
										2.0 * fabs(DiffBperp_dashRx)) :
						0.0;
		const PS::F64 DeltaBperpRy_mono =
				(sgn(DiffBperpRy) == sgn(DeltaBperpRy)
						&& sgn(DeltaBperpRy) == sgn(DiffBperp_dashRy)) ?
						sgn(DeltaBperpRy)
								* trimin(2.0 * fabs(DiffBperpRy),
										fabs(DeltaBperpRy),
										2.0 * fabs(DiffBperp_dashRy)) :
						0.0;
		const PS::F64 DeltaBperpRz_mono =
				(sgn(DiffBperpRz) == sgn(DeltaBperpRz)
						&& sgn(DeltaBperpRz) == sgn(DiffBperp_dashRz)) ?
						sgn(DeltaBperpRz)
								* trimin(2.0 * fabs(DiffBperpRz),
										fabs(DeltaBperpRz),
										2.0 * fabs(DiffBperp_dashRz)) :
						0.0;

		const PS::F64 DeltaVparaR_mono =
				(sgn(DiffVparaR) == sgn(DeltaVparaR)
						&& sgn(DeltaVparaR) == sgn(DeltaVpara_dashR)) ?
						sgn(DeltaVparaR)
								* trimin(2.0 * fabs(DiffVparaR),
										fabs(DeltaVparaR),
										2.0 * fabs(DeltaVpara_dashR)) :
						0.0;
		const PS::F64 DeltaDensR_mono =
				(sgn(DiffDensR) == sgn(DeltaDensR)
						&& sgn(DeltaDensR) == sgn(DeltaDens_dashR)) ?
						sgn(DeltaDensR)
								* trimin(2.0 * fabs(DiffDensR),
										fabs(DeltaDensR),
										2.0 * fabs(DeltaDens_dashR)) :
						0.0;
		const PS::F64 DeltaPresR_mono =
				(sgn(DiffPresR) == sgn(DeltaPresR)
						&& sgn(DeltaPresR) == sgn(DiffPres_dashR)) ?
						sgn(DeltaPresR)
								* trimin(2.0 * fabs(DiffPresR),
										fabs(DeltaPresR),
										2.0 * fabs(DiffPres_dashR)) :
						0.0;

		/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
////////////////No mono constraints//////////////////////////////
//		const PS::F64 DeltaBperpLx_mono = DeltaBperpLx;
//		const PS::F64 DeltaBperpLy_mono = DeltaBperpLy;
//		const PS::F64 DeltaBperpRx_mono = DeltaBperpRx;
//		const PS::F64 DeltaBperpRy_mono = DeltaBperpRy;

//		const PS::F64 DeltaVparaL_mono = DeltaVparaL;
//		const PS::F64 DeltaBperpL2_mono = DeltaBperpL2;
//		const PS::F64 DeltaPTL_mono =DeltaPTL;
////		const PS::F64 DeltaPresL_mono = DeltaPresL;
//		const PS::F64 DeltaDensL_mono = DeltaDensL;
//		const PS::F64 DeltaVparaR_mono = DeltaVparaR;
//		const PS::F64 DeltaBperpR2_mono = DeltaBperpR2;
//		const PS::F64 DeltaPTR_mono = DeltaPTR;
////		const PS::F64 DeltaPresR_mono = DeltaPresR;
//		const PS::F64 DeltaDensR_mono = DeltaDensR;
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//		//////////////////////////Monotone Constraints Finish////////////////////////////////////////////////////////////////////////////
//		////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		///////////////////////////////////////////////////////////////////////////////////////////////////////////
//		////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		///////////////////////////////////////////////////////////////////////////////////////////////////////////
//		////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//		//////////////////////////RP Variables////////////////////////////////
		const PS::F64 C_fastL = sqrt(
				(gamma * ep_L.pres + orig_MagneticBL * orig_MagneticBL)
						* orig_idensL);
		const PS::F64 dod_facLRP = .5 * (1.0 - C_fastL * dt / deltaS);
		const PS::F64 C_fastR = sqrt(
				(gamma * ep_R.pres + orig_MagneticBR * orig_MagneticBR)
						* orig_idensR);
		const PS::F64 dod_facRRP = .5 * (1.0 - C_fastR * dt / deltaS);
//
		const PS::F64 fac_velparaLRP = DeltaVparaL_mono * dod_facLRP;
		const PS::F64 fac_BperpLxRP = DeltaBperpLx_mono * dod_facLRP;
		const PS::F64 fac_BperpLyRP = DeltaBperpLy_mono * dod_facLRP;
		const PS::F64 fac_BperpLzRP = DeltaBperpLz_mono * dod_facLRP;
		const PS::F64 fac_densLRP = DeltaDensL_mono * dod_facLRP;
		const PS::F64 fac_PresLRP = DeltaPresL_mono * dod_facLRP;

		const PS::F64 fac_velparaRRP = DeltaVparaR_mono * dod_facRRP;
		const PS::F64 fac_BperpRxRP = DeltaBperpRx_mono * dod_facRRP;
		const PS::F64 fac_BperpRyRP = DeltaBperpRy_mono * dod_facRRP;
		const PS::F64 fac_BperpRzRP = DeltaBperpRz_mono * dod_facRRP;
		const PS::F64 fac_densRRP = DeltaDensR_mono * dod_facRRP;
		const PS::F64 fac_PresRRP = DeltaPresR_mono * dod_facRRP;







//		//////L -> j///////////////////////////
		PS::F64 vparaLRP = orig_vparaL + fac_velparaLRP;
		PS::F64 BperpLxRP = orig_BperpL.x + fac_BperpLxRP;
		PS::F64 BperpLyRP = orig_BperpL.y + fac_BperpLyRP;
		PS::F64 BperpLzRP = orig_BperpL.z + fac_BperpLzRP;
		PS::F64 BperpL2RP = BperpLxRP * BperpLxRP + BperpLyRP * BperpLyRP
				+ BperpLzRP * BperpLzRP;
		PS::F64 densLRP = orig_DensL + fac_densLRP;
		PS::F64 presLRP = orig_PresL + fac_PresLRP;

		const PS::F64 PTLRP = presLRP + .5 * BperpL2RP;
//
//		//////R -> i///////////////////////////
//
		PS::F64 vparaRRP = orig_vparaR - fac_velparaRRP;
		PS::F64 BperpRxRP = orig_BperpR.x - fac_BperpRxRP;
		PS::F64 BperpRyRP = orig_BperpR.y - fac_BperpRyRP;
		PS::F64 BperpRzRP = orig_BperpR.z - fac_BperpRzRP;
		PS::F64 BperpR2RP = BperpRxRP * BperpRxRP + BperpRyRP * BperpRyRP
				+ BperpRzRP * BperpRzRP;
		PS::F64 densRRP = orig_DensR - fac_densRRP;
		PS::F64 presRRP = orig_PresR - fac_PresRRP;

		const PS::F64 PTRRP = presRRP + .5 * BperpR2RP;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////L -> j///////////////////////////
//		PS::F64 vparaLRP = orig_vparaL;
//		PS::F64 BperpL2RP = orig_BperpL2;
//		PS::F64 densLRP = orig_DensL;

//		const PS::F64 PTLRP = PTLRP_star;
//		////////R -> i///////////////////////////
//		PS::F64 vparaRRP = orig_vparaR;
//		PS::F64 BperpR2RP = orig_BperpR2;
//		PS::F64 densRRP = orig_DensR;
//		const PS::F64 PTRRP = PTRRP_star;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		const PS::F64 PTLRP_star = orig_PresL + .5 * orig_BperpL2;
		const PS::F64 PTRRP_star = orig_PresR + .5 * orig_BperpR2;
		PS::F64 ptRRP = PTRRP;
		PS::F64 ptLRP = PTLRP;
		if (ptRRP < 0.0 || ptLRP < 0.0 || densLRP < 0.0 || densRRP < 0.0
				|| BperpL2RP < 0.0 || BperpR2RP < 0.0) {
			ptRRP = PTRRP_star;
			ptLRP = PTLRP_star;
			BperpR2RP = orig_BperpR2;
			BperpL2RP = orig_BperpL2;
			densRRP = orig_DensR;
			densLRP = orig_DensL;
			vparaRRP = orig_vparaR;
			vparaLRP = orig_vparaL;
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		PS::F64 ptStar = 0.5 * (ptLRP + ptRRP);
		PS::F64 ptStar_old = 0.0;
		PS::F64 loop = 1;
		while (loop) {
			loop += 1;
			ptStar_old = ptStar;
			PS::F64 DR2 = pow((gamma + 1) * ptRRP + (gamma - 1) * ptStar, 2.)
					- (gamma - 2) * BperpR2RP
							* (2. * (gamma - 3) * ptRRP
									+ 2. * (gamma + 3) * ptStar
									- (gamma - 2) * BperpR2RP);
			PS::F64 DL2 = pow((gamma + 1) * ptLRP + (gamma - 1) * ptStar, 2.)
					- (gamma - 2) * BperpL2RP
							* (2. * (gamma - 3) * ptLRP
									+ 2. * (gamma + 3) * ptStar
									- (gamma - 2) * BperpL2RP);
			PS::F64 M_R2 = .25 * densRRP
					* ((gamma - 3.0) * ptRRP + (gamma + 3) * ptStar
							- (gamma - 2) * BperpR2RP + sqrt(DR2));
			PS::F64 M_L2 = .25 * densLRP
					* ((gamma - 3.0) * ptLRP + (gamma + 3) * ptStar
							- (gamma - 2) * BperpL2RP + sqrt(DL2));
			const PS::F64 M_R = sqrt(M_R2);
			const PS::F64 M_L = sqrt(M_L2);
			ptStar = (ptRRP / M_R + ptLRP / M_L - (vparaRRP - vparaLRP))
					/ (1. / M_R + 1. / M_L);
			if (ptStar != ptStar) {
				ptRRP = PTRRP_star;
				ptLRP = PTLRP_star;
				BperpR2RP = orig_BperpR2;
				BperpL2RP = orig_BperpL2;
				densRRP = orig_DensR;
				densLRP = orig_DensL;
				vparaRRP = orig_vparaR;
				vparaLRP = orig_vparaL;
				ptStar = 0.5 * (PTRRP_star + PTLRP_star);

			}
			if (ptStar < 0.0) {
//				std::cout << loop << " otStar:" << ptStar << "otStar:" << ptRRP << " " << M_R2 << " " << ptLRP << " " << M_L2 << " " << sqrt(DL2) << " " << (gamma - 3.0) * ptLRP + (gamma + 3) * ptStar - (gamma - 2) * BperpL2RP << std::endl;
				ptStar = 0.5 * (PTRRP_star + PTLRP_star);
				V_RP = .5 * (orig_vparaL + orig_vparaR);
				Pt_RP = ptStar;

				break;
			}
			if (fabs(ptStar - ptStar_old) < 0.01) {
				V_RP = (vparaRRP * M_R + vparaLRP * M_L - (ptRRP - ptLRP))
						/ (M_R + M_L);
				Pt_RP = ptStar;

				break;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////MoC Start///////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		const PS::F64vec orig_vperpL = ep_L.vel - (ep_L.vel * ndir) * ndir;
		const PS::F64vec orig_vperpR = ep_R.vel - (ep_R.vel * ndir) * ndir;
		const PS::F64 drv_vperpxL = ep_L.gradvperp_x * ndir;
		const PS::F64 drv_vperpyL = ep_L.gradvperp_y * ndir;
		const PS::F64 drv_vperpzL = ep_L.gradvperp_z * ndir;
		const PS::F64 drv_BperpxL = ep_L.gradBperp_x * ndir;
		const PS::F64 drv_BperpyL = ep_L.gradBperp_y * ndir;
		const PS::F64 drv_BperpzL = ep_L.gradBperp_z * ndir;
		const PS::F64 drv_vperpxR = ep_R.gradvperp_x * ndir;
		const PS::F64 drv_vperpyR = ep_R.gradvperp_y * ndir;
		const PS::F64 drv_vperpzR = ep_R.gradvperp_z * ndir;
		const PS::F64 drv_BperpxR = ep_R.gradBperp_x * ndir;
		const PS::F64 drv_BperpyR = ep_R.gradBperp_y * ndir;
		const PS::F64 drv_BperpzR = ep_R.gradBperp_z * ndir;
//
		const PS::F64vec DiffVperp = orig_vperpR - orig_vperpL;
		const PS::F64vec DiffBperp = orig_BperpR - orig_BperpL;
////
//		const PS::F64 DeltaVperpxL = drv_vperpxL * deltaS;
//		const PS::F64 DeltaVperpyL = drv_vperpyL * deltaS;
//		const PS::F64 DeltaVperpzL = drv_vperpzL * deltaS;
//		const PS::F64 DeltaBperpxL = drv_BperpxL * deltaS;
//		const PS::F64 DeltaBperpyL = drv_BperpyL * deltaS;
//		const PS::F64 DeltaBperpzL = drv_BperpzL * deltaS;
//		const PS::F64 DeltaVperpxR = drv_vperpxR * deltaS;
//		const PS::F64 DeltaVperpyR = drv_vperpyR * deltaS;
//		const PS::F64 DeltaVperpzR = drv_vperpzR * deltaS;
//		const PS::F64 DeltaBperpxR = drv_BperpxR * deltaS;
//		const PS::F64 DeltaBperpyR = drv_BperpyR * deltaS;
//		const PS::F64 DeltaBperpzR = drv_BperpzR * deltaS;
//		const PS::F64vec DeltaVperpL(DeltaVperpxL, DeltaVperpyL, DeltaVperpzL);
//		const PS::F64vec DeltaBperpL(DeltaBperpxL, DeltaBperpyL, DeltaBperpzL);
//		const PS::F64vec DeltaVperpR(DeltaVperpxR, DeltaVperpyR, DeltaBperpzL);
//		const PS::F64vec DeltaBperpR(DeltaBperpxR, DeltaBperpyR, DeltaBperpzL);
//		const PS::F64vec DiffBperp_dashL = 2.0 * DeltaBperpL + DiffBperp;
//		const PS::F64vec DiffVperp_dashL = 2.0 * DeltaVperpL + DiffVperp;
//		const PS::F64vec DiffBperp_dashR = 2.0 * DeltaBperpR - DiffBperp;
//		const PS::F64vec DiffVperp_dashR = 2.0 * DeltaVperpR - DiffVperp;
//
//		const PS::F64 DeltaVperpxL_mono = (-sgn(DiffVperp.x) == sgn(DeltaVperpxL) && sgn(DeltaVperpxL) == sgn(DiffVperp_dashL.x)) ? sgn(DeltaVperpxL) * trimin(2.0 * fabs(DiffVperp.x), fabs(DeltaVperpxL), 2.0 * fabs(DiffVperp_dashL.x)) : 0.0;
//		const PS::F64 DeltaVperpyL_mono = (-sgn(DiffVperp.y) == sgn(DeltaVperpyL) && sgn(DeltaVperpyL) == sgn(DiffVperp_dashL.y)) ? sgn(DeltaVperpyL) * trimin(2.0 * fabs(DiffVperp.y), fabs(DeltaVperpyL), 2.0 * fabs(DiffVperp_dashL.y)) : 0.0;
//		const PS::F64 DeltaVperpzL_mono = (-sgn(DiffVperp.z) == sgn(DeltaVperpzL) && sgn(DeltaVperpzL) == sgn(DiffVperp_dashL.z)) ? sgn(DeltaVperpzL) * trimin(2.0 * fabs(DiffVperp.z), fabs(DeltaVperpzL), 2.0 * fabs(DiffVperp_dashL.z)) : 0.0;
//
//		const PS::F64 DeltaBperpxL_mono = (-sgn(DiffBperp.x) == sgn(DeltaBperpxL) && sgn(DeltaBperpxL) == sgn(DiffBperp_dashL.x)) ? sgn(DeltaBperpxL) * trimin(2.0 * fabs(DiffBperp.x), fabs(DeltaBperpxL), 2.0 * fabs(DiffBperp_dashL.x)) : 0.0;
//		const PS::F64 DeltaBperpyL_mono = (-sgn(DiffBperp.y) == sgn(DeltaBperpyL) && sgn(DeltaBperpyL) == sgn(DiffBperp_dashL.y)) ? sgn(DeltaBperpyL) * trimin(2.0 * fabs(DiffBperp.y), fabs(DeltaBperpyL), 2.0 * fabs(DiffBperp_dashL.y)) : 0.0;
//		const PS::F64 DeltaBperpzL_mono = (-sgn(DiffBperp.z) == sgn(DeltaBperpzL) && sgn(DeltaBperpzL) == sgn(DiffBperp_dashL.z)) ? sgn(DeltaBperpzL) * trimin(2.0 * fabs(DiffBperp.z), fabs(DeltaBperpzL), 2.0 * fabs(DiffBperp_dashL.z)) : 0.0;
//
//		const PS::F64 DeltaVperpxR_mono = (sgn(DiffVperp.x) == sgn(DeltaVperpxR) && sgn(DeltaVperpxR) == sgn(DiffVperp_dashR.x)) ? sgn(DeltaVperpxR) * trimin(2.0 * fabs(DiffVperp.x), fabs(DeltaVperpxR), 2.0 * fabs(DiffVperp_dashR.x)) : 0.0;
//		const PS::F64 DeltaVperpyR_mono = (sgn(DiffVperp.y) == sgn(DeltaVperpyR) && sgn(DeltaVperpyR) == sgn(DiffVperp_dashR.y)) ? sgn(DeltaVperpyR) * trimin(2.0 * fabs(DiffVperp.y), fabs(DeltaVperpyR), 2.0 * fabs(DiffVperp_dashR.y)) : 0.0;
//		const PS::F64 DeltaVperpzR_mono = (sgn(DiffVperp.z) == sgn(DeltaVperpzR) && sgn(DeltaVperpzR) == sgn(DiffVperp_dashR.z)) ? sgn(DeltaVperpzR) * trimin(2.0 * fabs(DiffVperp.z), fabs(DeltaVperpzR), 2.0 * fabs(DiffVperp_dashR.z)) : 0.0;
//
//		const PS::F64 DeltaBperpxR_mono = (sgn(DiffBperp.x) == sgn(DeltaBperpxR) && sgn(DeltaBperpxR) == sgn(DiffBperp_dashR.x)) ? sgn(DeltaBperpxR) * trimin(2.0 * fabs(DiffBperp.x), fabs(DeltaBperpxR), 2.0 * fabs(DiffBperp_dashR.x)) : 0.0;
//		const PS::F64 DeltaBperpyR_mono = (sgn(DiffBperp.y) == sgn(DeltaBperpyR) && sgn(DeltaBperpyR) == sgn(DiffBperp_dashR.y)) ? sgn(DeltaBperpyR) * trimin(2.0 * fabs(DiffBperp.y), fabs(DeltaBperpyR), 2.0 * fabs(DiffBperp_dashR.y)) : 0.0;
//		const PS::F64 DeltaBperpzR_mono = (sgn(DiffBperp.z) == sgn(DeltaBperpzR) && sgn(DeltaBperpzR) == sgn(DiffBperp_dashR.z)) ? sgn(DeltaBperpzR) * trimin(2.0 * fabs(DiffBperp.z), fabs(DeltaBperpzR), 2.0 * fabs(DiffBperp_dashR.z)) : 0.0;

//		const PS::F64 DeltaVperpxL_mono = DeltaVperpxL;
//		const PS::F64 DeltaVperpyL_mono = DeltaVperpyL;
//		const PS::F64 DeltaBperpxL_mono = DeltaBperpxL;
//		const PS::F64 DeltaBperpyL_mono = DeltaBperpyL;
//		const PS::F64 DeltaVperpxR_mono = DeltaVperpxR;
//		const PS::F64 DeltaVperpyR_mono = DeltaVperpyR;
//		const PS::F64 DeltaBperpxR_mono = DeltaBperpxR;
//		const PS::F64 DeltaBperpyR_mono = DeltaBperpyR;

		//ToDO Using Riemann Invariant
		const PS::F64vec orig_JPlusL = orig_vperpL
				- orig_BperpL * orig_idensLsqrt;
		const PS::F64vec orig_JMinusL = orig_vperpL
				+ orig_BperpL * orig_idensLsqrt;
		const PS::F64vec orig_JPlusR = orig_vperpR
				- orig_BperpR * orig_idensRsqrt;
		const PS::F64vec orig_JMinusR = orig_vperpR
				+ orig_BperpR * orig_idensRsqrt;
//
		const PS::F64 drv_JPlusxL = drv_vperpxL - drv_BperpxL * orig_idensLsqrt;
		const PS::F64 drv_JPlusyL = drv_vperpyL - drv_BperpyL * orig_idensLsqrt;
		const PS::F64 drv_JPluszL = drv_vperpzL - drv_BperpzL * orig_idensLsqrt;

		const PS::F64 drv_JMinusxL = drv_vperpxL
				+ drv_BperpxL * orig_idensLsqrt;
		const PS::F64 drv_JMinusyL = drv_vperpyL
				+ drv_BperpyL * orig_idensLsqrt;
		const PS::F64 drv_JMinuszL = drv_vperpzL
				+ drv_BperpzL * orig_idensLsqrt;

		const PS::F64 drv_JPlusxR = drv_vperpxR - drv_BperpxR * orig_idensRsqrt;
		const PS::F64 drv_JPlusyR = drv_vperpyR - drv_BperpyR * orig_idensRsqrt;
		const PS::F64 drv_JPluszR = drv_vperpzR - drv_BperpzR * orig_idensRsqrt;

		const PS::F64 drv_JMinusxR = drv_vperpxR
				+ drv_BperpxR * orig_idensRsqrt;
		const PS::F64 drv_JMinusyR = drv_vperpyR
				+ drv_BperpyR * orig_idensRsqrt;
		const PS::F64 drv_JMinuszR = drv_vperpzR
				+ drv_BperpzR * orig_idensRsqrt;

//
		const PS::F64vec DiffJPlus = orig_JPlusR - orig_JPlusL;
		const PS::F64vec DiffJMinus = orig_JMinusR - orig_JMinusL;
		const PS::F64 DeltaJPlusxL = drv_JPlusxL * deltaS;
		const PS::F64 DeltaJPlusyL = drv_JPlusyL * deltaS;
		const PS::F64 DeltaJPluszL = drv_JPluszL * deltaS;

		const PS::F64 DeltaJPlusxR = drv_JPlusxR * deltaS;
		const PS::F64 DeltaJPlusyR = drv_JPlusyR * deltaS;
		const PS::F64 DeltaJPluszR = drv_JPluszR * deltaS;

		const PS::F64 DeltaJMinusxL = drv_JMinusxL * deltaS;
		const PS::F64 DeltaJMinusyL = drv_JMinusyL * deltaS;
		const PS::F64 DeltaJMinuszL = drv_JMinuszL * deltaS;

		const PS::F64 DeltaJMinusxR = drv_JMinusxR * deltaS;
		const PS::F64 DeltaJMinusyR = drv_JMinusyR * deltaS;
		const PS::F64 DeltaJMinuszR = drv_JMinuszR * deltaS;

		const PS::F64vec DeltaJPlusL(DeltaJPlusxL, DeltaJPlusyL, DeltaJPluszL);
		const PS::F64vec DeltaJPlusR(DeltaJPlusxR, DeltaJPlusyR, DeltaJPluszR);
		const PS::F64vec DeltaJMinusL(DeltaJMinusxL, DeltaJMinusyL,
				DeltaJMinuszL);
		const PS::F64vec DeltaJMinusR(DeltaJMinusxR, DeltaJMinusyR,
				DeltaJMinuszR);
		const PS::F64vec DiffJPlus_dashL = 2.0 * DeltaJPlusL - DiffJPlus;
		const PS::F64vec DiffJMinus_dashL = 2.0 * DeltaJMinusL - DiffJMinus;
		const PS::F64vec DiffJPlus_dashR = 2.0 * DeltaJPlusR - DiffJPlus;
		const PS::F64vec DiffJMinus_dashR = 2.0 * DeltaJMinusR - DiffJMinus;

//		const PS::F64 DeltaJPlusxL_mono = (sgn(DiffJPlus.x) == sgn(DeltaJPlusxL) && sgn(DeltaJPlusxL) == sgn(DiffJPlus_dashL.x)) ? sgn(DeltaJPlusxL) * trimin(2.0 * fabs(DiffJPlus.x), fabs(DeltaJPlusxL), 2.0 * fabs(DiffJPlus_dashL.x)) : 0.0;
//		const PS::F64 DeltaJPlusyL_mono = (sgn(DiffJPlus.y) == sgn(DeltaJPlusyL) && sgn(DeltaJPlusyL) == sgn(DiffJPlus_dashL.y)) ? sgn(DeltaJPlusyL) * trimin(2.0 * fabs(DiffJPlus.y), fabs(DeltaJPlusyL), 2.0 * fabs(DiffJPlus_dashL.y)) : 0.0;
//		const PS::F64 DeltaJPluszL_mono = (sgn(DiffJPlus.z) == sgn(DeltaJPluszL) && sgn(DeltaJPluszL) == sgn(DiffJPlus_dashL.z)) ? sgn(DeltaJPluszL) * trimin(2.0 * fabs(DiffJPlus.z), fabs(DeltaJPlusyL), 2.0 * fabs(DiffJPlus_dashL.z)) : 0.0;
//
//		const PS::F64 DeltaJMinusxL_mono = (sgn(DiffJMinus.x) == sgn(DeltaJMinusxL) && sgn(DeltaJMinusxL) == sgn(DiffJMinus_dashL.x)) ? sgn(DeltaJMinusxL) * trimin(2.0 * fabs(DiffJMinus.x), fabs(DeltaJMinusxL), 2.0 * fabs(DiffJMinus_dashL.x)) : 0.0;
//		const PS::F64 DeltaJMinusyL_mono = (sgn(DiffJMinus.y) == sgn(DeltaJMinusyL) && sgn(DeltaJMinusyL) == sgn(DiffJMinus_dashL.y)) ? sgn(DeltaJMinusyL) * trimin(2.0 * fabs(DiffJMinus.y), fabs(DeltaJMinusyL), 2.0 * fabs(DiffJMinus_dashL.y)) : 0.0;
//		const PS::F64 DeltaJMinuszL_mono = (sgn(DiffJMinus.z) == sgn(DeltaJMinuszL) && sgn(DeltaJMinuszL) == sgn(DiffJMinus_dashL.z)) ? sgn(DeltaJMinuszL) * trimin(2.0 * fabs(DiffJMinus.z), fabs(DeltaJMinuszL), 2.0 * fabs(DiffJMinus_dashL.z)) : 0.0;
//
//		const PS::F64 DeltaJPlusxR_mono = (sgn(DiffJPlus.x) == sgn(DeltaJPlusxR) && sgn(DeltaJPlusxR) == sgn(DiffJPlus_dashR.x)) ? sgn(DeltaJPlusxR) * trimin(2.0 * fabs(DiffJPlus.x), fabs(DeltaJPlusxR), 2.0 * fabs(DiffJPlus_dashR.x)) : 0.0;
//		const PS::F64 DeltaJPlusyR_mono = (sgn(DiffJPlus.y) == sgn(DeltaJPlusyR) && sgn(DeltaJPlusyR) == sgn(DiffJPlus_dashR.y)) ? sgn(DeltaJPlusyR) * trimin(2.0 * fabs(DiffJPlus.y), fabs(DeltaJPlusyR), 2.0 * fabs(DiffJPlus_dashR.y)) : 0.0;
//		const PS::F64 DeltaJPluszR_mono = (sgn(DiffJPlus.z) == sgn(DeltaJPluszR) && sgn(DeltaJPluszR) == sgn(DiffJPlus_dashR.z)) ? sgn(DeltaJPluszR) * trimin(2.0 * fabs(DiffJPlus.z), fabs(DeltaJPluszR), 2.0 * fabs(DiffJPlus_dashR.z)) : 0.0;
//
//		const PS::F64 DeltaJMinusxR_mono = (sgn(DiffJMinus.x) == sgn(DeltaJMinusxR) && sgn(DeltaJMinusxR) == sgn(DiffJMinus_dashR.x)) ? sgn(DeltaJMinusxR) * trimin(2.0 * fabs(DiffJMinus.x), fabs(DeltaJMinusxR), 2.0 * fabs(DiffJMinus_dashR.x)) : 0.0;
//		const PS::F64 DeltaJMinusyR_mono = (sgn(DiffJMinus.y) == sgn(DeltaJMinusyR) && sgn(DeltaJMinusyR) == sgn(DiffJMinus_dashR.y)) ? sgn(DeltaJMinusyR) * trimin(2.0 * fabs(DiffJMinus.y), fabs(DeltaJMinusyR), 2.0 * fabs(DiffJMinus_dashR.y)) : 0.0;
//		const PS::F64 DeltaJMinuszR_mono = (sgn(DiffJMinus.z) == sgn(DeltaJMinuszR) && sgn(DeltaJMinuszR) == sgn(DiffJMinus_dashR.z)) ? sgn(DeltaJMinuszR) * trimin(2.0 * fabs(DiffJMinus.z), fabs(DeltaJMinuszR), 2.0 * fabs(DiffJMinus_dashR.z)) : 0.0;

		const PS::F64 DeltaJPlusxL_mono = DeltaJPlusxL;
		const PS::F64 DeltaJPlusyL_mono = DeltaJPlusyL;
		const PS::F64 DeltaJPluszL_mono = DeltaJPluszL;

		const PS::F64 DeltaJMinusxL_mono = DeltaJMinusxL;
		const PS::F64 DeltaJMinusyL_mono = DeltaJMinusyL;
		const PS::F64 DeltaJMinuszL_mono = DeltaJMinuszL;

		const PS::F64 DeltaJPlusxR_mono = DeltaJPlusxR;
		const PS::F64 DeltaJPlusyR_mono = DeltaJPlusyR;
		const PS::F64 DeltaJPluszR_mono = DeltaJPluszR;

		const PS::F64 DeltaJMinusxR_mono = DeltaJMinusxR;
		const PS::F64 DeltaJMinusyR_mono = DeltaJMinusyR;
		const PS::F64 DeltaJMinuszR_mono = DeltaJMinuszR;

		const PS::F64 C_alvenL = fabs(orig_BparaL) * orig_idensLsqrt;
		const PS::F64 C_alvenR = fabs(orig_BparaR) * orig_idensRsqrt;
		const PS::F64 dod_fac_AlfvenL = .5 * (1.0 - C_alvenL * dt / deltaS);
		const PS::F64 dod_fac_AlfvenR = .5 * (1.0 - C_alvenR * dt / deltaS);
		const PS::F64 fac_densLMoC = DeltaDensL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_densRMoC = DeltaDensR_mono * dod_fac_AlfvenR;
		const PS::F64 densL_guardMoC = orig_DensL + fac_densLMoC;
		const PS::F64 densR_guardMoC = orig_DensR - fac_densRMoC;
		const PS::F64 densLMoC =
				densL_guardMoC < 0.0 ? orig_DensL : densL_guardMoC;
		const PS::F64 densRMoC =
				densR_guardMoC < 0.0 ? orig_DensR : densR_guardMoC;

		const PS::F64 DeltaVperpxL_mono = .5
				* (DeltaJPlusxL_mono + DeltaJMinusxL_mono);
		const PS::F64 DeltaVperpyL_mono = .5
				* (DeltaJPlusyL_mono + DeltaJMinusyL_mono);
		const PS::F64 DeltaVperpzL_mono = .5
				* (DeltaJPluszL_mono + DeltaJMinuszL_mono);

		const PS::F64 DeltaBperpxL_mono = .5
				* (DeltaJMinusxL_mono - DeltaJPlusxL_mono) * sqrt(orig_DensL);
		const PS::F64 DeltaBperpyL_mono = .5
				* (DeltaJMinusyL_mono - DeltaJPlusyL_mono) * sqrt(orig_DensL);
		const PS::F64 DeltaBperpzL_mono = .5
				* (DeltaJMinuszL_mono - DeltaJPluszL_mono) * sqrt(orig_DensL);

		const PS::F64 DeltaVperpxR_mono = .5
				* (DeltaJPlusxR_mono + DeltaJMinusxR_mono);
		const PS::F64 DeltaVperpyR_mono = .5
				* (DeltaJPlusyR_mono + DeltaJMinusyR_mono);
		const PS::F64 DeltaVperpzR_mono = .5
				* (DeltaJPluszR_mono + DeltaJMinuszR_mono);

		const PS::F64 DeltaBperpxR_mono = .5
				* (DeltaJMinusxR_mono - DeltaJPlusxR_mono) * sqrt(orig_DensR);
		const PS::F64 DeltaBperpyR_mono = .5
				* (DeltaJMinusyR_mono - DeltaJPlusyR_mono) * sqrt(orig_DensR);
		const PS::F64 DeltaBperpzR_mono = .5
				* (DeltaJMinuszR_mono - DeltaJPluszR_mono) * sqrt(orig_DensR);

		////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////

		const PS::F64 fac_vperpxLMoC = DeltaVperpxL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_vperpyLMoC = DeltaVperpyL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_vperpzLMoC = DeltaVperpzL_mono * dod_fac_AlfvenL;

		const PS::F64 fac_BperpxLMoC = DeltaBperpxL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_BperpyLMoC = DeltaBperpyL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_BperpzLMoC = DeltaBperpzL_mono * dod_fac_AlfvenL;

		const PS::F64 fac_vperpxRMoC = DeltaVperpxR_mono * dod_fac_AlfvenR;
		const PS::F64 fac_vperpyRMoC = DeltaVperpyR_mono * dod_fac_AlfvenR;
		const PS::F64 fac_vperpzRMoC = DeltaVperpzR_mono * dod_fac_AlfvenR;

		const PS::F64 fac_BperpxRMoC = DeltaBperpxR_mono * dod_fac_AlfvenR;
		const PS::F64 fac_BperpyRMoC = DeltaBperpyR_mono * dod_fac_AlfvenR;
		const PS::F64 fac_BperpzRMoC = DeltaBperpzR_mono * dod_fac_AlfvenR;

//		const PS::F64 densL_guardMoC = orig_DensL;
//		const PS::F64 densR_guardMoC = orig_DensR;
////		const PS::F64 densLMoC = densL_guardMoC < 0.0 ? orig_DensL : densL_guardMoC;
////		const PS::F64 densRMoC = densR_guardMoC < 0.0 ? orig_DensR : densR_guardMoC;
//		const PS::F64 densLMoC = orig_DensL;
//		const PS::F64 densRMoC = orig_DensR;
		const PS::F64 vperpxLMoC = orig_vperpL.x;
		const PS::F64 vperpyLMoC = orig_vperpL.y;
		const PS::F64 vperpzLMoC = orig_vperpL.z;
		const PS::F64 vperpxRMoC = orig_vperpR.x;
		const PS::F64 vperpyRMoC = orig_vperpR.y;
		const PS::F64 vperpzRMoC = orig_vperpR.z;
		const PS::F64 BperpxLMoC = orig_BperpL.x;
		const PS::F64 BperpyLMoC = orig_BperpL.y;
		const PS::F64 BperpzLMoC = orig_BperpL.z;
		const PS::F64 BperpxRMoC = orig_BperpR.x;
		const PS::F64 BperpyRMoC = orig_BperpR.y;
		const PS::F64 BperpzRMoC = orig_BperpR.z;

//		const PS::F64 vperpxLMoC = orig_vperpL.x + fac_vperpxLMoC;
//		const PS::F64 vperpyLMoC = orig_vperpL.y + fac_vperpyLMoC;
//		const PS::F64 vperpzLMoC = orig_vperpL.z + fac_vperpzLMoC;
//
//		const PS::F64 BperpxLMoC = orig_BperpL.x + fac_BperpxLMoC;
//		const PS::F64 BperpyLMoC = orig_BperpL.y + fac_BperpyLMoC;
//		const PS::F64 BperpzLMoC = orig_BperpL.z + fac_BperpzLMoC;
//
//		const PS::F64 vperpxRMoC = orig_vperpR.x - fac_vperpxRMoC;
//		const PS::F64 vperpyRMoC = orig_vperpR.y - fac_vperpyRMoC;
//		const PS::F64 vperpzRMoC = orig_vperpR.z - fac_vperpzRMoC;
//
//		const PS::F64 BperpxRMoC = orig_BperpR.x - fac_BperpxRMoC;
//		const PS::F64 BperpyRMoC = orig_BperpR.y - fac_BperpyRMoC;
//		const PS::F64 BperpzRMoC = orig_BperpR.z - fac_BperpzRMoC;

		const PS::F64vec vperpLMoC(vperpxLMoC, vperpyLMoC, vperpzLMoC);
		const PS::F64vec BperpLMoC(BperpxLMoC, BperpyLMoC, BperpzLMoC);
		const PS::F64vec vperpRMoC(vperpxRMoC, vperpyRMoC, vperpzRMoC);
		const PS::F64vec BperpRMoC(BperpxRMoC, BperpyRMoC, BperpzRMoC);

		const PS::F64 idensLMoC = 1. / densLMoC;
		const PS::F64 densLMoCsqrt = sqrt(densLMoC);
		const PS::F64 idensLMoCsqrt = 1. / sqrt(densLMoC);
		const PS::F64 idensRMoC = 1. / densRMoC;
		const PS::F64 densRMoCsqrt = sqrt(densRMoC);
		const PS::F64 idensRMoCsqrt = 1. / sqrt(densRMoC);

		const PS::F64 _Bpara = .5 * (orig_BparaL + orig_BparaR);
		BperpMoC = (BperpRMoC * idensRMoCsqrt + BperpLMoC * idensLMoCsqrt
				+ sgn(_Bpara) * (-vperpLMoC + vperpRMoC))
				/ (idensLMoCsqrt + idensRMoCsqrt);
		VperpMoC = (densRMoCsqrt * vperpRMoC + densLMoCsqrt * vperpLMoC
				+ sgn(_Bpara) * (-BperpLMoC + BperpRMoC))
				/ (densLMoCsqrt + densRMoCsqrt);
//		BperpMoC = .5 * (BperpRMoC + BperpLMoC);

//		if(BperpMoC.x!=BperpMoC.x){
//							std::cout<<idensRMoCsqrt<<" "<<densRMoC<<" "<<idensLMoCsqrt<<" "<<densLMoC<<std::endl;
//						}
		BparaStar = _Bpara;
	}

};

template<class TParticleJ> class CalcGravityForce {
public:
	void operator ()(const EPI::Grav* const ep_i, const PS::S32 Nip,
			const TParticleJ* const ep_j, const PS::S32 Njp,
			RESULT::Grav* const grav) {
		for (PS::S32 i = 0; i < Nip; ++i) {
			const EPI::Grav& ith = ep_i[i];
			for (PS::S32 j = 0; j < Njp; ++j) {
				const TParticleJ& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64 dr2 = dr * dr;
				const PS::F64 dr_inv = 1.0 / sqrt(dr2 + ith.getEps2());
				const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
				grav[i].acc += -m_dr3_inv * dr;
				grav[i].pot += -jth.mass * dr_inv;
			}
		}
	}
};

