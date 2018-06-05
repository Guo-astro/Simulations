#pragma once

class CalcDensity {
	kernel_t kernel;
public:
	void operator ()(const EPI::Dens* const ep_i, const PS::S32 Nip, const EPJ::Dens* const ep_j, const PS::S32 Njp, RESULT::Dens* const dens) {
		for (PS::S32 i = 0; i < Nip; ++i) {
			dens[i].clear();
			const EPI::Dens& ith = ep_i[i];
			for (PS::S32 j = 0; j < Njp; ++j) {
				const EPJ::Dens& jth = ep_j[j];

				const PS::F64vec dr = ith.pos - jth.pos;
				dens[i].dens += jth.mass * kernel.W(dr, ith.smth);

			}
			dens[i].smth = PARAM::SMTH * pow(ith.mass / dens[i].dens, 1.0 / (PS::F64)(PARAM::Dim));
		}
	}
};

class CalcDerivative {
	kernel_t kernel;
public:
	void operator ()(const EPI::Drvt* ep_i, const PS::S32 Nip, const EPJ::Drvt* ep_j, const PS::S32 Njp, RESULT::Drvt* const drvt) {
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
				const PS::F64vec MagneticBPerpj = jth.MagneticB - (jth.MagneticB * ndir) * ndir;
				const PS::F64vec MagneticBPerpi = ith.MagneticB - (ith.MagneticB * ndir) * ndir;
				const PS::F64 Bperpj2 = MagneticBPerpj * MagneticBPerpj;
				const PS::F64 Bperpi2 = MagneticBPerpi * MagneticBPerpi;
				const PS::F64 PTj = jth.pres + .5 * MagneticBPerpj * MagneticBPerpj;
				const PS::F64 PTi = ith.pres + .5 * MagneticBPerpi * MagneticBPerpi;
				const PS::F64vec vperpj = jth.vel - (jth.vel * ndir) * ndir;
				const PS::F64vec vperpi = ith.vel - (ith.vel * ndir) * ndir;
				drvt[i].grad_dens += jth.mass * kernel.gradW(dr, ith.smth);
				drvt[i].gradP += jth.mass * (jth.pres - ith.pres) * kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradVpara += jth.mass * ((jth.vel - ith.vel) * ndir) * kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradvperp_x += jth.mass * (vperpj.x - vperpi.x) * kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradvperp_y += jth.mass * (vperpj.y - vperpi.y) * kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradBperp_x += jth.mass * (MagneticBPerpj.x - MagneticBPerpi.x) * kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradBperp_y += jth.mass * (MagneticBPerpj.y - MagneticBPerpi.y) * kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradBperp2 += jth.mass * (Bperpj2 - Bperpi2) * kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradPT += jth.mass * (PTj - PTi) * kernel.gradW(dr, ith.smth) / jth.dens;

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
	void operator ()(const EPI::Hydro* const ep_i, const PS::S32 Nip, const EPJ::Hydro* const ep_j, const PS::S32 Njp, RESULT::Hydro* const hydro) {
		for (PS::S32 i = 0; i < Nip; ++i) {
			hydro[i].clear();
			PS::F64 Gamma = 5. / 3.;
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
				const PS::F64 hbar_ij = .5 * (ith.smth + jth.smth);
				const PS::F64 F_ij = (1. / (ith.dens * ith.dens) + 1. / (jth.dens * jth.dens)) * kernel.gradW_scoord(s_ij, hbar_ij);

				PS::F64 Pt_RP = 0.0;
				PS::F64 VparaRP = 0.0;
				PS::F64vec BperpMoC = 0.0;
				PS::F64vec VperpMoC = 0.0;
				PS::F64 P_Bpara = 0.0;
				PS::F64 BparaStar = 0.0;

				calcRPMoC(ith.dt, Pt_RP, VparaRP, BperpMoC, VperpMoC, P_Bpara, BparaStar, Gamma, ndir, ith, jth);
				const PS::F64vec fac_acc_correc = BparaStar * ith.MagneticB;
				const PS::F64 fac_eng_dot_correc = BparaStar * ith.vel_half * ith.MagneticB;

				const PS::F64vec fac_acc = -(Pt_RP - P_Bpara) * ndir + BparaStar * BperpMoC - fac_acc_correc;
				const PS::F64 fac_eng_dot = -(Pt_RP - P_Bpara) * VparaRP + BparaStar * (BperpMoC * VperpMoC) - fac_eng_dot_correc;
				const PS::F64vec fac_BoverDens_dot = BparaStar * (VparaRP * ndir + VperpMoC - ith.vel_half);
				hydro[i].acc += jth.mass * F_ij * fac_acc;
				hydro[i].eng_dot += jth.mass * F_ij * fac_eng_dot;
				hydro[i].BoverDens_dot += jth.mass * F_ij * fac_BoverDens_dot;

			}
			const PS::F64 C_fast = sqrt((Gamma * ith.pres + ith.MagneticB * ith.MagneticB) / ith.dens);
			const PS::F64 dt = 0.1 * PARAM::C_CFL * ith.smth / C_fast;
			hydro[i].dt = dt;

		}
	}
	double trimin(double x, double y, double z) {
		return x < y ? (x < z ? x : z) : (y < z ? y : z);
	}
	//i -> R; j->L
	void calcRPMoC(const PS::F64 dt, PS::F64 &Pt_RP, PS::F64 &V_RP, PS::F64vec &BperpMoC, PS::F64vec &VperpMoC, PS::F64 &P_Bpara, PS::F64 &BparaStar, const PS::F64 gamma, const PS::F64vec ndir, const EPI::Hydro ep_R, const EPJ::Hydro ep_L) {
		const PS::F64 deltaS = sqrt((ep_R.pos - ep_L.pos) * (ep_R.pos - ep_L.pos));
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////Original Variables Start/////////////////////////////////////////////////////////////////////////////
		/**Left j**/
		const PS::F64 orig_vparaL = ep_L.vel * ndir;
		const PS::F64vec orig_MagneticBL = ep_L.MagneticB;
		const PS::F64vec orig_BperpL = ep_L.MagneticB - (ep_L.MagneticB * ndir) * ndir;
		const PS::F64 orig_BparaL = ep_L.MagneticB * ndir;
		const PS::F64 orig_PresL = ep_L.pres;
		const PS::F64 orig_DensL = ep_L.dens;
		/**Right i**/
		const PS::F64 orig_vparaR = ep_R.vel * ndir;
		const PS::F64vec orig_MagneticBR = ep_R.MagneticB;
		const PS::F64vec orig_BperpR = ep_R.MagneticB - (ep_R.MagneticB * ndir) * ndir;
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
		const PS::F64 drv_densL = ep_L.grad_dens * ndir;
		const PS::F64 drv_PresL = ep_L.gradP * ndir;
		const PS::F64 drv_PTL = ep_L.gradPT * ndir;
		const PS::F64 drv_velparaL = ep_L.gradVpara * ndir;
		const PS::F64 drv_BperpL2 = ep_L.gradBperp2 * ndir;
		const PS::F64 drv_BperpLx = ep_L.gradBperp_x * ndir;
		const PS::F64 drv_BperpLy = ep_L.gradBperp_y * ndir;

		/**Right i**/
		const PS::F64 drv_densR = ep_R.grad_dens * ndir;
		const PS::F64 drv_PresR = ep_R.gradP * ndir;
		const PS::F64 drv_PTR = ep_R.gradPT * ndir;
		const PS::F64 drv_velparaR = ep_R.gradVpara * ndir;
		const PS::F64 drv_BperpR2 = ep_R.gradBperp2 * ndir;
		const PS::F64 drv_BperpRx = ep_R.gradBperp_x * ndir;
		const PS::F64 drv_BperpRy = ep_R.gradBperp_y * ndir;
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
		const PS::F64 DiffVpara = orig_vparaR - orig_vparaL;
		const PS::F64 DiffBperp2 = orig_BperpR2 - orig_BperpL2;
		const PS::F64 orig_PTL = orig_PresL + .5 * orig_BperpL2;
		const PS::F64 orig_PTR = orig_PresR + .5 * orig_BperpR2;
		const PS::F64 DiffPT = orig_PTR - orig_PTL;
		const PS::F64 DiffPres = orig_PresR - orig_PresL;
		const PS::F64 DiffDens = orig_DensR - orig_DensL;
		/**Left j**/
		const PS::F64 DeltaVparaL = drv_velparaL * deltaS;
		const PS::F64 DeltaBperpL2 = drv_BperpL2 * deltaS;
		const PS::F64 DeltaPTL = drv_PTL * deltaS;
		const PS::F64 DeltaPresL = drv_PresL * deltaS;
		const PS::F64 DeltaDensL = drv_densL * deltaS;
		const PS::F64 DeltaBperpLx = drv_BperpLx * deltaS;
		const PS::F64 DeltaBperpLy = drv_BperpLy * deltaS;

		const PS::F64 DiffVpara_dashL = 2.0 * DeltaVparaL + DiffVpara;
		const PS::F64 DiffBperp_dashL2 = 2.0 * DeltaBperpL2 + DiffBperp2;
		const PS::F64 DiffPT_dashL = 2.0 * DeltaPTL + DiffPT;
		const PS::F64 DiffPres_dashL = 2.0 * DeltaPresL + DiffPres;
		const PS::F64 DiffDens_dashL = 2.0 * DeltaDensL + DiffDens;
		/**Right i**/
		const PS::F64 DeltaVparaR = drv_velparaR * deltaS;
		const PS::F64 DeltaBperpR2 = drv_BperpR2 * deltaS;
		const PS::F64 DeltaPTR = drv_PTR * deltaS;
		const PS::F64 DeltaPresR = drv_PresR * deltaS;
		const PS::F64 DeltaDensR = drv_densR * deltaS;
		const PS::F64 DeltaBperpRx = drv_BperpRx * deltaS;
		const PS::F64 DeltaBperpRy = drv_BperpRy * deltaS;

		const PS::F64 DiffVpara_dashR = 2.0 * DeltaVparaR - DiffVpara;
		const PS::F64 DiffBperp_dashR2 = 2.0 * DeltaBperpR2 - DiffBperp2;
		const PS::F64 DiffPT_dashR = 2.0 * DeltaPTR - DiffPT;
		const PS::F64 DiffPres_dashR = 2.0 * DeltaPresR - DiffPres;
		const PS::F64 DiffDens_dashR = 2.0 * DeltaDensR - DiffDens;
		//////////////////////////Monotone Variables Finish////////////////////////////////////////////////////////////////////////////
		//////////////////////////Monotone Constraints Start////////////////////////////////////////////////////////////////////////////
		const PS::F64 DeltaVparaL_mono = (-sgn(DiffVpara) == sgn(DeltaVparaL) && sgn(DeltaVparaL) == sgn(DiffVpara_dashL)) ? sgn(DeltaVparaL) * trimin(2.0 * fabs(DiffVpara), fabs(DeltaVparaL), 2.0 * fabs(DiffVpara_dashL)) : 0.0;
		const PS::F64 DeltaBperpL2_mono = (-sgn(DiffBperp2) == sgn(DeltaBperpL2) && sgn(DeltaBperpL2) == sgn(DiffBperp_dashL2)) ? sgn(DeltaBperpL2) * trimin(2.0 * fabs(DiffBperp2), fabs(DeltaBperpL2), 2.0 * fabs(DiffBperp_dashL2)) : 0.0;
		const PS::F64 DeltaPTL_mono = (-sgn(DiffPT) == sgn(DeltaPTL) && sgn(DeltaPTL) == sgn(DiffPT_dashL)) ? sgn(DeltaPTL) * trimin(2.0 * fabs(DiffPT), fabs(DeltaPTL), 2.0 * fabs(DiffPT_dashL)) : 0.0;
		const PS::F64 DeltaPresL_mono = (-sgn(DiffPres) == sgn(DeltaPresL) && sgn(DeltaPresL) == sgn(DiffPres_dashL)) ? sgn(DeltaPresL) * trimin(2.0 * fabs(DiffPres), fabs(DeltaPresL), 2.0 * fabs(DiffPres_dashL)) : 0.0;
		const PS::F64 DeltaDensL_mono = (-sgn(DiffDens) == sgn(DeltaDensL) && sgn(DeltaDensL) == sgn(DiffDens_dashL)) ? sgn(DeltaDensL) * trimin(2.0 * fabs(DiffDens), fabs(DeltaDensL), 2.0 * fabs(DiffDens_dashL)) : 0.0;

		const PS::F64 DeltaVparaR_mono = (sgn(DiffVpara) == sgn(DeltaVparaR) && sgn(DeltaVparaR) == sgn(DiffVpara_dashR)) ? sgn(DeltaVparaR) * trimin(2.0 * fabs(DiffVpara), fabs(DeltaVparaR), 2.0 * fabs(DiffVpara_dashR)) : 0.0;
		const PS::F64 DeltaBperpR2_mono = (sgn(DiffBperp2) == sgn(DeltaBperpR2) && sgn(DeltaBperpR2) == sgn(DiffBperp_dashR2)) ? sgn(DeltaBperpR2) * trimin(2.0 * fabs(DiffBperp2), fabs(DeltaBperpR2), 2.0 * fabs(DiffBperp_dashR2)) : 0.0;
		const PS::F64 DeltaPTR_mono = (sgn(DiffPT) == sgn(DeltaPTR) && sgn(DeltaPTR) == sgn(DiffPT_dashR)) ? sgn(DeltaPTR) * trimin(2.0 * fabs(DiffPT), fabs(DeltaPTR), 2.0 * fabs(DiffPT_dashR)) : 0.0;
		const PS::F64 DeltaPresR_mono = (sgn(DiffPres) == sgn(DeltaPresR) && sgn(DeltaPresR) == sgn(DiffPres_dashR)) ? sgn(DeltaPresR) * trimin(2.0 * fabs(DiffPres), fabs(DeltaPresR), 2.0 * fabs(DiffPres_dashR)) : 0.0;
		const PS::F64 DeltaDensR_mono = (sgn(DiffDens) == sgn(DeltaDensR) && sgn(DeltaDensR) == sgn(DiffDens_dashR)) ? sgn(DeltaDensR) * trimin(2.0 * fabs(DiffDens), fabs(DeltaDensR), 2.0 * fabs(DiffDens_dashR)) : 0.0;
//		const PS::F64 DeltaBperpLx_mono = DeltaBperpLx;
//		const PS::F64 DeltaBperpLy_mono = DeltaBperpLy;
//		const PS::F64 DeltaBperpRx_mono = DeltaBperpRx;
//		const PS::F64 DeltaBperpRy_mono = DeltaBperpRy;

//		const PS::F64 DeltaVparaL_mono = DeltaVparaL;
//		const PS::F64 DeltaBperpL2_mono = DeltaBperpL2;
//		const PS::F64 DeltaPTL_mono = DeltaPTL;
//		const PS::F64 DeltaPresL_mono = DeltaPresL;
//		const PS::F64 DeltaDensL_mono = DeltaDensL;
//		const PS::F64 DeltaVparaR_mono = DeltaVparaR;
//		const PS::F64 DeltaBperpR2_mono = DeltaBperpR2;
//		const PS::F64 DeltaPTR_mono = DeltaPTR;
//		const PS::F64 DeltaPresR_mono = DeltaPresR;
//		const PS::F64 DeltaDensR_mono = DeltaDensR;

		//////////////////////////Monotone Constraints Finish////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////RP Variables////////////////////////////////
		const PS::F64 C_fastL = sqrt((gamma * ep_L.pres + orig_MagneticBL * orig_MagneticBL) * orig_idensL);
		const PS::F64 dod_facLRP = .5 * (1.0 - C_fastL * dt / deltaS);
		const PS::F64 C_fastR = sqrt((gamma * ep_R.pres + orig_MagneticBR * orig_MagneticBR) * orig_idensR);
		const PS::F64 dod_facRRP = .5 * (1.0 - C_fastR * dt / deltaS);

		const PS::F64 fac_velparaLRP = DeltaVparaL_mono * dod_facLRP;
		const PS::F64 fac_BperpL2RP = DeltaBperpL2_mono * dod_facLRP;
		const PS::F64 fac_PTLRP = DeltaPTL_mono * dod_facLRP;
		const PS::F64 fac_PresLRP = DeltaPresL_mono * dod_facLRP;
		const PS::F64 fac_densLRP = DeltaDensL_mono * dod_facLRP;
//		const PS::F64 fac_BperpLxRP = DeltaBperpLx_mono * dod_facLRP;
//		const PS::F64 fac_BperpLyRP = DeltaBperpLy_mono * dod_facLRP;

		const PS::F64 fac_velparaRRP = DeltaVparaR_mono * dod_facRRP;
		const PS::F64 fac_BperpR2RP = DeltaBperpR2_mono * dod_facRRP;
		const PS::F64 fac_PTRRP = DeltaPTR_mono * dod_facRRP;
		const PS::F64 fac_PresRRP = DeltaPresR_mono * dod_facRRP;
		const PS::F64 fac_densRRP = DeltaDensR_mono * dod_facRRP;
//		const PS::F64 fac_BperpRxRP = DeltaBperpRx_mono * dod_facRRP;
//		const PS::F64 fac_BperpRyRP = DeltaBperpRy_mono * dod_facRRP;
		////////L -> j///////////////////////////
		const PS::F64 vparaLRP = orig_vparaL + fac_velparaLRP;
		const PS::F64 BperpL2RP = orig_BperpL2 + fac_BperpL2RP;
		const PS::F64 presLRP = orig_PresL + fac_PresLRP;
		const PS::F64 densLRP = orig_DensL + fac_densLRP;
//		const PS::F64 BperpLxRP = orig_BperpL.x + fac_BperpLxRP;
//		const PS::F64 BperpLyRP = orig_BperpL.y + fac_BperpLyRP;
//		const PS::F64 PTLRP = presLRP + .5 * (BperpLxRP * BperpLxRP + BperpLyRP * BperpLyRP);
		const PS::F64 PTLRP = orig_PresL + .5 * orig_BperpL2 + fac_PTLRP;

		//////R -> i///////////////////////////

		const PS::F64 vparaRRP = orig_vparaR - fac_velparaRRP;
		const PS::F64 BperpR2RP = orig_BperpR2 - fac_BperpR2RP;
		const PS::F64 presRRP = orig_PresR - fac_PresRRP;
		const PS::F64 densRRP = orig_DensR - fac_densRRP;
//		const PS::F64 BperpRxRP = orig_BperpR.x - fac_BperpRxRP;
//		const PS::F64 BperpRyRP = orig_BperpR.y - fac_BperpRyRP;
//		const PS::F64 PTRRP = presRRP + .5 * (BperpRxRP * BperpRxRP + BperpRyRP * BperpRyRP);
		const PS::F64 PTRRP = orig_PresR + .5 * orig_BperpR2 - fac_PTRRP;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////L -> j///////////////////////////
//		const PS::F64 vparaLRP = orig_vparaL;
//		const PS::F64 BperpL2RP = orig_BperpL2;
//		const PS::F64 presLRP = orig_PresL;
//		const PS::F64 densLRP = orig_DensL;
//		const PS::F64 PTLRP = orig_PresL + .5 * orig_BperpL2;
//		////////R -> i///////////////////////////
//		const PS::F64 vparaRRP = orig_vparaR;
//		const PS::F64 BperpR2RP = orig_BperpR2;
//		const PS::F64 presRRP = orig_PresR;
//		const PS::F64 densRRP = orig_DensR;
//		const PS::F64 PTRRP = orig_PresR + .5 * orig_BperpR2;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		const PS::F64 ptRRP = PTRRP;
		const PS::F64 ptLRP = PTLRP;
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
			PS::F64 DR2 = pow((gamma + 1) * ptRRP + (gamma - 1) * ptStar, 2.) - (gamma - 2) * BperpR2RP * (2. * (gamma - 3) * ptRRP + 2. * (gamma + 3) * ptStar - (gamma - 2) * BperpR2RP);
			PS::F64 DL2 = pow((gamma + 1) * ptLRP + (gamma - 1) * ptStar, 2.) - (gamma - 2) * BperpL2RP * (2. * (gamma - 3) * ptLRP + 2. * (gamma + 3) * ptStar - (gamma - 2) * BperpL2RP);
			DR2 = DR2 < 0.0 ? 0.0 : DR2;
			DL2 = DL2 < 0.0 ? 0.0 : DL2;
			if (DR2 == 0.0 || DL2 == 0.0) {
				std::cout << "otStar:" << loop << "otStar:" << ptRRP << " " << DR2 << " " << ptLRP << " " << DL2 << " " << sqrt(DL2) << " " << (gamma - 3.0) * ptLRP + (gamma + 3) * ptStar - (gamma - 2) * BperpL2RP << std::endl;

			}
			PS::F64 M_R2 = .25 * densRRP * ((gamma - 3.0) * ptRRP + (gamma + 3) * ptStar - (gamma - 2) * BperpR2RP + sqrt(DR2));
			PS::F64 M_L2 = .25 * densLRP * ((gamma - 3.0) * ptLRP + (gamma + 3) * ptStar - (gamma - 2) * BperpL2RP + sqrt(DL2));
			if (M_R2 < 0.0 || M_L2 < 0.0) {
				std::cout << loop << " otStar:" << ptStar << "otStar:" << ptRRP << " " << M_R2 << " " << ptLRP << " " << M_L2 << " " << sqrt(DL2) << " " << (gamma - 3.0) * ptLRP + (gamma + 3) * ptStar - (gamma - 2) * BperpL2RP << std::endl;

			}
			M_R2 = M_R2 < 0.0 ? 0.0 : M_R2;
			M_L2 = M_L2 < 0.0 ? 0.0 : M_L2;

			if (M_R2 == 0.0 || M_L2 == 0.0) {
				M_L2 = M_R2 = 0.0001;
				std::cout << loop << " otStar:" << ptStar << "otStar:" << ptRRP << " " << M_R2 << " " << ptLRP << " " << M_L2 << " " << sqrt(DL2) << " " << (gamma - 3.0) * ptLRP + (gamma + 3) * ptStar - (gamma - 2) * BperpL2RP << std::endl;

//				ptStar = 0.5 * (ptLRP + ptRRP);
//				V_RP = .5 * (vparaRRP + vparaLRP);
//				Pt_RP = ptStar;
//
//				break;
//			while(1){}
			}
			const PS::F64 M_R = sqrt(M_R2);
			const PS::F64 M_L = sqrt(M_L2);
			ptStar = (ptRRP / M_R + ptLRP / M_L - (vparaRRP - vparaLRP)) / (1. / M_R + 1. / M_L);
			if (ptStar < 0.0) {
//				std::cout << loop << " otStar:" << ptStar << "otStar:" << ptRRP << " " << M_R2 << " " << ptLRP << " " << M_L2 << " " << sqrt(DL2) << " " << (gamma - 3.0) * ptLRP + (gamma + 3) * ptStar - (gamma - 2) * BperpL2RP << std::endl;
				ptStar = 0.5 * (ptLRP + ptRRP);
				V_RP = .5 * (vparaRRP + vparaLRP);
				Pt_RP = ptStar;

				break;
			}
			ptStar = ptStar < 0.0 ? 0.0 : ptStar;
			if (fabs(ptStar - ptStar_old) < 0.01) {
				V_RP = (vparaRRP * M_R + vparaLRP * M_L - (ptRRP - ptLRP)) / (M_R + M_L);
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
		const PS::F64 drv_BperpxL = ep_L.gradBperp_x * ndir;
		const PS::F64 drv_BperpyL = ep_L.gradBperp_y * ndir;
		const PS::F64 drv_vperpxR = ep_R.gradvperp_x * ndir;
		const PS::F64 drv_vperpyR = ep_R.gradvperp_y * ndir;
		const PS::F64 drv_BperpxR = ep_R.gradBperp_x * ndir;
		const PS::F64 drv_BperpyR = ep_R.gradBperp_y * ndir;
//
		const PS::F64vec DiffVperp = orig_vperpR - orig_vperpL;
		const PS::F64vec DiffBperp = orig_BperpR - orig_BperpL;
//
		const PS::F64 DeltaVperpxL = drv_vperpxL * deltaS;
		const PS::F64 DeltaVperpyL = drv_vperpyL * deltaS;
		const PS::F64 DeltaBperpxL = drv_BperpxL * deltaS;
		const PS::F64 DeltaBperpyL = drv_BperpyL * deltaS;
		const PS::F64 DeltaVperpxR = drv_vperpxR * deltaS;
		const PS::F64 DeltaVperpyR = drv_vperpyR * deltaS;
		const PS::F64 DeltaBperpxR = drv_BperpxR * deltaS;
		const PS::F64 DeltaBperpyR = drv_BperpyR * deltaS;
		const PS::F64vec DeltaVperpL(DeltaVperpxL, DeltaVperpyL);
		const PS::F64vec DeltaBperpL(DeltaBperpxL, DeltaBperpyL);
		const PS::F64vec DeltaVperpR(DeltaVperpxR, DeltaVperpyR);
		const PS::F64vec DeltaBperpR(DeltaBperpxR, DeltaBperpyR);
		const PS::F64vec DiffBperp_dashL = 2.0 * DeltaBperpL + DiffBperp;
		const PS::F64vec DiffVperp_dashL = 2.0 * DeltaVperpL + DiffVperp;
		const PS::F64vec DiffBperp_dashR = 2.0 * DeltaBperpR - DiffBperp;
		const PS::F64vec DiffVperp_dashR = 2.0 * DeltaVperpR - DiffVperp;

		const PS::F64 DeltaVperpxL_mono = (-sgn(DiffVperp.x) == sgn(DeltaVperpxL) && sgn(DeltaVperpxL) == sgn(DiffVperp_dashL.x)) ? sgn(DeltaVperpxL) * trimin(2.0 * fabs(DiffVperp.x), fabs(DeltaVperpxL), 2.0 * fabs(DiffVperp_dashL.x)) : 0.0;
		const PS::F64 DeltaVperpyL_mono = (-sgn(DiffVperp.y) == sgn(DeltaVperpyL) && sgn(DeltaVperpyL) == sgn(DiffVperp_dashL.y)) ? sgn(DeltaVperpyL) * trimin(2.0 * fabs(DiffVperp.y), fabs(DeltaVperpyL), 2.0 * fabs(DiffVperp_dashL.y)) : 0.0;
		const PS::F64 DeltaBperpxL_mono = (-sgn(DiffBperp.x) == sgn(DeltaBperpxL) && sgn(DeltaBperpxL) == sgn(DiffBperp_dashL.x)) ? sgn(DeltaBperpxL) * trimin(2.0 * fabs(DiffBperp.x), fabs(DeltaBperpxL), 2.0 * fabs(DiffBperp_dashL.x)) : 0.0;
		const PS::F64 DeltaBperpyL_mono = (-sgn(DiffBperp.y) == sgn(DeltaBperpyL) && sgn(DeltaBperpyL) == sgn(DiffBperp_dashL.y)) ? sgn(DeltaBperpyL) * trimin(2.0 * fabs(DiffBperp.y), fabs(DeltaBperpyL), 2.0 * fabs(DiffBperp_dashL.y)) : 0.0;

		const PS::F64 DeltaVperpxR_mono = (sgn(DiffVperp.x) == sgn(DeltaVperpxR) && sgn(DeltaVperpxR) == sgn(DiffVperp_dashR.x)) ? sgn(DeltaVperpxR) * trimin(2.0 * fabs(DiffVperp.x), fabs(DeltaVperpxR), 2.0 * fabs(DiffVperp_dashR.x)) : 0.0;
		const PS::F64 DeltaVperpyR_mono = (sgn(DiffVperp.y) == sgn(DeltaVperpyR) && sgn(DeltaVperpyR) == sgn(DiffVperp_dashR.y)) ? sgn(DeltaVperpyR) * trimin(2.0 * fabs(DiffVperp.y), fabs(DeltaVperpyR), 2.0 * fabs(DiffVperp_dashR.y)) : 0.0;
		const PS::F64 DeltaBperpxR_mono = (sgn(DiffBperp.x) == sgn(DeltaBperpxR) && sgn(DeltaBperpxR) == sgn(DiffBperp_dashR.x)) ? sgn(DeltaBperpxR) * trimin(2.0 * fabs(DiffBperp.x), fabs(DeltaBperpxR), 2.0 * fabs(DiffBperp_dashR.x)) : 0.0;
		const PS::F64 DeltaBperpyR_mono = (sgn(DiffBperp.y) == sgn(DeltaBperpyR) && sgn(DeltaBperpyR) == sgn(DiffBperp_dashR.y)) ? sgn(DeltaBperpyR) * trimin(2.0 * fabs(DiffBperp.y), fabs(DeltaBperpyR), 2.0 * fabs(DiffBperp_dashR.y)) : 0.0;

//		const PS::F64 DeltaVperpxL_mono = DeltaVperpxL;
//		const PS::F64 DeltaVperpyL_mono = DeltaVperpyL;
//		const PS::F64 DeltaBperpxL_mono = DeltaBperpxL;
//		const PS::F64 DeltaBperpyL_mono = DeltaBperpyL;
//		const PS::F64 DeltaVperpxR_mono = DeltaVperpxR;
//		const PS::F64 DeltaVperpyR_mono = DeltaVperpyR;
//		const PS::F64 DeltaBperpxR_mono = DeltaBperpxR;
//		const PS::F64 DeltaBperpyR_mono = DeltaBperpyR;

		//ToDO Using Riemann Invariant
//		const PS::F64vec orig_JPlusL = orig_vperpL - orig_BperpL * orig_idensLsqrt;
//		const PS::F64vec orig_JMinusL = orig_vperpL + orig_BperpL * orig_idensLsqrt;
//		const PS::F64vec orig_JPlusR = orig_vperpR - orig_BperpR * orig_idensRsqrt;
//		const PS::F64vec orig_JMinusR = orig_vperpR + orig_BperpR * orig_idensRsqrt;
//
//		const PS::F64 drv_JPlusxL = drv_vperpxL - drv_BperpxL * orig_idensLsqrt;
//		const PS::F64 drv_JPlusyL = drv_vperpyL - drv_BperpyL * orig_idensLsqrt;
//		const PS::F64 drv_JMinusxL = drv_vperpxL + drv_BperpxL * orig_idensLsqrt;
//		const PS::F64 drv_JMinusyL = drv_vperpyL + drv_BperpyL * orig_idensLsqrt;
//
//		const PS::F64 drv_JPlusxR = drv_vperpxR - drv_BperpxR * orig_idensRsqrt;
//		const PS::F64 drv_JPlusyR = drv_vperpyR - drv_BperpyR * orig_idensRsqrt;
//		const PS::F64 drv_JMinusxR = drv_vperpxR + drv_BperpxR * orig_idensRsqrt;
//		const PS::F64 drv_JMinusyR = drv_vperpyR + drv_BperpyR * orig_idensRsqrt;
//
//		const PS::F64vec DiffJPlus = orig_JPlusR - orig_JPlusL;
//		const PS::F64vec DiffJMinus = orig_JMinusR - orig_JMinusL;
//		const PS::F64 DeltaJPlusxL = drv_JPlusxL * deltaS;
//		const PS::F64 DeltaJPlusyL = drv_JPlusyL * deltaS;
//		const PS::F64 DeltaJPlusxR = drv_JPlusxR * deltaS;
//		const PS::F64 DeltaJPlusyR = drv_JPlusyR * deltaS;
//		const PS::F64 DeltaJMinusxL = drv_JMinusxL * deltaS;
//		const PS::F64 DeltaJMinusyL = drv_JMinusyL * deltaS;
//		const PS::F64 DeltaJMinusxR = drv_JMinusxR * deltaS;
//		const PS::F64 DeltaJMinusyR = drv_JMinusyR * deltaS;
//		const PS::F64vec DeltaJPlusL(DeltaJPlusxL, DeltaJPlusyL);
//		const PS::F64vec DeltaJPlusR(DeltaJPlusxR, DeltaJPlusyR);
//		const PS::F64vec DeltaJMinusL(DeltaJMinusxL, DeltaJMinusyL);
//		const PS::F64vec DeltaJMinusR(DeltaJMinusxR, DeltaJMinusyR);
//		const PS::F64vec DiffJPlus_dashL = 2.0 * DeltaJPlusL + DiffJPlus;
//		const PS::F64vec DiffJPlus_dashR = 2.0 * DeltaJPlusR - DiffJPlus;
//		const PS::F64vec DiffJMinus_dashL = 2.0 * DeltaJMinusL + DiffJMinus;
//		const PS::F64vec DiffJMinus_dashR = 2.0 * DeltaJMinusR - DiffJMinus;
//		const PS::F64 DeltaJPlusxL_mono = (-sgn(DiffJPlus.x) == sgn(DeltaJPlusxL) && sgn(DeltaJPlusxL) == sgn(DiffJPlus_dashL.x)) ? sgn(DeltaJPlusxL) * trimin(2.0 * fabs(DiffJPlus_dashL.x), fabs(DeltaJPlusxL), 2.0 * fabs(DiffJPlus_dashL.x)) : 0.0;
//		const PS::F64 DeltaJPlusyL_mono = (-sgn(DiffJPlus.y) == sgn(DeltaJPlusyL) && sgn(DeltaJPlusyL) == sgn(DiffJPlus_dashL.y)) ? sgn(DeltaJPlusyL) * trimin(2.0 * fabs(DiffJPlus_dashL.y), fabs(DeltaJPlusyL), 2.0 * fabs(DiffJPlus_dashL.y)) : 0.0;
//		const PS::F64 DeltaJMinusxL_mono = (-sgn(DiffJMinus.x) == sgn(DeltaJMinusxL) && sgn(DeltaJMinusxL) == sgn(DiffJMinus_dashL.x)) ? sgn(DeltaJMinusxL) * trimin(2.0 * fabs(DiffJMinus_dashL.x), fabs(DeltaJMinusxL), 2.0 * fabs(DiffJMinus_dashL.x)) : 0.0;
//		const PS::F64 DeltaJMinusyL_mono = (-sgn(DiffJMinus.y) == sgn(DeltaJMinusyL) && sgn(DeltaJMinusyL) == sgn(DiffJMinus_dashL.y)) ? sgn(DeltaJMinusyL) * trimin(2.0 * fabs(DiffJMinus_dashL.y), fabs(DeltaJMinusyL), 2.0 * fabs(DiffJMinus_dashL.y)) : 0.0;
//		const PS::F64 DeltaJPlusxR_mono = (sgn(DiffJPlus.x) == sgn(DeltaJPlusxR) && sgn(DeltaJPlusxR) == sgn(DiffJPlus_dashR.x)) ? sgn(DeltaJPlusxR) * trimin(2.0 * fabs(DiffJPlus_dashR.x), fabs(DeltaJPlusxR), 2.0 * fabs(DiffJPlus_dashR.x)) : 0.0;
//		const PS::F64 DeltaJPlusyR_mono = (sgn(DiffJPlus.y) == sgn(DeltaJPlusyR) && sgn(DeltaJPlusyR) == sgn(DiffJPlus_dashR.y)) ? sgn(DeltaJPlusyR) * trimin(2.0 * fabs(DiffJPlus_dashR.y), fabs(DeltaJPlusyR), 2.0 * fabs(DiffJPlus_dashR.y)) : 0.0;
//		const PS::F64 DeltaJMinusxR_mono = (sgn(DiffJMinus.x) == sgn(DeltaJMinusxR) && sgn(DeltaJMinusxR) == sgn(DiffJMinus_dashR.x)) ? sgn(DeltaJMinusxR) * trimin(2.0 * fabs(DiffJMinus_dashR.x), fabs(DeltaJMinusxR), 2.0 * fabs(DiffJMinus_dashR.x)) : 0.0;
//		const PS::F64 DeltaJMinusyR_mono = (sgn(DiffJMinus.y) == sgn(DeltaJMinusyR) && sgn(DeltaJMinusyR) == sgn(DiffJMinus_dashR.y)) ? sgn(DeltaJMinusyR) * trimin(2.0 * fabs(DiffJMinus_dashR.y), fabs(DeltaJMinusyR), 2.0 * fabs(DiffJMinus_dashR.y)) : 0.0;
//		const PS::F64 DeltaVperpxL_mono = .5 * (DeltaJPlusxL_mono + DeltaJMinusxL_mono);
//		const PS::F64 DeltaVperpyL_mono = .5 * (DeltaJPlusyL_mono + DeltaJMinusyL_mono);
//		const PS::F64 DeltaBperpxL_mono = .5 * (DeltaJMinusxL_mono - DeltaJPlusxL_mono) * orig_densLsqrt;
//		const PS::F64 DeltaBperpyL_mono = .5 * (DeltaJMinusyL_mono - DeltaJPlusyL_mono) * orig_densLsqrt;
//		const PS::F64 DeltaVperpxR_mono = .5 * (DeltaJPlusxR_mono + DeltaJMinusxR_mono);
//		const PS::F64 DeltaVperpyR_mono = .5 * (DeltaJPlusyR_mono + DeltaJMinusyR_mono);
//		const PS::F64 DeltaBperpxR_mono = .5 * (DeltaJMinusxR_mono - DeltaJPlusxR_mono) * orig_densRsqrt;
//		const PS::F64 DeltaBperpyR_mono = .5 * (DeltaJMinusyR_mono - DeltaJPlusyR_mono) * orig_densRsqrt;

		////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////
		const PS::F64 C_alvenL = fabs(orig_BparaL) * orig_idensLsqrt;
		const PS::F64 C_alvenR = fabs(orig_BparaR) * orig_idensRsqrt;
		const PS::F64 dod_fac_AlfvenL = .5 * (1.0 - C_alvenL * dt / deltaS);
		const PS::F64 dod_fac_AlfvenR = .5 * (1.0 - C_alvenR * dt / deltaS);
		const PS::F64 fac_densLMoC = DeltaDensL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_vperpxLMoC = DeltaVperpxL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_vperpyLMoC = DeltaVperpyL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_BperpxLMoC = DeltaBperpxL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_BperpyLMoC = DeltaBperpyL_mono * dod_fac_AlfvenL;
		const PS::F64 fac_densRMoC = DeltaDensR_mono * dod_fac_AlfvenR;
		const PS::F64 fac_vperpxRMoC = DeltaVperpxR_mono * dod_fac_AlfvenR;
		const PS::F64 fac_vperpyRMoC = DeltaVperpyR_mono * dod_fac_AlfvenR;
		const PS::F64 fac_BperpxRMoC = DeltaBperpxR_mono * dod_fac_AlfvenR;
		const PS::F64 fac_BperpyRMoC = DeltaBperpyR_mono * dod_fac_AlfvenR;

		const PS::F64 densL_guardMoC = orig_DensL;
		const PS::F64 densR_guardMoC = orig_DensR;
//		const PS::F64 densLMoC = densL_guardMoC < 0.0 ? orig_DensL : densL_guardMoC;
//		const PS::F64 densRMoC = densR_guardMoC < 0.0 ? orig_DensR : densR_guardMoC;
		const PS::F64 densLMoC = orig_DensL;
		const PS::F64 densRMoC = orig_DensR;
		const PS::F64 vperpxLMoC = orig_vperpL.x;
		const PS::F64 vperpyLMoC = orig_vperpL.y;
		const PS::F64 vperpxRMoC = orig_vperpR.x;
		const PS::F64 vperpyRMoC = orig_vperpR.y;
		const PS::F64 BperpxLMoC = orig_BperpL.x;
		const PS::F64 BperpyLMoC = orig_BperpL.y;
		const PS::F64 BperpxRMoC = orig_BperpR.x;
		const PS::F64 BperpyRMoC = orig_BperpR.y;

//		const PS::F64 densL_guardMoC = orig_DensL - fac_densLMoC;
//		const PS::F64 densR_guardMoC = orig_DensR + fac_densRMoC;
//		const PS::F64 densLMoC = densL_guardMoC < 0.0 ? orig_DensL : densL_guardMoC;
//		const PS::F64 densRMoC = densR_guardMoC < 0.0 ? orig_DensR : densR_guardMoC;
//		const PS::F64 vperpxLMoC = orig_vperpL.x - fac_vperpxLMoC;
//		const PS::F64 vperpyLMoC = orig_vperpL.y - fac_vperpyLMoC;
//		const PS::F64 BperpxLMoC = orig_BperpL.x - fac_BperpxLMoC;
//		const PS::F64 BperpyLMoC = orig_BperpL.y - fac_BperpyLMoC;
//		const PS::F64 vperpxRMoC = orig_vperpR.x + fac_vperpxRMoC;
//		const PS::F64 vperpyRMoC = orig_vperpR.y + fac_vperpyRMoC;
//		const PS::F64 BperpxRMoC = orig_BperpR.x + fac_BperpxRMoC;
//		const PS::F64 BperpyRMoC = orig_BperpR.y + fac_BperpyRMoC;

		const PS::F64vec vperpLMoC(vperpxLMoC, vperpyLMoC);
		const PS::F64vec BperpLMoC(BperpxLMoC, BperpyLMoC);
		const PS::F64vec vperpRMoC(vperpxRMoC, vperpyRMoC);
		const PS::F64vec BperpRMoC(BperpxRMoC, BperpyRMoC);

		const PS::F64 idensLMoC = 1. / densLMoC;
		const PS::F64 densLMoCsqrt = sqrt(densLMoC);
		const PS::F64 idensLMoCsqrt = 1. / sqrt(densLMoC);
		const PS::F64 idensRMoC = 1. / densRMoC;
		const PS::F64 densRMoCsqrt = sqrt(densRMoC);
		const PS::F64 idensRMoCsqrt = 1. / sqrt(densRMoC);

		const PS::F64 _Bpara = .5 * (orig_BparaL + orig_BparaR);
		BperpMoC = (BperpRMoC * idensRMoCsqrt + BperpLMoC * idensLMoCsqrt + sgn(_Bpara) * (-vperpLMoC + vperpRMoC)) / (idensLMoCsqrt + idensRMoCsqrt);
		VperpMoC = (densRMoCsqrt * vperpRMoC + densLMoCsqrt * vperpLMoC + sgn(_Bpara) * (-BperpLMoC + BperpRMoC)) / (densLMoCsqrt + densRMoCsqrt);
//		BperpMoC = .5 * (BperpRMoC + BperpLMoC);

//		if(BperpStar!=BperpStar){
//							std::cout<<idensRMoCsqrt<<" "<<densRMoC<<" "<<idensLMoCsqrt<<" "<<densLMoC<<std::endl;
//						}
		BparaStar = _Bpara;
	}

};

template<class TParticleJ> class CalcGravityForce {
public:
	void operator ()(const EPI::Grav* const ep_i, const PS::S32 Nip, const TParticleJ* const ep_j, const PS::S32 Njp, RESULT::Grav* const grav) {
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

