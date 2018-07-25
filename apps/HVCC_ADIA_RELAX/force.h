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
				if (dens[i].dens != dens[i].dens) {
					std::cout << dens[i].dens << " --- " << ith.smth
							<< std::endl;
				}
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
				drvt[i].gradvel_x += jth.mass * (vperpj.x - vperpi.x)
						* kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradvel_y += jth.mass * (vperpj.y - vperpi.y)
						* kernel.gradW(dr, ith.smth) / jth.dens;
				drvt[i].gradvel_z += jth.mass * (vperpj.z - vperpi.z)
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

			}

		}
	}
};

class CalcHydroForce {
	const kernel_t kernel;
	const PS::S32 NSIZE = 2;
	const PS::S32 mesh_size = 256;
	const PS::F64 xi = 3.65375374;
	const PS::F64 dxi = -0.2033013;
	const PS::F64 ddxi = -xi * xi * dxi;
public:
	void calc_Vij2_and_ss(const EPI::Hydro ep_i, const EPJ::Hydro ep_j,
			PS::F64&Vij2_hi, PS::F64&Vij2_hj, PS::F64& ssP, const PS::F64 delta,
			const PS::F64vec eij) {
		//		const PS::F64 dVi = ep_i.gradV * eij;
		//		const PS::F64 dVj = ep_j.gradV * eij;
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
		//		Aij = -2 * (Vi - Vj) * pow(delta, -3) + (dVi + dVj) * pow(delta, -2);
		//		Bij = 0.5 * (dVi - dVj) / delta;
		//		Cij = 1.5 * (Vi - Vj) / delta - 0.25 * (dVi + dVj);
		//		Dij = 0.5 * (Vi + Vj) - 0.125 * (dVi - dVj) * delta;
		//		Vij2_hi = ((0.234375 * hi2 * Aij * Aij + 0.1875 * (2 * Aij * Cij + Bij * Bij)) * hi2 + 0.25 * (2 * Bij * Dij + Cij * Cij)) * hi2 + Dij * Dij;
		//		Vij2_hj = ((0.234375 * hj2 * Aij * Aij + 0.1875 * (2 * Aij * Cij + Bij * Bij)) * hj2 + 0.25 * (2 * Bij * Dij + Cij * Cij)) * hj2 + Dij * Dij;
		//		ssP = 0.5 * (((0.46875 * hi2 * Aij * Bij + 0.375 * (Aij * Dij + Bij * Cij)) * hi2 + 0.5 * Cij * Dij) * hi2 / Vij2_hi + ((0.46875 * hj2 * Aij * Bij + 0.375 * (Aij * Dij + Bij * Cij)) * hj2 + 0.5 * Cij * Dij) * hj2 / Vij2_hj);
	}
	void calc_riemann_solver(const EPI::Hydro ep_R, const EPJ::Hydro ep_L,
			const PS::F64 ss, const PS::F64 delta, const PS::F64vec eij,
			const PS::F64 dt, PS::F64 &pstar, PS::F64 &vstar) {
		PS::F64 si = .5 * delta;
		PS::F64 sj = .5 * delta;
		PS::F64 drdsR = ep_R.grad_dens * eij;
		PS::F64 dpdsR = ep_R.gradP * eij;
		PS::F64 vR = ep_R.vel * eij;
		PS::F64 drdsL = ep_L.grad_dens * eij;
		PS::F64 dpdsL = ep_L.gradP * eij;
		PS::F64 vL = ep_L.vel * eij;

		PS::F64 dvdsR = ep_R.gradvel_x * eij + ep_R.gradvel_y * eij
				+ ep_R.gradvel_z * eij;
		PS::F64 dvdsL = ep_L.gradvel_x * eij + ep_L.gradvel_y * eij
				+ ep_L.gradvel_z * eij;
		PS::F64 ddensi = drdsR * (-ep_R.snds * 0.5 * dt + si);
		PS::F64 dpresi = dpdsR * (-ep_R.snds * 0.5 * dt + si);
		PS::F64 dveli = dvdsR * (-ep_R.snds * 0.5 * dt + si);
		PS::F64 ddensj = drdsL * (-ep_L.snds * 0.5 * dt + sj);
		PS::F64 dpresj = dpdsL * (-ep_L.snds * 0.5 * dt + sj);
		PS::F64 dvelj = dvdsL * (-ep_L.snds * 0.5 * dt + sj);

		PS::F64 _rhoR = ep_R.dens - ddensi;
		PS::F64 _pR = ep_R.pres - dpresi;
		PS::F64 _vR = vR - dveli;
		PS::F64 _rhoL = ep_L.dens + ddensj;
		PS::F64 _pL = ep_L.pres + dpresj;
		PS::F64 _vL = vL + dvelj;

		PS::F64 W1, W2;
		const PS::F64 alpha = (2.0 * PARAM::GAMMA) / (PARAM::GAMMA - 1.0);
		double critif = 1.0;
		double pRs = _pR;
		+.5 * eij * ep_R.grav * sqrt(PARAM::GAMMA * _pR * _rhoR) * dt;
		double pLs = _pL;
		-.5 * eij * ep_L.grav * sqrt(PARAM::GAMMA * _pL * _rhoL) * dt;
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////

		if (3.0 * fabs(vL - vR)
				> ((ep_R.snds < ep_L.snds) ? ep_R.snds : ep_L.snds) ||_pR<0||_pL<0||_rhoR<0||_rhoL<0) {

			PS::F64 _rhoR = ep_R.dens;
			PS::F64 _pR = ep_R.pres;
			PS::F64 _vR = vR;
			PS::F64 _rhoL = ep_L.dens;
			PS::F64 _pL = ep_L.pres;
			PS::F64 _vL = vL;

			PS::F64 W1, W2;
			const PS::F64 alpha = (2.0 * PARAM::GAMMA) / (PARAM::GAMMA - 1.0);

			double critif = 1.0;
			double pRs = _pR;
			+.5 * eij * ep_R.grav * sqrt(PARAM::GAMMA * _pR * _rhoR) * dt;
			double pLs = _pL;
			-.5 * eij * ep_L.grav * sqrt(PARAM::GAMMA * _pL * _rhoL) * dt;
			PS::F64 p = .5 * (pRs + pLs);
			pstar = p;
			PS::F64 ppre = p;
			vstar = .5 * (_vR + _vL);
			for (PS::U32 loop = 0; loop < 10; loop++) {
				ppre = p;
				if ((fabs(p / pRs) + 1.0e-9 < critif)) {
					W1 =
							sqrt(pRs * _rhoR)
									* ((PARAM::GAMMA - 1.0)
											/ (2.0 * sqrt(PARAM::GAMMA)))
									* (1.0 - p / pRs)
									/ (1.0 - pow(p / pRs, 1.0 / alpha));
				} else {
					W1 = sqrt(pRs * _rhoR)
							* sqrt(
									0.5 * (PARAM::GAMMA + 1.0) * p / (pRs)
											+ 0.5 * (PARAM::GAMMA - 1.0));
				}
				if ((fabs(p / pLs) + 1.0e-9 < critif)) {
					W2 =
							sqrt(pLs * _rhoL)
									* ((PARAM::GAMMA - 1.0)
											/ (2.0 * sqrt(PARAM::GAMMA)))
									* (1.0 - p / pLs)
									/ (1.0 - pow(p / pLs, 1.0 / alpha));
				} else {
					W2 = sqrt(pLs * _rhoL)
							* sqrt(
									0.5 * (PARAM::GAMMA + 1.0) * p / (pLs)
											+ 0.5 * (PARAM::GAMMA - 1.0));
				}
				p = ((pLs / W2 + pRs / W1) + _vL - _vR) / (1.0 / W2 + 1.0 / W1);
				if (p < 0.0) {
					p = 0.5 * ppre;
					break;
				} else if (1.e-12 + fabs(1.0 - p / ppre) < 1.e-3) {
					break;
				}
			}
			vstar = ((W1 * _vR + W2 * _vL) + pLs - pRs) / (W1 + W2);
			pstar = p;
			if ((pstar != pstar) || vstar != vstar || pstar <= 0.0) {
				vstar = .5 * (vR + vL);
				pstar = .5 * (_pR + _pL);
				std::cout << _pR << "/// " << _pL << "/// " << pRs << "/// "
						<< "/// " << pLs << "/// " << W1 << "/// " << W2
						<< std::endl;
			}
		} else {
			PS::F64 p = .5 * (pRs + pLs);
			pstar = p;
			PS::F64 ppre = p;
			vstar = .5 * (_vR + _vL);
			for (PS::U32 loop = 0; loop < 10; loop++) {
				ppre = p;
				if ((fabs(p / pRs) + 1.0e-9 < critif)) {
					W1 =
							sqrt(pRs * _rhoR)
									* ((PARAM::GAMMA - 1.0)
											/ (2.0 * sqrt(PARAM::GAMMA)))
									* (1.0 - p / pRs)
									/ (1.0 - pow(p / pRs, 1.0 / alpha));
				} else {
					W1 = sqrt(pRs * _rhoR)
							* sqrt(
									0.5 * (PARAM::GAMMA + 1.0) * p / (pRs)
											+ 0.5 * (PARAM::GAMMA - 1.0));
				}
				if ((fabs(p / pLs) + 1.0e-9 < critif)) {
					W2 =
							sqrt(pLs * _rhoL)
									* ((PARAM::GAMMA - 1.0)
											/ (2.0 * sqrt(PARAM::GAMMA)))
									* (1.0 - p / pLs)
									/ (1.0 - pow(p / pLs, 1.0 / alpha));
				} else {
					W2 = sqrt(pLs * _rhoL)
							* sqrt(
									0.5 * (PARAM::GAMMA + 1.0) * p / (pLs)
											+ 0.5 * (PARAM::GAMMA - 1.0));
				}
				p = ((pLs / W2 + pRs / W1) + _vL - _vR) / (1.0 / W2 + 1.0 / W1);
				if (p < 0.0) {
					p = 0.5 * ppre;
					break;
				} else if (1.e-12 + fabs(1.0 - p / ppre) < 1.e-3) {
					break;
				}
			}
			pstar = p;
			if ((pstar != pstar) || vstar != vstar || pstar <= 0.0) {
				vstar = .5 * (vR + vL);
				pstar = .5 * (_pR + _pL);
				std::cout <<dt<<"^^^"<< _pR << "in secound " << _pL << "/// " << pRs
						<< "/// " << "/// " << pLs << "/// " << W1 << "/// "
						<< W2 << std::endl;
			}
			vstar = .5 * (vR + vL + dveli + dvelj);
		}

	}
	double getPhi(double x) {
		double a = 9.999995857999977E-01;
		double b = -3.215355146247138E-04;
		double c = -1.667864474453594E-01;
		double d = 5.252845006191754E-04;
		double e = 1.155630978128144E-02;
		double f = 1.050795977407310E-03;
		double g = -1.389661348847463E-03;
		double h = 2.648078982351008E-04;
		double i = -1.692465157938468E-05;
		double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3)
				+ e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7)
				+ i * pow(x, 8);

		return phi;

	}
	double getPhi_dash(double x) {
		double a = -3.257098683322130E-04;
		double b = -3.335390631705378E-01;
		double c = 1.494801635423298E-03;
		double d = 4.629812975807299E-02;
		double e = 5.242083986884144E-03;
		double f = -8.358773124387748E-03;
		double g = 1.867615193618750E-03;
		double h = -1.388085203319444E-04;
		double i = 2.998565502661576E-07;
		double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3)
				+ e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7)
				+ i * pow(x, 8);

		return phi;

	}
	void operator ()(const EPI::Hydro* const ep_i, const PS::S32 Nip,
			const EPJ::Hydro* const ep_j, const PS::S32 Njp,
			RESULT::Hydro* const hydro) {

		for (PS::S32 i = 0; i < Nip; ++i) {
			hydro[i].clear();
			const EPI::Hydro& ith = ep_i[i];

			for (PS::S32 j = 0; j < Njp; ++j) {
				const EPJ::Hydro& jth = ep_j[j];

				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64 delta = sqrt(dr * dr);
				if (dr == 0.0) {
					continue;
				}
				PS::F64 PSTAR, VSTAR;
				PS::F64 sstar = 0.0;
				PS::F64vec interpo = 0.0;
				const PS::F64vec dv = ith.vel - jth.vel;

				const PS::F64 dr_norm = sqrt(dr * dr);
				const PS::F64vec ndir = dr / dr_norm;

				const PS::F64 s_ij = dr_norm;
				const PS::F64 hbar_ij = .5 * (ith.smth + jth.smth);
				const PS::F64 F_ij = (kernel.gradW_scoord(s_ij, hbar_ij)
						/ (ith.dens * ith.dens)
						+ kernel.gradW_scoord(s_ij, hbar_ij)
								/ (jth.dens * jth.dens));

				const PS::F64vec gradW_hi = kernel.gradW(dr, ith.smth);
				const PS::F64vec gradW_hj = kernel.gradW(dr,
						sqrt(2.0) * jth.smth);

				calc_riemann_solver(ith, jth, sstar, delta, ndir, ith.dt, PSTAR,
						VSTAR);
				const PS::F64vec fac_acc = -PSTAR * F_ij * ndir;
//				const PS::F64 fac_eng_dot = -PSTAR * (VSTAR - ith.vel_half * ndir) * F_ij;
				PS::F64 velref = 0.0;
//				velref = (jth.vel - ith.vel) * ndir;
//				if ((ith.vel - jth.vel) * ndir < 0.0) {
//					velref = (jth.vel - ith.vel) * ndir;
//				} else {
				velref = (VSTAR - ith.vel_half * ndir);
//				}
				const PS::F64 fac_eng_dot = -PSTAR * velref * F_ij;
				hydro[i].acc += jth.mass * fac_acc;
				hydro[i].eng_dot += jth.mass * fac_eng_dot;
//				if (hydro[i].acc.x != hydro[i].acc.x || hydro[i].eng_dot != hydro[i].eng_dot) {
//					std::cout << hydro[i].acc << std::endl;
//				}
//				std::cout << hydro[i].acc << std::endl;

			}
			const PS::F64 _r = sqrt(ep_i[i].pos * ep_i[i].pos);
			double current_phi_dash = getPhi_dash(_r);
			PS::F64vec acc_factor = ep_i[i].pos / _r;

			hydro[i].acc += current_phi_dash * acc_factor;
//			PS::F64vec bh_loc(0, 0, 0);
//
//			const PS::F64 _r2 = (ep_i[i].pos - bh_loc) * (ep_i[i].pos - bh_loc);
//			const PS::F64 _r = sqrt(_r2);
//
//			const PS::F64 _r3 = _r * _r * _r;
//
//			PS::F64 xtilde = _r / ith.smth;
//			PS::F64 xtilde2 = xtilde * xtilde;
//			PS::F64 formfac = 0.0;
//
//			formfac = erf(xtilde) - exp(-xtilde2) * xtilde * 2 / sqrt(M_PI);
//
//			//formfac=1.0;
//			PS::F64vec BHGrav = (ep_i[i].pos - bh_loc) * formfac * PARAM::G * PARAM::sM_BH / (_r3);
//			hydro[i].extF = -BHGrav;
			hydro[i].dt = 0.001 * ith.smth / ith.snds;
			hydro[i].dt = fmin(hydro[i].dt,
					0.02 * sqrt(ith.smth / sqrt(ith.grav * ith.grav)));
//			hydro[i].dt = fmin(hydro[i].dt, 0.02 * sqrt(ith.smth / sqrt(BHGrav * BHGrav)));

		}
	}

}
;

template<class TParticleJ> class CalcGravityForce {
public:
	void operator ()(const EPI::Grav* const ep_i, const PS::S32 Nip,
			const TParticleJ* const ep_j, const PS::S32 Njp,
			RESULT::Grav* const grav) {
		for (PS::S32 i = 0; i < Nip; ++i) {
			const EPI::Grav& ith = ep_i[i];
			for (PS::S32 j = 0; j < Njp; ++j) {
				const TParticleJ& jth = ep_j[j];
				if (ith.pos != jth.pos) {
					const PS::F64vec dr = ith.pos - jth.pos;
					const PS::F64 dr2 = dr * dr;

					const PS::F64 rij = sqrt(dr * dr);
					const PS::F64 dr_inv = 1.0 / sqrt(dr2);

					if (typeid(TParticleJ) == typeid(EPJ::Grav)) {

						const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
						PS::F64 hmax = fmax(ith.smth, jth.smth);
						PS::F64 hmin = fmin(ith.smth, jth.smth);
						PS::F64 x = hmin / hmax;
						PS::F64 fit = 1.0 / (1.0 + x * x);
						PS::F64 xtilde = rij * fit / hmax;
						PS::F64 xtilde2 = xtilde * xtilde;

						PS::F64 formfac = 0.0;

						formfac = erf(xtilde)
								- exp(-xtilde2) * xtilde * 2 / sqrt(M_PI);
						grav[i].acc += -PARAM::G * m_dr3_inv * dr * formfac;
					} else if (typeid(TParticleJ) == typeid(PS::SPJMonopole)) {
						//				std::cout<<"in"<<std::endl;
						const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
						PS::F64 formfac = 1.0;
						grav[i].acc += -PARAM::G * m_dr3_inv * dr * formfac;
						grav[i].pot += -PARAM::G * jth.mass * dr_inv;
					}
				}
			}
		}
	}
};

