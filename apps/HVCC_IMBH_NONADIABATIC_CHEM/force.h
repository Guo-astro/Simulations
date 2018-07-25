#pragma once
class CalcDensity {
	kernel_t kernel;
	double hmin = 1e30;
public:
	void operator ()(const EPI::Dens* const ep_i, const PS::S32 Nip,
			const EPJ::Dens* const ep_j, const PS::S32 Njp,
			RESULT::Dens* const dens) {

		for (PS::S32 i = 0; i < Nip; ++i) {
			dens[i].clear();
			double tmp_dens = 0.0;
			const EPI::Dens& ith = ep_i[i];

			for (PS::S32 j = 0; j < Njp; ++j) {
				const EPJ::Dens& jth = ep_j[j];
				const PS::F64vec dr = jth.pos - ith.pos;
				dens[i].dens += jth.mass * kernel.W(dr, ith.smth);
				tmp_dens += jth.mass * kernel.W(dr, ith.smth);

			}
			dens[i].smth = PARAM::SMTH
					* pow(ith.mass / tmp_dens, 1.0 / (PS::F64)(PARAM::Dim));

		}

	}
};

class CalcDerivative {
	kernel_t kernel;
public:
	void operator ()(const EPI::Drvt* const ep_i, const PS::S32 Nip,
			const EPJ::Drvt* const ep_j, const PS::S32 Njp,
			RESULT::Drvt* const gradients) {
		for (PS::S32 i = 0; i < Nip; ++i) {
			gradients[i].clear();
			const PS::F64 rhoi_inv2 = 1.0 / (ep_i[i].dens * ep_i[i].dens);
			for (PS::S32 j = 0; j < Njp; ++j) {

				const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;
				const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;

				const PS::F64vec gradW = kernel.gradW(dr, ep_i[i].smth);
				gradients[i].gradV += -ep_j[j].mass * gradW;
				gradients[i].graddens += ep_j[j].mass * gradW;
				gradients[i].gradpres += ep_j[j].mass
						* (ep_j[j].pres - ep_i[i].pres) * gradW / ep_j[j].dens;
				gradients[i].gradvel_x += (ep_j[j].vel.x - ep_i[i].vel.x)
						* gradW * ep_j[j].mass / ep_j[j].dens;
				gradients[i].gradvel_y += (ep_j[j].vel.y - ep_i[i].vel.y)
						* gradW * ep_j[j].mass / ep_j[j].dens;
				gradients[i].gradvel_z += (ep_j[j].vel.z - ep_i[i].vel.z)
						* gradW * ep_j[j].mass / ep_j[j].dens;

			}
			gradients[i].gradV *= rhoi_inv2;

		}
	}
};

class CalcHydroForce {
	const kernel_t kernel;
	const PS::S32 NSIZE = 2;
	const PS::S32 mesh_size = 256;
	const PS::F64 xi = 3.65375374;
	const PS::F64 dxi = -0.2033013;
	//const PS::F64 xi = pi;
	//const PS::F64 dxi =-1/pi;
	//const PS::F64 ddxi = 2.71406;
	const PS::F64 ddxi = -xi * xi * dxi;
public:
	double f1(double t, double x, double v) {
		return v;
	}

	double f2(double t, double x, double v) {
		return -pow(x, 1.5) - 2.0 / t * v;
	}

	void rk4(PS::F64 y[], PS::F64 dydx[], PS::S32 n, PS::F64 x, PS::F64 h,
			PS::F64 yout[], PS::F64 N)
			/* Given values for the variables y[1..n] and their derivatives dydx[1..n]
			 known at x, use the fourth-order Runge-Kutta method to advance the solution
			 over an PS::S32 erval h and return the incremented variables as yout[1..n], which
			 need not be a distinct array from y. The user supplies the routine
			 derivs(x,y,dydx), which returns derivatives dydx at x. */
			{
		PS::S32 i;
		PS::F64 xh, hh, h6, dym[NSIZE], dyt[NSIZE], yt[NSIZE];

		hh = h * 0.5;
		h6 = h / 6.0;
		xh = x + hh;

		for (i = 0; i < n; i++) {
			yt[i] = y[i] + hh * dydx[i];			// First step.
		}
		dyt[0] = yt[1] / (xh) / (xh);
		dyt[1] = -(xh) * (xh) * pow(yt[0], N);
		// Second step.
		for (i = 0; i < n; i++) {
			yt[i] = y[i] + hh * dyt[i];
		}
		dym[0] = yt[1] / (xh) / (xh);
		dym[1] = -(xh) * (xh) * pow(yt[0], N);
		for (i = 0; i < n; i++) {
			yt[i] = y[i] + h * dym[i];
			dym[i] += dyt[i];
		}
		dyt[0] = yt[1] / (x + h) / (x + h);
		dyt[1] = -(x + h) * (x + h) * pow(yt[0], N);

		for (i = 0; i < n; i++) {
			yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
		}
	}
	void calc_Vij2_and_ss(const EPI::Hydro ep_i, const EPJ::Hydro ep_j,
			PS::F64&Vij2_hi, PS::F64&Vij2_hj, PS::F64& ssP, const PS::F64 delta,
			const PS::F64vec eij) {
		//		const PS::F64 dVi = ep_i.gradV * eij;
		//		const PS::F64 dVj = ep_j.gradV * eij;
		const PS::F64 Vi = 1.0 / ep_i.dens;
		const PS::F64 Vj = 1.0 / ep_j.dens;
		const PS::F64 hi2 = ep_i.smth * ep_i.smth;
		const PS::F64 hj2 = ep_j.smth * ep_j.smth;
		//		const PS::F64 rhoij = 0.5 * (ep_i.dens + ep_j.dens);
		//		const PS::F64 Vij2crit = 10.0 * (.0 / (rhoij * rhoij));
		const PS::F64 Pij = 0.5 * (ep_i.pres + ep_j.pres);
		PS::F64 Aij, Bij, Cij, Dij;

		Cij = (Vi - Vj) / delta;
		Dij = 0.5 * (Vi + Vj);
		Vij2_hi = 0.25 * hi2 * Cij * Cij + Dij * Dij;
		Vij2_hj = 0.25 * hj2 * Cij * Cij + Dij * Dij;
		ssP = .5
				* (hi2 * Cij * Dij / (2 * Vij2_hi)
						+ hj2 * Cij * Dij / (2 * Vij2_hj));
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
			const PS::F64 dt, PS::F64 &pstar, PS::F64 &vstar,
			) {
		PS::F64 si = .5 * delta;
		PS::F64 sj = .5 * delta;
		PS::F64 drdsR = ep_R.grad_dens * eij;
		PS::F64 drdsL = ep_L.grad_dens * eij;
		PS::F64 dpdsR = ep_R.gradP * eij;
		PS::F64 dpdsL = ep_L.gradP * eij;
		PS::F64 vR = ep_R.vel * eij;
		PS::F64 vL = ep_L.vel * eij;

		PS::F64 dvdsR = ep_R.gradvel_x * eij + ep_R.gradvel_y * eij
				+ ep_R.gradvel_z * eij;
		PS::F64 dvdsL = ep_L.gradvel_x * eij + ep_L.gradvel_y * eij
				+ ep_L.gradvel_z * eij;
		PS::F64 ddensi = drdsR * (ss + ep_R.snds * 0.5 * dt - si);
		PS::F64 dpresi = dpdsR * (ss + ep_R.snds * 0.5 * dt - si);
		PS::F64 dveli = dvdsR * (ss + ep_R.snds * 0.5 * dt - si);
		PS::F64 ddensj = drdsL * (ss - ep_L.snds * 0.5 * dt + sj);
		PS::F64 dpresj = dpdsL * (ss - ep_L.snds * 0.5 * dt + sj);
		PS::F64 dvelj = dvdsL * (ss - ep_L.snds * 0.5 * dt + sj);

		PS::F64 _rhoR = ep_R.dens + ddensi;
		PS::F64 _pR = ep_R.pres + dpresi;
		PS::F64 _vR = vR + dveli;
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

		PS::F64 p = .5 * (pRs + pLs);
		pstar = p;
		PS::F64 ppre = p;
		vstar = .5 * (_vR + _vL);
		for (PS::U32 loop = 0; loop < 10; loop++) {
			ppre = p;
			if ((fabs(p / pRs) + 1.0e-9 < critif)) {
				W1 = sqrt(pRs * _rhoR)
						* ((PARAM::GAMMA - 1.0) / (2.0 * sqrt(PARAM::GAMMA)))
						* (1.0 - p / pRs) / (1.0 - pow(p / pRs, 1.0 / alpha));
			} else {
				W1 = sqrt(pRs * _rhoR)
						* sqrt(
								0.5 * (PARAM::GAMMA + 1.0) * p / (pRs)
										+ 0.5 * (PARAM::GAMMA - 1.0));
			}
			if ((fabs(p / pLs) + 1.0e-9 < critif)) {
				W2 = sqrt(pLs * _rhoL)
						* ((PARAM::GAMMA - 1.0) / (2.0 * sqrt(PARAM::GAMMA)))
						* (1.0 - p / pLs) / (1.0 - pow(p / pLs, 1.0 / alpha));
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
		vstar = .5 * (vR + vL + dveli + dvelj);
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////

		if (3.0 * fabs(vL - vR)
				> ((ep_R.snds < ep_L.snds) ? ep_R.snds : ep_L.snds)) {

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
		}

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
				if (ith.id != jth.id) {
					PS::F64 VIJ_I, VIJ_J, PSTAR, VSTAR;
					PS::F64 sstar = 0.0;
					PS::F64vec interpo = 0.0;
					const PS::F64vec dv = ith.vel - jth.vel;

					const PS::F64vec eij = dr / delta;
					const PS::F64vec gradW_hi = kernel.gradW(dr,
							sqrt(2.0) * ith.smth);
					const PS::F64vec gradW_hj = kernel.gradW(dr,
							sqrt(2.0) * jth.smth);

					calc_Vij2_and_ss(ith, jth, VIJ_I, VIJ_J, sstar, delta, eij);

					calc_riemann_solver(ith, jth, sstar, delta, eij, ith.dt,
							PSTAR, VSTAR);
					interpo = (PSTAR) * (gradW_hi * VIJ_I + gradW_hj * VIJ_J);
					hydro[i].acc -= jth.mass * interpo;
				}

			}
			double dt = ep_i[i].dt;
			PS::F64vec vel_half = ith.vel + 0.5 * hydro[i].acc * dt;
			for (PS::S32 j = 0; j < Njp; ++j) {
				const EPJ::Hydro& jth = ep_j[j];
				const PS::F64vec dr = ith.pos - jth.pos;
				const PS::F64 delta = sqrt(dr * dr);

				if (ith.id != jth.id) {
					PS::F64 VIJ_I, VIJ_J, PSTAR, VSTAR;
					PS::F64 sstar = 0.0;
					PS::F64vec interpo = 0.0;
					PS::F64vec eij, evij;
					const PS::F64vec dv = ith.vel - jth.vel;
					const PS::F64 deltav = sqrt(dv * dv);
					const PS::F64vec gradW_hi = kernel.gradW(dr,
							sqrt(2.0) * ith.smth);
					const PS::F64vec gradW_hj = kernel.gradW(dr,
							sqrt(2.0) * jth.smth);

					evij = deltav == 0.0 ? 0.0 : dv / deltav;
					eij = dr / delta;
					PS::F64vec vij = 0.0;
					calc_Vij2_and_ss(ith, jth, VIJ_I, VIJ_J, sstar, delta, eij);

					calc_riemann_solver(ith, jth, sstar, delta, eij, ith.dt,
							PSTAR, VSTAR);
					interpo = (PSTAR) * (gradW_hi * VIJ_I + gradW_hj * VIJ_J);

					vij = VSTAR * eij;
					hydro[i].eng_dot -= jth.mass * (vij - vel_half) * interpo;
					hydro[i].dt =
							0.1 * ep_i[i].pres
									/ fabs(
											.5 * ep_i[i].grav * eij
													* sqrt(
															PARAM::GAMMA
																	* ep_i[i].pres
																	* ep_i[i].dens)
													+ .5 * eij * ep_i[i].extF
															* sqrt(
																	PARAM::GAMMA
																			* ep_i[i].pres
																			* ep_i[i].dens)
															* dt);
					hydro[j].dt =
							0.1 * ep_j[j].pres
									/ fabs(
											.5 * ep_j[j].grav * eij
													* sqrt(
															PARAM::GAMMA
																	* ep_j[j].pres
																	* ep_j[j].dens)
													+ .5 * ep_j[j].extF * eij
															* sqrt(
																	PARAM::GAMMA
																			* ep_j[j].pres
																			* ep_j[j].dens));

				}

			}

			PS::F64vec bh_loc(0, 0, 0);

			const PS::F64 _r2 = (ep_i[i].pos - bh_loc) * (ep_i[i].pos - bh_loc);
			const PS::F64 _r = sqrt(_r2 + ith.smth * ith.smth);

			const PS::F64 _r3 = _r * _r * _r;

			PS::F64 xtilde = _r / ith.smth;
			PS::F64 xtilde2 = xtilde * xtilde;
			PS::F64 formfac = 0.0;

			formfac = erf(xtilde) - exp(-xtilde2) * xtilde * 2 / sqrt(M_PI);

			formfac = 1.0;
			PS::F64vec BHGrav = (ep_i[i].pos - bh_loc) * formfac * PARAM::G
					* PARAM::sM_BH / (_r3);
			hydro[i].extF = -BHGrav;

		}

	}

}
;
//template<class TParticleJ> class CalcGravityForce {
//public:
//	void operator ()(const EPI::Grav* const ep_i, const PS::S32 Nip,
//			const TParticleJ* const ep_j, const PS::S32 Njp,
//			RESULT::Grav* const grav) {
//		for (PS::S32 i = 0; i < Nip; ++i) {
//			const EPI::Grav& ith = ep_i[i];
//			for (PS::S32 j = 0; j < Njp; ++j) {
//				const TParticleJ& jth = ep_j[j];
//				const PS::F64vec dr = ith.pos - jth.pos;
//				const PS::F64 dr2 = dr * dr;
//				const PS::F64 dr_inv = PARAM::G / sqrt(dr2 + ith.getEps2());
//				const PS::F64 m_dr3_inv = jth.mass * math::pow3(dr_inv);
//				grav[i].acc += -m_dr3_inv * dr;
//				grav[i].pot += -jth.mass * dr_inv;
//			}
//		}
//	}
//};

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
					const PS::F64 dr_inv = 1.0 / sqrt(dr2 + 1e-18);

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
//
