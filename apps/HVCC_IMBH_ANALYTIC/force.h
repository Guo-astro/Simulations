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
			hmin = fmin(hmin, dens[i].smth);

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
struct interP {
	PS::F64 VIJ_I;
	PS::F64 VIJ_J;
	PS::F64 sstar;
	PS::F64 PSTAR;
	PS::F64 VSTAR;

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
	void calc_riemann_solver(const EPI::Hydro ep_i, const EPJ::Hydro ep_j,
			const PS::F64 ss, const PS::F64 delta, const PS::F64vec eij,
			const PS::F64 dt, PS::F64 &pstar, PS::F64 &vstar) {
		PS::F64 rho1, rho2, v1, v2, p1, p2;
		PS::F64 vi = ep_i.vel * eij;
		PS::F64 vj = ep_j.vel * eij;
		PS::F64 drdsi, drdsj, dpdsi, dpdsj, dvdsi, dvdsj;

		PS::F64 si = .5 * delta;
		PS::F64 sj = .5 * delta;
		drdsi = ep_i.graddens * eij;
		drdsj = ep_j.graddens * eij;
		dpdsi = ep_i.gradpres * eij;
		dpdsj = ep_j.gradpres * eij;
		dvdsi = ep_i.gradvel_x * eij * eij.x + ep_i.gradvel_y * eij * eij.y
				+ ep_i.gradvel_z * eij * eij.z;
		dvdsj = ep_j.gradvel_x * eij * eij.x + ep_j.gradvel_y * eij * eij.y
				+ ep_j.gradvel_z * eij * eij.z;
		if (dvdsi * dvdsj < 0.0 || drdsi * drdsj < 0.0 || dpdsi * dpdsj < 0.0) {
			dvdsi = 0.0;
			dvdsj = 0.0;
		}
//		if (dvdsi * dvdsj < 0.0) {
//			dvdsi = dvdsj = 0.0;
//			drdsi = drdsj = 0.0;
//			dpdsi = dpdsj = 0.0;
//		}
		if (3.0 * fabs(vj - vi)
				> ((ep_i.snds < ep_j.snds) ? ep_i.snds : ep_j.snds)) {
			//		//monotonisity constraint

			drdsi = drdsj = 0.0;
			dpdsi = dpdsj = 0.0;
			dvdsi = dvdsj = 0.0;
		}
		double ddensi = drdsi * (ss + ep_i.snds * 0.5 * dt - si);
		double dpresi = dpdsi * (ss + ep_i.snds * 0.5 * dt - si);
		double dveli = dvdsi * (ss + ep_i.snds * 0.5 * dt - si);
		double ddensj = drdsj * (ss - ep_j.snds * 0.5 * dt + sj);
		double dpresj = dpdsj * (ss - ep_j.snds * 0.5 * dt + sj);
		double dvelj = dvdsj * (ss - ep_j.snds * 0.5 * dt + sj);

		// ddensi =  0.0;
		// dpresi =  0.0;
		// ddensj =  0.0;
		// dpresj =  0.0;
		// dveli =   0.0;
		// dvelj =   0.0;

		rho1 = ep_i.dens + ddensi;
		p1 = ep_i.pres + dpresi;
		v1 = vi + dveli;
		rho2 = ep_j.dens + ddensj;
		p2 = ep_j.pres + dpresj;
		v2 = vj + dvelj;
////				std::cout<<dpresi<<std::endl;
		if (rho1 < 0.0 || p1 < 0.0 || rho2 < 0.0 || p2 < 0.0) {
			rho1 = 0.0;
			p1 = 0.0;
			v1 = 0.0;
			rho2 = 0.0;
			p2 = 0.0;
			v2 = 0.0;
		}
//		rho1 = ep_i.dens+ ddensj;
//		p1 = ep_i.pres+ dpresj;
//
//		v1 = vi + dvelj;
//		rho2 = ep_j.dens + ddensi;
//		p2 = ep_j.pres+ dpresi;
//		v2 = vj+dveli;
//		if (rho1 < 0.0 || p1 < 0.0 || rho2 < 0.0 || p2 < 0.0) {
//			rho1 = ep_i.dens;
//			p1 = ep_i.pres ;
//			v1 = vi;
//			rho2 = ep_j.dens;
//			p2 = ep_j.pres ;
//			v2 = vj;
//		}
//		PS::F64vec dr = ep_i.pos - ep_j.pos;
//						rho1 = ep_i.dens;
//						p1 = ep_i.pres-ep_i.grav*eij * sqrt(PARAM::GAMMA * p1 * rho1) * dt;
//
//						v1 = vi;
//						rho2 = ep_j.dens;
//						p2 = ep_j.pres  + ep_j.grav*eij * sqrt(PARAM::GAMMA * p2 * rho2) * dt;
//						v2 = vj;
//
		PS::F64 ppre, p, v;
		PS::F64 W1, W2;
		const PS::F64 alpha = (2.0 * PARAM::GAMMA) / (PARAM::GAMMA - 1.0);
		p = .5 * (p1 + p2);
		pstar = p;
		vstar = .5 * (v1 + v2);
		double critif = 1.0 - 1.0e-6;
		double p1s = p1
				- .5 * eij * ep_i.grav * sqrt(PARAM::GAMMA * p1 * rho1) * dt
				- .5 * eij * ep_i.extF * sqrt(PARAM::GAMMA * p1 * rho1) * dt;
		double p2s = p2
				+ .5 * eij * ep_j.grav * sqrt(PARAM::GAMMA * p2 * rho2) * dt
				+ .5 * eij * ep_j.extF * sqrt(PARAM::GAMMA * p2 * rho2) * dt;
//		double p1s = p1;
//		double p2s = p2;
//		double p1s = p1 - .5 * eij * ep_i.grav * sqrt(PARAM::GAMMA * p1 * rho1) * dt;
//		double p2s = p2 + .5 * eij * ep_j.grav * sqrt(PARAM::GAMMA * p2 * rho2) * dt;
//p1 -=.5 * eij * ep_i.grav * sqrt(PARAM::GAMMA * p1 * rho1) * dt;
//p2+=.5 * eij * ep_j.grav * sqrt(PARAM::GAMMA * p2 * rho2) * dt;
//		double p1s = p1
//				- .5 * ep_i.grav * eij * sqrt(PARAM::GAMMA * p1 * rho1) * dt;
//		double p2s = p2
//				+ .5 * ep_j.grav * eij * sqrt(PARAM::GAMMA * p2 * rho2) * dt;
//		double p1s = p1 - ep_i.grav * eij * sqrt(PARAM::GAMMA * p1 * rho1) * dt;
//		double p2s = p2 + ep_j.grav * eij * sqrt(PARAM::GAMMA * p2 * rho2) * dt;
//						double p1s = p1 -ep_i.grav*eij *6.0 * sqrt(PARAM::GAMMA * p1 * rho1) * dt;
//							double p2s = p2 + ep_j.grav*eij * 6.0 * sqrt(PARAM::GAMMA * p2 * rho2) * dt;
		p1 += 1e-8;
		p2 += 1e-8;
		if (p1s <= 0.0 || p2s <= 0.0) {
			pstar = .5 * (p1 + p2);
			vstar = .5 * (v1 + v2);
		} else {
			for (PS::U32 loop = 0; loop < 10; loop++) {
				ppre = p;
				if ((p / p1 < critif)) {
					W1 =
							sqrt(p1 * rho1)
									* ((PARAM::GAMMA - 1.0)
											/ (2.0 * sqrt(PARAM::GAMMA)))
									* (1.0 - p / p1)
									/ (1.0 - pow(p / p1, 1.0 / alpha));
				} else {
					W1 = sqrt(p1 * rho1)
							* sqrt(
									0.5 * (PARAM::GAMMA + 1.0) * p / (p1s)
											+ 0.5 * (PARAM::GAMMA - 1.0));
				}
				if ((p / p2 < critif)) {
					W2 =
							sqrt(p2 * rho2)
									* ((PARAM::GAMMA - 1.0)
											/ (2.0 * sqrt(PARAM::GAMMA)))
									* (1.0 - p / p2)
									/ (1.0 - pow(p / p2, 1.0 / alpha));
				} else {
					W2 = sqrt(p2 * rho2)
							* sqrt(
									0.5 * (PARAM::GAMMA + 1.0) * p / (p2s)
											+ 0.5 * (PARAM::GAMMA - 1.0));
				}
				p = ((p2s / W2 + p1s / W1) + v2 - v1) / (1.0 / W2 + 1.0 / W1);
				if (p < 0.0)
					p = 0.5 * ppre;
				if (fabs(1.0 - p / ppre) < 1e-3)
					break;
			}
			vstar = ((W1 * v1 + W2 * v2) + p2s - p1s) / (W1 + W2);
			pstar = p;
			if ((pstar != pstar) || vstar != vstar || pstar <= 0.0) {
				vstar = .5 * (v1 + v2);
				pstar = .5 * (p1 + p2);
				std::cout << p1 << "/// " << p2 << "/// " << p1s << "/// "
						<< "/// " << p2s << "/// " << W1 << "/// " << W2
						<< std::endl;
			}
		}
	}
	void operator ()(const EPI::Hydro* const ep_i, const PS::S32 Nip,
			const EPJ::Hydro* const ep_j, const PS::S32 Njp,
			RESULT::Hydro* const hydro) {

		for (PS::S32 i = 0; i < Nip; ++i) {
			hydro[i].clear();
			PS::F64 v_sig_max = 0.0;
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
			//	hydro[i].acc += hydro[i].extF + ep_i[i].grav;
			double dt = ep_i[i].dt;
			hydro[i].vel_half_old = ith.vel
					+ 0.5 * (hydro[i].acc + ep_i[i].grav + ep_i[i].extF) * dt;
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
					hydro[i].eng_dot -= jth.mass * (vij - hydro[i].vel_half_old)
							* interpo;

				}

			}
//hydro[i].eng_dot -=  hydro[i].vel_half_old* (hydro[i].extF + ep_i[i].grav);

			PS::F64vec bh_loc(0, 0, 0);

			const PS::F64 _r2 = (ep_i[i].pos - bh_loc) * (ep_i[i].pos - bh_loc);
			const PS::F64 _r = sqrt(_r2);

			const PS::F64 _r3 = _r * _r * _r;

			PS::F64 xtilde = _r/ith.smth;
			PS::F64 xtilde2 = xtilde * xtilde;
			PS::F64 formfac = 0.0;

			formfac = erf(xtilde) - exp(-xtilde2) * xtilde * 2 / sqrt(M_PI);

//formfac=1.0;
			PS::F64vec BHGrav = (ep_i[i].pos - bh_loc) *formfac* PARAM::G * PARAM::sM_BH / (_r3);
			hydro[i].extF = -BHGrav;


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
