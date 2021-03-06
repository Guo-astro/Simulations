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
					std::cout << dens[i].dens << std::endl;
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
	void calc_riemann_solver(const EPI::Hydro ep_i, const EPJ::Hydro ep_j,
			const PS::F64 ss, const PS::F64 delta, const PS::F64vec eij,
			const PS::F64 dt, PS::F64 &pstar, PS::F64 &vstar) {
		PS::F64 rho1, rho2, v1, v2, p1, p2;
		PS::F64 vi = ep_i.vel * eij;
		PS::F64 vj = ep_j.vel * eij;
		PS::F64 drdsi, drdsj, dpdsi, dpdsj, dvdsi, dvdsj;

		PS::F64 si = .5 * delta;
		PS::F64 sj = .5 * delta;
		drdsi = ep_i.grad_dens * eij;
		drdsj = ep_j.grad_dens * eij;
		dpdsi = ep_i.gradP * eij;
		dpdsj = ep_j.gradP * eij;
		if (dvdsi * dvdsj < 0.0 || drdsi * drdsj < 0.0 || dpdsi * dpdsj < 0.0) {
			dvdsi = 0.0;
			dvdsj = 0.0;
		}
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

		rho1 = ep_i.dens + ddensi;
		p1 = ep_i.pres + dpresi;
		v1 = vi + dveli;
		rho2 = ep_j.dens + ddensj;
		p2 = ep_j.pres + dpresj;
		v2 = vj + dvelj;
		if (rho1 < 0.0 || p1 < 0.0 || rho2 < 0.0 || p2 < 0.0) {
			rho1 = 0.0;
			p1 = 0.0;
			v1 = 0.0;
			rho2 = 0.0;
			p2 = 0.0;
			v2 = 0.0;
		}
		//
		PS::F64 ppre, p, v;
		PS::F64 W1, W2;
		const PS::F64 alpha = (2.0 * PARAM::GAMMA) / (PARAM::GAMMA - 1.0);
		p = .5 * (p1 + p2);
		pstar = p;
		vstar = .5 * (v1 + v2);
		double critif = 1.0 - 1.0e-6;

		double p1s = p1
				- .5 * eij * ep_i.grav * sqrt(PARAM::GAMMA * p1 * rho1) * dt;
		double p2s = p2
				+ .5 * eij * ep_j.grav * sqrt(PARAM::GAMMA * p2 * rho2) * dt;

		if (p1 <= 1.0e-6 && p2 <= 1.0e-6) {
			pstar = 0.0;
			vstar = 0.0;

		} else {
			for (PS::U32 loop = 0; loop < 40; loop++) {
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
				if (p < 0.0) {
					p = 0.5 * ppre;
					break;
				}
				if (fabs(1.0 - p / ppre) < 1e-3)
					break;
			}
			vstar = ((W1 * v1 + W2 * v2) + p2s - p1s) / (W1 + W2);
			pstar = p;
			if ((pstar != pstar) || vstar != vstar || pstar <= 0.0) {
				vstar = .5 * (v1 + v2);
				pstar = .5 * (p1 + p2);
				//						std::cout <<p1<<"=="<<p2<<"=="<<pstar<<std::endl;
			}
		}
	}
	double getPhi(double x) {

		double a = 1.000744397;

		double b = 1.50738116e-2;

		double c = -2.334693e-1;

		double d = 7.8109735e-2;

		double e = -1.021152772e-2;

		double f = 4.833905487e-4;

		double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3)
				+ e * pow(x, 4) + f * pow(x, 5);

		return phi;

	}

	double getPhi_dash(double x) {

		double a = 8.841062613e-3;

		double b = -3.86087802e-1;

		double c = 6.92279238e-2;

		double d = 6.159126e-3;

		double e = -2.44979111e-3;

		double f = 1.6204056e-4;

		double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3)
				+ e * pow(x, 4) + f * pow(x, 5);

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
				const PS::F64 _r = sqrt(ep_i[i].pos * ep_i[i].pos);

				double current_phi = getPhi(_r);
				double current_phi_dash = getPhi_dash(_r);
				PS::F64vec acc_factor = ep_i[i].pos / _r;

				hydro[i].acc += (PARAM::GAMMA - 1) * ith.eng * current_phi_dash
						* acc_factor;
				hydro[i].dt = 0.1 * ith.smth / ith.snds;
				hydro[i].dt = fmin(hydro[i].dt,
						0.02 * sqrt(ith.smth / sqrt(ith.grav * ith.grav)));

			}
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

