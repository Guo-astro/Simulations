#pragma once

/*
 struct boundary{
 PS::F64 x, y, z;
 };
 */

class FileHeader {
public:
	int Nbody;
	double time;
	double phy_time;

	int readAscii(FILE* fp) {
		fscanf(fp, "%lf\n", &time);
		fscanf(fp, "%d\n", &Nbody);
		return Nbody;
	}
	void writeAscii(FILE* fp) const {
		fprintf(fp, "%e\n", time);
		fprintf(fp, "%d\n", Nbody);
	}
#ifdef OUTPUT_VTK
	void writeVTKAscii(FILE* fp) const {

		fprintf(fp, "# vtk DataFile Version 3.0\n");
		fprintf(fp, "Example 3D regular grid VTK file.\n");
		fprintf(fp, "ASCII\n");
		fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
		fprintf(fp, "FIELD FieldData 1\n");
		fprintf(fp, "TIME 1 1 double\n");
		fprintf(fp, "%e\n", phy_time);
		fprintf(fp, "POINTS %d double\n", Nbody);
//		fprintf(fp, "%e\n", time);
	}
	void writeAscii_scalars_pre_header(FILE* fp) const {

		fprintf(fp, "POINT_DATA %d\n", Nbody);
	}
	void writeAscii_scalars_header(const char* name,FILE* fp) const {

		fprintf(fp, "SCALARS %s double 1\n",name);
		fprintf(fp, "LOOKUP_TABLE default\n");
	}
	void writeAscii_vectors_header(const char* name,FILE* fp) const {

		fprintf(fp, "VECTORS %s double\n",name);
//			fprintf(fp, "LOOKUP_TABLE default\n");
	}
	void writeAsciiTxt(FILE* fp) const {

		fprintf(fp, "%d\n", Nbody);
		fprintf(fp, "%e\n", time);
	}
#endif
};

namespace RESULT {
//Density summation
class Dens {
public:
	PS::F64 dens;
	PS::F64 smth;
	void clear() {
		dens = 0;
	}
};
class Drvt {
public:
	PS::F64vec grad_dens;
	PS::F64vec gradP;
	PS::F64vec gradPT;
	PS::F64vec gradVpara;
	PS::F64vec gradBperp2;
	PS::F64vec gradvperp_x, gradvperp_y, gradvperp_z;
	PS::F64vec gradBperp_x, gradBperp_y, gradBperp_z;
	PS::F64 div_v;
	PS::F64vec rot_v;
	PS::F64vec gradV;
	PS::F64vec graddens;
	PS::F64vec gradpres;
	PS::F64vec gradvel_x;
	PS::F64vec gradvel_y;
	PS::F64vec gradvel_z;
	PS::F64vec gradvel;
	void clear() {
		grad_dens = 0.0;
		gradP = 0.0;
		gradPT = 0.0;
		gradVpara = 0.0;
		gradBperp2 = 0.0;
		gradvperp_x = gradvperp_y = gradvperp_z = 0.0;
		gradBperp_x = gradBperp_y = gradBperp_z = 0.0;
		div_v = 0.0;
		rot_v = 0.0;
		gradV = 0.0;
		graddens = 0.0;
		gradpres = 0.0;
		gradvel_x = 0.;
		gradvel_y = 0.;
		gradvel_z = 0.;
		gradvel = 0.0;
	}
};
//Hydro force

class Hydro {
public:
	PS::F64vec vel_half_old;

	PS::F64vec acc;
	PS::F64 eng_dot;
	PS::F64vec BoverDens_dot;
	PS::F64 dt;
	PS::F64vec extF;
	void clear() {
		extF = 0.;
		acc = 0.;
		eng_dot = 0.;
		BoverDens_dot = 0.;
		dt = 1.0e+30;
	}
};

//Self gravity
class Grav {
public:
	PS::F64vec acc;
	PS::F64 pot;
	PS::F64 dt;

	void clear() {
		acc = 0.0;
		pot = 0.0;
		dt = 1.0e+30;
	}
};
}

class RealPtcl {
public:
	PS::F64 mass;
	PS::F64vec pos, vel, acc;
	PS::F64vec extF;
	PS::F64 dens; //DENSity
	PS::F64 eng; //ENerGy
	PS::F64 pres; //PRESsure
	PS::F64 smth; //SMooTHing length
	PS::F64 snds; //SouND Speed
	PS::F64 Lambda;
	PS::F64 Gamma;
	PS::F64 cooling_timescale, sg_timescale, bh_timescale, hydro_timescale,
			eng_timescale;
	PS::F64 div_v;
	PS::F64vec rot_v;
	PS::F64 Bal; //Balsala switch
	PS::F64vec grad_dens;
	PS::F64vec gradP;
	PS::F64vec gradPT;
	PS::F64vec gradvperp_x, gradvperp_y, gradvperp_z;
	PS::F64vec gradBperp_x, gradBperp_y, gradBperp_z;
	PS::F64vec gradBperp2;
	PS::F64vec gradVpara;
	PS::F64vec xstar;
	PS::F64vec gradV; //gradient of specific volume
	PS::F64vec graddens; //gradient of density
	PS::F64vec gradpres; //gradient of pressure
	PS::F64vec gradvel_x;
	PS::F64vec gradvel_y;
	PS::F64vec gradvel_z;
	PS::F64vec gradvel;
	PS::F64vec MagneticB;
	PS::F64vec BoverDens_dot;
	PS::F64vec BoverDens;
	PS::F64vec BoverDens_half;
	PS::F64 eng_dot;
	PS::F64vec vel_half;
	PS::F64 eng_half;
	PS::F64vec magneticB_half;
	PS::F64 dt;
	PS::S64 id;
	PS::S64 minus_eng;
	PS::F64vec grav;
	PS::F64 pot;

	const char* scalars[PARAM::COMP + 11] = { "mass", "pres", "dens", "vx",
			"vy", "vz", "T", "eng", "ab_e", "ab_HI", "ab_HeI", "ab_OI", "ab_CI",
			"ab_H2", "ab_CO", "ab_HII", "ab_CII", "ab_FeII", "ab_SiII",
			"cooling", "heating", "cooling_timescale" };
	const char* vectors[3] = { "vel", "acc", "pos" };

	PS::F64 abundances[PARAM::COMP];

//Copy functions
	void copyFromForce(const RESULT::Dens& dens) {
		this->dens = dens.dens;
		this->smth = dens.smth;
	}
	void copyFromForce(const RESULT::Drvt& drvt) {

		this->gradV = drvt.gradV;
		this->graddens = drvt.graddens;
		this->gradpres = drvt.gradpres;
		this->gradvel_x = drvt.gradvel_x;
		this->gradvel_y = drvt.gradvel_y;
		this->gradvel_z = drvt.gradvel_z;
		this->gradvel = drvt.gradvel;

		this->gradPT = drvt.gradPT;
		this->gradBperp2 = drvt.gradBperp2;
		this->gradVpara = drvt.gradVpara;
		this->gradP = drvt.gradP;
		this->grad_dens = drvt.grad_dens;
		this->gradBperp_x = drvt.gradBperp_x;
		this->gradBperp_y = drvt.gradBperp_y;
		this->gradBperp_z = drvt.gradBperp_z;
		this->gradvperp_x = drvt.gradvperp_x;
		this->gradvperp_y = drvt.gradvperp_y;
		this->gradvperp_z = drvt.gradvperp_z;
		this->div_v = drvt.div_v;
		this->rot_v = drvt.rot_v;
		this->Bal = fabs(drvt.div_v)
				/ (fabs(drvt.div_v) + sqrt(drvt.rot_v * drvt.rot_v)
						+ 1.0e-4 * this->snds / this->smth); //Balsala switch
	}
	void copyFromForce(const RESULT::Hydro& force) {
		this->acc = force.acc;
		this->eng_dot = force.eng_dot;
		this->extF = force.extF;
		this->dt = force.dt;
		this->BoverDens_dot = force.BoverDens_dot;
	}
	void copyFromForce(const RESULT::Grav& force) {
		this->grav = force.acc;
		this->pot = force.pot;
		//not copy dt
	}
//Give necessary values to FDPS
	PS::F64 getCharge() const {
		return this->mass;
	}
	PS::F64vec getPos() const {
		return this->pos;
	}
	PS::F64 getRSearch() const {
		return kernel_t::supportRadius() * this->smth;
	}
	void setPos(const PS::F64vec& pos) {
		this->pos = pos;
	}
	void writeAscii(FILE* fp) const {
		const PS::F64 radial_mag = sqrt(this->pos * this->pos);
		const PS::F64vec radial = this->pos / radial_mag;
		const PS::F64 radialForce_mag = this->acc * (radial);
		const PS::F64 radialForce_grav_mag = this->grav * radial;
		const PS::F64 radialForce_extF_mag = this->extF * radial;

		double mu = this->abundances[1] + this->abundances[7]
				+ 2.0 * this->abundances[5] + 4.0 * this->abundances[2];
		mu /= (this->abundances[1] + this->abundances[7] + this->abundances[5]
				+ this->abundances[2] + this->abundances[0]);
		PS::F64 NUMDENS_CGS = this->dens * PARAM::SMassDens
				/ (mu * PARAM::PROTONMASS);
		PS::F64 tem = mu * PARAM::PROTONMASS_CGS * this->pres
				* PARAM::SEng_per_Mass
				/ (this->dens * PARAM::KBOLTZ_cgs
						* (1.1 + PARAM::i_abundance_e - PARAM::i_abundance_H2));
		//14
		//15,16,17 ,18->grav_force
		//23 sg,24 bh,25 hydro 26 eng
		//38 minus_eng
		fprintf(fp,
				"%lld\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e"
						"\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n",
				this->id, this->mass, this->smth, this->pos.x, this->pos.y,
				this->pos.z, this->vel.x, this->vel.y, this->vel.z, this->dens,
				this->eng, this->pres, this->MagneticB.x, this->MagneticB.y,
				this->MagneticB.z, radial_mag, radialForce_mag,
				radialForce_grav_mag, radialForce_extF_mag, tem, NUMDENS_CGS,
				this->Lambda, this->cooling_timescale, this->sg_timescale,
				this->bh_timescale, this->hydro_timescale, this->eng_timescale,
				this->abundances[0], this->abundances[1], this->abundances[2],
				this->abundances[3], this->abundances[4], this->abundances[5],
				this->abundances[6], this->abundances[7], this->abundances[8],
				this->abundances[9], mu, this->minus_eng);
	}
	void readAscii(FILE* fp) {
#ifdef RESTART
		PS::F64 radial_mag = sqrt(this->pos * this->pos);
		PS::F64vec radial = this->pos / radial_mag;
		PS::F64 radialForce_mag = this->acc * (radial);
		PS::F64 radialForce_grav_mag = this->grav * radial;
		PS::F64 mu, tem, NUMDENS_CGS;
		PS::S64 minus_eng;
		fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf"
				"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &this->id, &this->mass, &this->smth, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y,
				&this->vel.z, &this->dens, &this->eng, &this->pres, &this->MagneticB.x, &this->MagneticB.y, &this->MagneticB.z, &radial_mag, &radialForce_mag, &radialForce_grav_mag, &tem, &NUMDENS_CGS, &this->Lambda, &this->cooling_timescale ,&this->sg_timescale,&this->bh_timescale,&this->hydro_timescale,&this->eng_timescale,&this->abundances[0], &this->abundances[1], &this->abundances[2],
				&this->abundances[3], &this->abundances[4], &this->abundances[5], &this->abundances[6], &this->abundances[7], &this->abundances[8], &this->abundances[9], &mu,&minus_eng);
		if (this->eng < 0.0 || this->pres < 0.0) {
			std::cout << pres << "in class.h readAscii" << eng<<" "<<dens << std::endl;

		}

#endif
//		std::cout << this->mass << std::endl;
#ifdef INIT_FROM_RELAXED_PROFILE
		PS::F64 radial_mag = sqrt(this->pos * this->pos);
		PS::F64vec radial = this->pos / radial_mag;
		PS::F64 radialForce_mag = this->acc * (radial);
		PS::F64 radialForce_grav_mag = this->grav * radial;
		fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &this->id, &this->mass,&this->smth, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z, &this->dens, &this->eng, &this->pres, &this->MagneticB.x, &this->MagneticB.y,
				&this->MagneticB.z, &radial_mag, &radialForce_mag, &radialForce_grav_mag);
		double mu = 2.34;
		double mu_new = 2.34;
		PS::F64 NUMDENS_CGS = 0;
		for (int i = 0; i < 100; i++) {
			mu = mu_new;
			NUMDENS_CGS = dens * PARAM::SMassDens / (mu * PARAM::PROTONMASS);

			this->abundances[0] = fitAB_E(NUMDENS_CGS);
			this->abundances[1] = fitAB_HI(NUMDENS_CGS);
			this->abundances[2] = PARAM::i_abundance_HeI;
			this->abundances[3] = PARAM::i_abundance_OI;
			this->abundances[4] = PARAM::i_abundance_CI;
			this->abundances[5] = fitAB_H2(NUMDENS_CGS);
			this->abundances[6] = fitAB_CO(NUMDENS_CGS);
			this->abundances[7] = fmax(1. - 2. * this->abundances[5] - this->abundances[1], 0);
			this->abundances[8] = fmax(this->abundances[0] - this->abundances[7] - PARAM::i_abundance_SiII, 0);
			this->abundances[9] = PARAM::i_abundance_FeII;
			this->abundances[10] = PARAM::i_abundance_SiII;
			mu_new = this->abundances[1] + this->abundances[7] + 2.0 * this->abundances[5] + 4.0 * this->abundances[2];
			mu_new /= (this->abundances[1] + this->abundances[7] + this->abundances[5] + this->abundances[2] + this->abundances[0]);
			if (fabs(mu - mu_new) < 1e-3) {
				break;
			}
		}

		double GAMMA = 5. * (1.1 * this->abundances[1] + this->abundances[0] + this->abundances[5]) / (3. * (1.1 * this->abundances[1] + this->abundances[0] + this->abundances[5]));
		GAMMA = 5./3.;
		mu = 2.34;
		eng = fitTEMP(NUMDENS_CGS) * PARAM::KBOLTZ_cgs * (1.1 + this->abundances[0] - this->abundances[5]) / ((GAMMA - 1.0) * mu * PARAM::PROTONMASS_CGS) / PARAM::SEng_per_Mass;
//		std::cout << NUMDENS_CGS << "hhh" << eng << std::endl;
		pres = (GAMMA - 1.0) * dens * eng;
#endif
	}

#ifdef OUTPUT_VTK
	void writeAscii_pos(FILE* fp) const {

//		fprintf(fp, "%lf\t%lf\t%lf\t\n", this->pos.x*PARAM::SL/PARAM::pc_cgs,
//				this->pos.y*PARAM::SL/PARAM::pc_cgs, this->pos.z*PARAM::SL/PARAM::pc_cgs);
//
//		fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", this->id, this->mass, this->pos.x,
//				this->pos.y, this->pos.z, sqrt(pow(this->pos.x, 2) + pow(this->pos.y, 2) + pow(this->pos.z, 2)),
//				this->vel.x, this->vel.y, this->vel.z, this->dens, this->eng, this->pres);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		fprintf(fp, "%lf\t%lf\n", this->pos.x,
				this->pos.y);
#else
		fprintf(fp, "%lf\t%lf\t%lf\t\n", this->pos.x,
				this->pos.y, this->pos.z);
#endif
	}
	void writeAscii_vectors(const char* name,FILE* fp) const {
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		if(name == "vel") {
			fprintf(fp, "%lf\t%lf\t%lf\t", this->vel.x*PARAM::SVel/1e5, this->vel.y*PARAM::SVel/1e5);
		}
#else
		if(std::strcmp(name, "vel")==0 ) {
			fprintf(fp, "%lf\t%lf\t%lf\t\t", this->vel.x*PARAM::SVel/1e5, this->vel.y*PARAM::SVel/1e5,this->vel.z*PARAM::SVel/1e5);
		}
		if(std::strcmp(name, "acc")==0 ) {
			fprintf(fp, "%lf\t%lf\t%lf\t\t", this->acc.x, this->acc.y,this->acc.z);
		}
		if(std::strcmp(name, "pos")==0 ) {
			fprintf(fp, "%lf\t%lf\t%lf\t\t", this->pos.x, this->pos.y,this->pos.z);
		}

#endif
	}
//scalar name in char
	void writeAscii_scalars(const char* name,FILE* fp) const {
		//Convert them to CGS

//	std::cout << "returned temp :" <<" "<<name<< std::endl;
//	PS::F64 tem = this->pres / ( PARAM::SKB*NUMDENS_CGS *PARAM::SL*PARAM::SL*PARAM::SL*(1.1+this->abundances[0]-2.*this->abundances[5]));
//		std::cout<<this->pres<<std::endl;
		double mu = this->abundances[1]+this->abundances[7] + 2.0 * this->abundances[5] + 4.0 *this->abundances[2];
		mu /= (this->abundances[1]+ this->abundances[7] +this->abundances[5] + this->abundances[2] + this->abundances[0]);
		PS::F64 NUMDENS_CGS =this-> dens * PARAM::SMassDens / (mu*PARAM::PROTONMASS);
		PS::F64 tem = mu * PARAM::PROTONMASS_CGS * this->pres * PARAM::SEng_per_Mass / (this->dens * PARAM::KBOLTZ_cgs * (1.1 + PARAM::i_abundance_e - PARAM::i_abundance_H2));

		if(std::strcmp(name , "mass")==0) {
			fprintf(fp, "%lf\n", this->mass*PARAM::SM/PARAM::M_SUN_cgs);
		}
		if(std::strcmp(name,"dens")==0) {

			fprintf(fp, "%lf\n", NUMDENS_CGS );
		}
		if(std::strcmp(name,"pres")==0) {
			fprintf(fp, "%lf\n", this->pres);
		}
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		if(name == "vx") {
			fprintf(fp, "%lf\n", this->vel.x*PARAM::SVel/1e5);
		}
		if(name == "vy") {
			fprintf(fp, "%lf\n", this->vel.y*PARAM::SVel/1e5);
		}
#else
		if(std::strcmp(name, "vx")==0 ) {
			fprintf(fp, "%lf\n", this->vel.x*PARAM::SVel/1e5);
		}
		if(std::strcmp(name, "vy")==0 ) {
			fprintf(fp, "%lf\n", this->vel.y*PARAM::SVel/1e5);
		}if(std::strcmp(name, "vz")==0) {
			fprintf(fp, "%lf\n", this->vel.z*PARAM::SVel/1e5);
		}
#endif
		if(std::strcmp(name, "T")==0 ) {
			fprintf(fp, "%lf\n", tem);
		}if(std::strcmp(name, "eng")==0 ) {
			fprintf(fp, "%lf\n", this->eng);
		}

		if(std::strcmp(name, "ab_e")==0 ) {
			fprintf(fp, "%lf\n", this->abundances[0]);
		}
		if(std::strcmp(name, "ab_HI")==0) {
			fprintf(fp, "%lf\n", this->abundances[1]);
		}
		if(std::strcmp(name, "ab_HeI")==0) {
			fprintf(fp, "%lf\n", this->abundances[2]);
		}
		if(std::strcmp(name, "ab_OI")==0) {
			fprintf(fp, "%lf\n", this->abundances[3]);
		}
		if(std::strcmp(name, "ab_CI")==0) {
			fprintf(fp, "%lf\n", this->abundances[4]);
		}
		if(std::strcmp(name, "ab_H2")==0) {
			fprintf(fp, "%lf\n", this->abundances[5]);
		}
		if(std::strcmp(name,"ab_CO")==0) {
			fprintf(fp, "%lf\n", this->abundances[6]);
		}
		if(std::strcmp(name, "ab_HII")==0) {
			fprintf(fp, "%lf\n", this->abundances[7]);
		}
		if(std::strcmp(name, "ab_CII")==0) {
			fprintf(fp, "%lf\n", this->abundances[8]);
		}
		if(std::strcmp(name, "ab_FeII")==0) {
			fprintf(fp, "%lf\n", this->abundances[9]);
		}
		if(std::strcmp(name,"ab_SiII")==0) {
			fprintf(fp, "%lf\n", this->abundances[9]);
		}

		if(std::strcmp(name,"cooling")==0) {
			fprintf(fp, "%lf\n", this->Lambda);
		}
		if(std::strcmp(name, "heating")==0) {
			fprintf(fp, "%lf\n",this->snds);
		}

		if(std::strcmp(name,"cooling_timescale")==0) {
			fprintf(fp, "%lf\n", this->cooling_timescale);
		}

		//		fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", this->id, this->mass, this->pos.x,
		//				this->pos.y, this->pos.z, sqrt(pow(this->pos.x, 2) + pow(this->pos.y, 2) + pow(this->pos.z, 2)),
		//				this->vel.x, this->vel.y, this->vel.z, this->dens, this->eng, this->pres);
	}

	void writeAsciiTxt(FILE* fp) const {
		PS::F64 massdens = this->dens* PARAM::SMassDens;
		PS::F64 numdens = this->dens* PARAM::SMassDens/PARAM::PROTONMASS;
		PS::F64 tem = this->dens;
//		PS::F64 tem = this->eng * this->dens * (PARAM::GAMMA - 1.0) / (numdens * PARAM::sNumDens * PARAM::sKB);
		//		std::cout <<tem << std::endl;
		//		PS::F64 tem = this->eng;
//		fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", this->id,
//				this->mass * PARAM::SM / PARAM::M_SUN_cgs, this->pos.x * PARAM::SL / PARAM::pc_cgs,
//				this->pos.y * PARAM::SL / PARAM::pc_cgs, this->pos.z * PARAM::SL / PARAM::pc_cgs,
//				sqrt(pow(this->pos.x, 2) + pow(this->pos.y, 2) + pow(this->pos.z, 2)) * PARAM::SL / PARAM::pc_cgs,
//				this->vel.x * PARAM::sVel / 1e6, this->vel.y * PARAM::sVel / 1e6, this->vel.z * PARAM::sVel / 1e6,
//				this->snds * PARAM::sVel / 1e6, tem, this->pres);

		PS::F64 radial_mag = sqrt(this->pos * this->pos);
		PS::F64vec radial = this->pos / radial_mag;
		PS::F64 acc_mag = PARAM::GAMMA * this->pres / this->dens;
		PS::F64 radialForce_mag_ref =this-> gradP *radial;
		//		PS::F64vec perpForce = acc*radial - radialForce;=
		//		PS::F64 radialForce_mag =acc_mag*(this->graddens.x*radial_mag/this->pos.x + this->graddens.y*radial_mag/this->pos.y +this->graddens.z*radial_mag/this->pos.z );
		PS::F64 radialForce_mag = (this->acc ) * radial;
		PS::F64 vol = pow(radial_mag, 3) * M_PI * 4. / 3.;

		PS::F64 radialForce_grav_mag =this->grav*radial;
//		PS::F64 radialForce_mag_ref = -(cos(radial_mag) / radial_mag - sin(radial_mag) / pow(radial_mag, 2));
		//		PS::F64 perpForce_mag = sqrt(perpForce * perpForce);
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", this->id, this->mass,
				this->pos.x, this->pos.y, radial_mag, this->vel.x, this->vel.y, this->dens,
				this->eng, this->pres, radialForce_mag, radialForce_grav_mag, radialForce_mag_ref,this->smth);
#else
		fprintf(fp, "%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", this->id, this->mass,
				this->pos.x, this->pos.y, this->pos.z, radial_mag, this->vel.x, this->vel.y, this->vel.z, this->dens,
				this->eng, this->pres, radialForce_mag, radialForce_grav_mag, radialForce_mag_ref,this->smth);
#endif

	}
#endif

	void setPressure() {
		double GAMMA = 5.
				* (1.1 * this->abundances[1] + this->abundances[0]
						+ this->abundances[5])
				/ (3.
						* (1.1 * this->abundances[1] + this->abundances[0]
								+ this->abundances[5]));

//		if (eng < 0.0 || pres < 0.0) {
//			std::cout << pres << "in class.h " << eng << " " << dens << std::endl;
//			double NUMDENS_CGS =dens * PARAM::SMassDens / (2.34 * PARAM::PROTONMASS);
//
//			double T = fitTEMP(NUMDENS_CGS);
//			double abundance_e = fitAB_E(NUMDENS_CGS);
//			double abundance_H2 = fitAB_H2(NUMDENS_CGS);
//			double energy = T * (PARAM::KBOLTZ_CGS * (1.1 + abundance_e - abundance_H2)) / (2.34 * PARAM::PROTONMASS_CGS * (5. / 3. - 1.0));
//			eng = energy / PARAM::SEng_per_Mass;
//
//		}

		//PS::F64 hcr = 1.4;//heat capacity ratio
		PS::F64 hcr = 5. / 3.; //heat capacity ratio

//		std::cout << "press"<<eng - .5 * vel * vel - .5 * MagneticB * MagneticB / dens << std::endl;
		pres = (hcr - 1.0) * dens
				* (eng - .5 * vel * vel - .5 * MagneticB * MagneticB / dens);
		if (pres < 0.0) {
//			dt = 0.1*fabs( (eng -(.5 * vel * vel +.5 * MagneticB * MagneticB) / dens)/eng_dot);
//					std::cout << "press"<<(eng - .5 * vel * vel - .5 * MagneticB * MagneticB / dens)/eng_dot << std::endl;
			pres = 1.e-4;

		}

		snds = sqrt((hcr * pres) / dens);
	}

};

namespace EPI {
class Dens {
public:
	PS::F64vec pos;
	PS::F64 mass;
	PS::F64 smth;
	void copyFromFP(const RealPtcl& rp) {
		this->pos = rp.pos;
		this->mass = rp.mass;
		this->smth = rp.smth;
	}
	PS::F64vec getPos() const {
		return this->pos;
	}
	PS::F64 getRSearch() const {
		return kernel_t::supportRadius() * this->smth;
	}
	void setPos(const PS::F64vec& pos) {
		this->pos = pos;
	}
};
class Drvt {
public:
	PS::F64vec pos;
	PS::F64vec vel;
	PS::F64 smth;
	PS::F64 dens;
	PS::F64 mass;
	PS::F64 pres;
	PS::F64vec MagneticB;
	void copyFromFP(const RealPtcl& rp) {
		this->MagneticB = rp.MagneticB;
		this->pres = rp.pres;
		pos = rp.pos;
		vel = rp.vel;
		dens = rp.dens;
		smth = rp.smth;
	}
	PS::F64vec getPos() const {
		return this->pos;
	}
	PS::F64 getRSearch() const {
		return kernel_t::supportRadius() * this->smth;
	}
};
class Hydro {
public:
	PS::F64vec pos;
	PS::F64vec vel;
	PS::F64 smth;
	PS::F64 dens;
	PS::F64 pres;
	PS::F64 snds;
	PS::F64 Bal;
	PS::F64vec grad_dens;
	PS::F64vec gradP;
	PS::F64vec gradPT;
	PS::F64vec gradVpara;
	PS::F64vec gradBperp2;
	PS::F64vec gradvperp_x, gradvperp_y, gradvperp_z;
	PS::F64vec gradBperp_x, gradBperp_y, gradBperp_z;
	PS::F64vec grav;
	PS::F64vec gradV;
	PS::F64vec graddens;
	PS::F64vec gradpres;
	PS::F64vec gradvel_x;
	PS::F64vec gradvel_y;
	PS::F64vec gradvel_z;
	PS::F64vec gradvel;
	PS::F64vec xstar;
	PS::F64 dt;
	PS::F64vec MagneticB;
	PS::F64vec vel_half;
	PS::F64vec extF;

	PS::S64 id; ///DEBUG
	void copyFromFP(const RealPtcl& rp) {
		this->extF = rp.extF;
		this->grav = rp.grav;
		this->gradPT = rp.gradPT;
		this->gradBperp2 = rp.gradBperp2;
		this->gradVpara = rp.gradVpara;
		this->dt = rp.dt;
		this->vel_half = rp.vel_half;
		this->MagneticB = rp.MagneticB;
		this->pos = rp.pos;
		this->vel = rp.vel;
		this->smth = rp.smth;
		this->dens = rp.dens;
		this->pres = rp.pres;
		this->snds = rp.snds;
		this->Bal = rp.Bal;
		this->id = rp.id; ///DEBUG
		this->gradBperp_x = rp.gradBperp_x;
		this->gradBperp_y = rp.gradBperp_y;
		this->gradBperp_z = rp.gradBperp_z;
		this->gradvperp_x = rp.gradvperp_x;
		this->gradvperp_y = rp.gradvperp_y;
		this->gradvperp_z = rp.gradvperp_z;
		this->gradP = rp.gradP;
		this->grad_dens = rp.grad_dens;
		this->gradV = rp.gradV;
		this->graddens = rp.graddens;
		this->gradpres = rp.gradpres;
		this->gradvel_x = rp.gradvel_x;
		this->gradvel_y = rp.gradvel_y;
		this->gradvel_z = rp.gradvel_z;

	}
	PS::F64vec getPos() const {
		return this->pos;
	}
	PS::F64 getRSearch() const {
		return kernel_t::supportRadius() * this->smth;
	}
};
class Grav {
public:
	PS::F64vec pos;
	PS::F64 eps2;
	PS::S64 id;
	PS::F64 smth;

	PS::F64vec getPos() const {
		return this->pos;
	}
	PS::F64 getEps2(void) const {
		return 1.0e-4;
	}
	void copyFromFP(const RealPtcl& rp) {
		pos = rp.pos;
		id = rp.id;
		this->smth = rp.smth;
	}
};
}

namespace EPJ {
class Dens {
public:
	PS::F64 mass;
	PS::F64vec pos;
	void copyFromFP(const RealPtcl& rp) {
		this->mass = rp.mass;
		this->pos = rp.pos;
	}
	PS::F64vec getPos() const {
		return this->pos;
	}
	void setPos(const PS::F64vec& pos) {
		this->pos = pos;
	}
};
class Drvt {
public:
	PS::F64 mass;
	PS::F64 dens;
	PS::F64 pres;
	PS::F64vec MagneticB;
	PS::F64vec pos;
	PS::F64vec vel;

	void copyFromFP(const RealPtcl& rp) {
		this->MagneticB = rp.MagneticB;
		this->pres = rp.pres;
		this->dens = rp.dens;
		this->mass = rp.mass;
		this->pos = rp.pos;
		this->vel = rp.vel;
	}
	PS::F64vec getPos() const {
		return this->pos;
	}
	void setPos(const PS::F64vec& pos) {
		this->pos = pos;
	}
};
class Hydro {
public:
	PS::F64 mass;
	PS::F64vec pos;
	PS::F64vec vel;
	PS::F64 smth;
	PS::F64 dens;
	PS::F64 pres;
	PS::F64 snds;
	PS::F64 Bal;
	PS::F64vec grad_dens;
	PS::F64vec grav;
	PS::F64vec gradP;
	PS::F64vec gradPT;
	PS::F64vec gradBperp2;
	PS::F64vec gradvperp_x, gradvperp_y, gradvperp_z;
	PS::F64vec gradBperp_x, gradBperp_y, gradBperp_z;
	PS::F64vec gradVpara;
	PS::F64vec gradV;
	PS::F64vec graddens;
	PS::F64vec gradpres;
	PS::F64vec gradvel_x;
	PS::F64vec gradvel_y;
	PS::F64vec gradvel_z;
	PS::F64vec gradvel;
	PS::F64vec xstar;
	PS::F64vec MagneticB;
	PS::F64vec vel_half;
	PS::F64vec extF;

	PS::S64 id; ///DEBUG
	void copyFromFP(const RealPtcl& rp) {
		this->extF = rp.extF;
		this->grav = rp.grav;
		this->gradPT = rp.gradPT;
		this->gradBperp2 = rp.gradBperp2;
		this->gradVpara = rp.gradVpara;
		this->mass = rp.mass;
		this->vel_half = rp.vel_half;
		this->MagneticB = rp.MagneticB;
		this->pos = rp.pos;
		this->vel = rp.vel;
		this->smth = rp.smth;
		this->dens = rp.dens;
		this->pres = rp.pres;
		this->snds = rp.snds;
		this->Bal = rp.Bal;
		this->id = rp.id; ///DEBUG
		this->gradBperp_x = rp.gradBperp_x;
		this->gradBperp_y = rp.gradBperp_y;
		this->gradBperp_z = rp.gradBperp_z;
		this->gradvperp_x = rp.gradvperp_x;
		this->gradvperp_y = rp.gradvperp_y;
		this->gradvperp_z = rp.gradvperp_z;
		this->gradP = rp.gradP;
		this->grad_dens = rp.grad_dens;
		this->gradV = rp.gradV;
		this->graddens = rp.graddens;
		this->gradpres = rp.gradpres;
		this->gradvel_x = rp.gradvel_x;
		this->gradvel_y = rp.gradvel_y;
		this->gradvel_z = rp.gradvel_z;

	}
	PS::F64vec getPos() const {
		return this->pos;
	}
	PS::F64 getRSearch() const {
		return kernel_t::supportRadius() * this->smth;
	}
	void setPos(const PS::F64vec& pos) {
		this->pos = pos;
	}
};
class Grav {
public:
	PS::F64vec pos;
	PS::F64 mass;
	PS::S64 id;
	PS::F64 smth;
	PS::F64 getEps2(void) const {
		return 1.0e-4;
	}
	PS::F64vec getPos() const {
		return this->pos;
	}
	PS::F64 getCharge(void) const {
		return this->mass;
	}
	void copyFromFP(const RealPtcl& rp) {
		this->mass = rp.mass;
		this->pos = rp.pos;
		this->id = rp.id;
		this->smth = rp.smth;
	}
};
}

namespace SPJ {

}

