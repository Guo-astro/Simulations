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
//		pres = (GAMMA - 1.0) * dens * eng;
		pres = .4* pow(dens,GAMMA);
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

		pres = (PARAM::GAMMA - 1) * dens * eng;

//		pres = .4 * pow(dens, PARAM::GAMMA);
		snds = sqrt(PARAM::GAMMA * pres / dens);
	}
public:
	double fitTEMP(double numdens) {
		double x = log10(numdens);
		double fit_T_NUMDENS_a = 3.626416807912980E+00;
		double fit_T_NUMDENS_b = -1.224933568952209E+00;
		double fit_T_NUMDENS_c = -1.274578242438200E+00;
		double fit_T_NUMDENS_d = 5.219886021076703E-01;
		double fit_T_NUMDENS_e = 1.053208709509926E+00;
		double fit_T_NUMDENS_f = -3.879459648821958E-01;
		double fit_T_NUMDENS_g = -4.660615655301823E-01;
		double fit_T_NUMDENS_h = 2.356638469442452E-01;
		double fit_T_NUMDENS_i = 8.724038787337413E-02;
		double fit_T_NUMDENS_j = -7.718373457853663E-02;
		double fit_T_NUMDENS_k = 4.013294530466450E-03;
		double fit_T_NUMDENS_l = 1.047836002918579E-02;
		double fit_T_NUMDENS_m = -3.511934998554032E-03;
		double fit_T_NUMDENS_n = 3.528079479227492E-05;
		double fit_T_NUMDENS_o = 2.594946122742665E-04;
		double fit_T_NUMDENS_p = -8.283047406203526E-05;
		double fit_T_NUMDENS_q = 1.383358597907899E-05;
		double fit_T_NUMDENS_r = -1.423434735457324E-06;
		double fit_T_NUMDENS_s = 9.102911931652432E-08;
		double fit_T_NUMDENS_t = -3.340355558070870E-09;
		double fit_T_NUMDENS_u = 5.400936208142200E-11;
		double T_DENS = fit_T_NUMDENS_a + +fit_T_NUMDENS_b * pow(x, 1)
				+ fit_T_NUMDENS_c * pow(x, 2) + fit_T_NUMDENS_d * pow(x, 3)
				+ fit_T_NUMDENS_e * pow(x, 4) + fit_T_NUMDENS_f * pow(x, 5)
				+ fit_T_NUMDENS_g * pow(x, 6) + fit_T_NUMDENS_h * pow(x, 7)
				+ fit_T_NUMDENS_i * pow(x, 8) + fit_T_NUMDENS_j * pow(x, 9)
				+ fit_T_NUMDENS_k * pow(x, 10) + fit_T_NUMDENS_l * pow(x, 11)
				+ fit_T_NUMDENS_m * pow(x, 12) + fit_T_NUMDENS_n * pow(x, 13)
				+ fit_T_NUMDENS_o * pow(x, 14) + fit_T_NUMDENS_p * pow(x, 15)
				+ fit_T_NUMDENS_q * pow(x, 16) + fit_T_NUMDENS_r * pow(x, 17)
				+ fit_T_NUMDENS_s * pow(x, 18) + fit_T_NUMDENS_t * pow(x, 19)
				+ fit_T_NUMDENS_u * pow(x, 20);
		return pow(10.0, T_DENS);
	}

	double fitAB_E(double numdens) {
		double x = log10(numdens);
		double fit_AB_E_NUMDENS_a = -1.872011804950556E+00;
		double fit_AB_E_NUMDENS_b = -8.937317358549671E-01;
		double fit_AB_E_NUMDENS_c = -4.049467520228735E-01;
		double fit_AB_E_NUMDENS_d = 2.340570540550370E-01;
		double fit_AB_E_NUMDENS_e = 3.399024284192451E-01;
		double fit_AB_E_NUMDENS_f = -1.688149498189914E-01;
		double fit_AB_E_NUMDENS_g = -1.423033769865771E-01;
		double fit_AB_E_NUMDENS_h = 9.325026512578551E-02;
		double fit_AB_E_NUMDENS_i = 2.200234025037475E-02;
		double fit_AB_E_NUMDENS_j = -2.777434952732058E-02;
		double fit_AB_E_NUMDENS_k = 3.150027423127273E-03;
		double fit_AB_E_NUMDENS_l = 3.351038047856944E-03;
		double fit_AB_E_NUMDENS_m = -1.344002678105646E-03;
		double fit_AB_E_NUMDENS_n = 6.892793938024374E-05;
		double fit_AB_E_NUMDENS_o = 8.392156974890075E-05;
		double fit_AB_E_NUMDENS_p = -2.999536252357328E-05;
		double fit_AB_E_NUMDENS_q = 5.300810713536076E-06;
		double fit_AB_E_NUMDENS_r = -5.689318414940446E-07;
		double fit_AB_E_NUMDENS_s = 3.771620461037251E-08;
		double fit_AB_E_NUMDENS_t = -1.429771105375033E-09;
		double fit_AB_E_NUMDENS_u = 2.382784941187684E-11;
		double AB_E_DENS = fit_AB_E_NUMDENS_a + +fit_AB_E_NUMDENS_b * pow(x, 1)
				+ fit_AB_E_NUMDENS_c * pow(x, 2)
				+ fit_AB_E_NUMDENS_d * pow(x, 3)
				+ fit_AB_E_NUMDENS_e * pow(x, 4)
				+ fit_AB_E_NUMDENS_f * pow(x, 5)
				+ fit_AB_E_NUMDENS_g * pow(x, 6)
				+ fit_AB_E_NUMDENS_h * pow(x, 7)
				+ fit_AB_E_NUMDENS_i * pow(x, 8)
				+ fit_AB_E_NUMDENS_j * pow(x, 9)
				+ fit_AB_E_NUMDENS_k * pow(x, 10)
				+ fit_AB_E_NUMDENS_l * pow(x, 11)
				+ fit_AB_E_NUMDENS_m * pow(x, 12)
				+ fit_AB_E_NUMDENS_n * pow(x, 13)
				+ fit_AB_E_NUMDENS_o * pow(x, 14)
				+ fit_AB_E_NUMDENS_p * pow(x, 15)
				+ fit_AB_E_NUMDENS_q * pow(x, 16)
				+ fit_AB_E_NUMDENS_r * pow(x, 17)
				+ fit_AB_E_NUMDENS_s * pow(x, 18)
				+ fit_AB_E_NUMDENS_t * pow(x, 19)
				+ fit_AB_E_NUMDENS_u * pow(x, 20);
		return pow(10.0, AB_E_DENS);
	}
	double fitAB_H2(double numdens) {
		double x = log10(numdens);
		double fit_AB_E_NUMDENS_a = -7.360098032317123E+00;
		double fit_AB_E_NUMDENS_b = -3.423714872708186E-01;
		double fit_AB_E_NUMDENS_c = -2.137800396209498E-01;
		double fit_AB_E_NUMDENS_d = 1.695543107865461E+00;
		double fit_AB_E_NUMDENS_e = 6.265122212733395E-01;
		double fit_AB_E_NUMDENS_f = -1.382649164788568E+00;
		double fit_AB_E_NUMDENS_g = -3.758096546109677E-01;
		double fit_AB_E_NUMDENS_h = 6.720367129380326E-01;
		double fit_AB_E_NUMDENS_i = 3.241466672648317E-02;
		double fit_AB_E_NUMDENS_j = -1.862316720489950E-01;
		double fit_AB_E_NUMDENS_k = 3.500919192325033E-02;
		double fit_AB_E_NUMDENS_l = 2.190695186294473E-02;
		double fit_AB_E_NUMDENS_m = -1.046778351894664E-02;
		double fit_AB_E_NUMDENS_n = 7.087659697485427E-04;
		double fit_AB_E_NUMDENS_o = 6.495692527397854E-04;
		double fit_AB_E_NUMDENS_p = -2.446944640902132E-04;
		double fit_AB_E_NUMDENS_q = 4.416837349827385E-05;
		double fit_AB_E_NUMDENS_r = -4.790959830696257E-06;
		double fit_AB_E_NUMDENS_s = 3.190719716023346E-07;
		double fit_AB_E_NUMDENS_t = -1.210306966686387E-08;
		double fit_AB_E_NUMDENS_u = 2.012517239806765E-10;
		double AB_H2_DENS = fit_AB_E_NUMDENS_a + +fit_AB_E_NUMDENS_b * pow(x, 1)
				+ fit_AB_E_NUMDENS_c * pow(x, 2)
				+ fit_AB_E_NUMDENS_d * pow(x, 3)
				+ fit_AB_E_NUMDENS_e * pow(x, 4)
				+ fit_AB_E_NUMDENS_f * pow(x, 5)
				+ fit_AB_E_NUMDENS_g * pow(x, 6)
				+ fit_AB_E_NUMDENS_h * pow(x, 7)
				+ fit_AB_E_NUMDENS_i * pow(x, 8)
				+ fit_AB_E_NUMDENS_j * pow(x, 9)
				+ fit_AB_E_NUMDENS_k * pow(x, 10)
				+ fit_AB_E_NUMDENS_l * pow(x, 11)
				+ fit_AB_E_NUMDENS_m * pow(x, 12)
				+ fit_AB_E_NUMDENS_n * pow(x, 13)
				+ fit_AB_E_NUMDENS_o * pow(x, 14)
				+ fit_AB_E_NUMDENS_p * pow(x, 15)
				+ fit_AB_E_NUMDENS_q * pow(x, 16)
				+ fit_AB_E_NUMDENS_r * pow(x, 17)
				+ fit_AB_E_NUMDENS_s * pow(x, 18)
				+ fit_AB_E_NUMDENS_t * pow(x, 19)
				+ fit_AB_E_NUMDENS_u * pow(x, 20);

		return pow(10.0, AB_H2_DENS);
	}

	double fitAB_CO(double numdens) {
		double x = log10(numdens);
		double fit_AB_E_NUMDENS_a = -2.051274717031706E+01;
		double fit_AB_E_NUMDENS_b = 1.656155893322455E+00;
		double fit_AB_E_NUMDENS_c = -2.263458598301474E-01;
		double fit_AB_E_NUMDENS_d = 1.704379890823728E+00;
		double fit_AB_E_NUMDENS_e = 6.540640533501590E-01;
		double fit_AB_E_NUMDENS_f = -1.401277757381409E+00;
		double fit_AB_E_NUMDENS_g = -3.956320154899900E-01;
		double fit_AB_E_NUMDENS_h = 6.889568533047997E-01;
		double fit_AB_E_NUMDENS_i = 3.647545368628379E-02;
		double fit_AB_E_NUMDENS_j = -1.930445428728992E-01;
		double fit_AB_E_NUMDENS_k = 3.602797485651373E-02;
		double fit_AB_E_NUMDENS_l = 2.291844991727394E-02;
		double fit_AB_E_NUMDENS_m = -1.093717655323075E-02;
		double fit_AB_E_NUMDENS_n = 7.384884969774629E-04;
		double fit_AB_E_NUMDENS_o = 6.831606115135190E-04;
		double fit_AB_E_NUMDENS_p = -2.577595889622206E-04;
		double fit_AB_E_NUMDENS_q = 4.662860736266569E-05;
		double fit_AB_E_NUMDENS_r = -5.069249559221394E-06;
		double fit_AB_E_NUMDENS_s = 3.383564575974674E-07;
		double fit_AB_E_NUMDENS_t = -1.286218414053287E-08;
		double fit_AB_E_NUMDENS_u = 2.143159560810884E-10;
		double AB_H2_DENS = fit_AB_E_NUMDENS_a + +fit_AB_E_NUMDENS_b * pow(x, 1)
				+ fit_AB_E_NUMDENS_c * pow(x, 2)
				+ fit_AB_E_NUMDENS_d * pow(x, 3)
				+ fit_AB_E_NUMDENS_e * pow(x, 4)
				+ fit_AB_E_NUMDENS_f * pow(x, 5)
				+ fit_AB_E_NUMDENS_g * pow(x, 6)
				+ fit_AB_E_NUMDENS_h * pow(x, 7)
				+ fit_AB_E_NUMDENS_i * pow(x, 8)
				+ fit_AB_E_NUMDENS_j * pow(x, 9)
				+ fit_AB_E_NUMDENS_k * pow(x, 10)
				+ fit_AB_E_NUMDENS_l * pow(x, 11)
				+ fit_AB_E_NUMDENS_m * pow(x, 12)
				+ fit_AB_E_NUMDENS_n * pow(x, 13)
				+ fit_AB_E_NUMDENS_o * pow(x, 14)
				+ fit_AB_E_NUMDENS_p * pow(x, 15)
				+ fit_AB_E_NUMDENS_q * pow(x, 16)
				+ fit_AB_E_NUMDENS_r * pow(x, 17)
				+ fit_AB_E_NUMDENS_s * pow(x, 18)
				+ fit_AB_E_NUMDENS_t * pow(x, 19)
				+ fit_AB_E_NUMDENS_u * pow(x, 20);

		return pow(10.0, AB_H2_DENS);
	}

	double fitAB_HI(double numdens) {
		double x = log10(numdens);
		double fit_AB_E_NUMDENS_a = -6.137627239621324E-03;
		double fit_AB_E_NUMDENS_b = 1.072591293581314E-02;
		double fit_AB_E_NUMDENS_c = -1.268358087584409E-03;
		double fit_AB_E_NUMDENS_d = 8.298703190565435E-04;
		double fit_AB_E_NUMDENS_e = -9.380209094011755E-03;
		double fit_AB_E_NUMDENS_f = 8.149209759993337E-06;
		double fit_AB_E_NUMDENS_g = 7.699114612721374E-03;
		double fit_AB_E_NUMDENS_h = -1.073911173255449E-03;
		double fit_AB_E_NUMDENS_i = -2.938966080876961E-03;
		double fit_AB_E_NUMDENS_j = 9.529549555644621E-04;
		double fit_AB_E_NUMDENS_k = 4.251472373259463E-04;
		double fit_AB_E_NUMDENS_l = -2.620579798064097E-04;
		double fit_AB_E_NUMDENS_m = 1.330515545506811E-05;
		double fit_AB_E_NUMDENS_n = 2.122109753312525E-05;
		double fit_AB_E_NUMDENS_o = -6.482216455042173E-06;
		double fit_AB_E_NUMDENS_p = 6.283216303757292E-07;
		double fit_AB_E_NUMDENS_q = 5.133851349512566E-08;
		double fit_AB_E_NUMDENS_r = -1.964694860379738E-08;
		double fit_AB_E_NUMDENS_s = 2.148351439547345E-09;
		double fit_AB_E_NUMDENS_t = -1.114467075592349E-10;
		double fit_AB_E_NUMDENS_u = 2.331356984330223E-12;
		double AB_H2_DENS = fit_AB_E_NUMDENS_a + +fit_AB_E_NUMDENS_b * pow(x, 1)
				+ fit_AB_E_NUMDENS_c * pow(x, 2)
				+ fit_AB_E_NUMDENS_d * pow(x, 3)
				+ fit_AB_E_NUMDENS_e * pow(x, 4)
				+ fit_AB_E_NUMDENS_f * pow(x, 5)
				+ fit_AB_E_NUMDENS_g * pow(x, 6)
				+ fit_AB_E_NUMDENS_h * pow(x, 7)
				+ fit_AB_E_NUMDENS_i * pow(x, 8)
				+ fit_AB_E_NUMDENS_j * pow(x, 9)
				+ fit_AB_E_NUMDENS_k * pow(x, 10)
				+ fit_AB_E_NUMDENS_l * pow(x, 11)
				+ fit_AB_E_NUMDENS_m * pow(x, 12)
				+ fit_AB_E_NUMDENS_n * pow(x, 13)
				+ fit_AB_E_NUMDENS_o * pow(x, 14)
				+ fit_AB_E_NUMDENS_p * pow(x, 15)
				+ fit_AB_E_NUMDENS_q * pow(x, 16)
				+ fit_AB_E_NUMDENS_r * pow(x, 17)
				+ fit_AB_E_NUMDENS_s * pow(x, 18)
				+ fit_AB_E_NUMDENS_t * pow(x, 19)
				+ fit_AB_E_NUMDENS_u * pow(x, 20);
		return pow(10.0, AB_H2_DENS);
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

