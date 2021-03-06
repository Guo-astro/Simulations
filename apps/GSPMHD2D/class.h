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
	int readAscii(FILE* fp) {
		fscanf(fp, "%lf\n", &time);
		fscanf(fp, "%d\n", &Nbody);
		return Nbody;
	}
	void writeAscii(FILE* fp) const {
		fprintf(fp, "%e\n", time);
		fprintf(fp, "%d\n", Nbody);
	}
};

namespace RESULT {
//Density summation
class Dens {
public:
	PS::F64 dens;
	PS::F64 smth;
	void clear() {
		dens = smth = 0;
	}
};
//for Balsara switch
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
	}
};
//Hydro force
class Hydro {
public:
	PS::F64vec acc;
	PS::F64 eng_dot;
	PS::F64vec BoverDens_dot;
	PS::F64 dt;
	void clear() {
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

	PS::F64 dens; //DENSity
	PS::F64 eng; //ENerGy
	PS::F64 pres; //PRESsure
	PS::F64 smth; //SMooTHing length
	PS::F64 snds; //SouND Speed
	PS::F64 div_v;
	PS::F64vec rot_v;
	PS::F64 Bal; //Balsala switch
	PS::F64vec grad_dens;
	PS::F64vec gradP;
	PS::F64vec gradPT;
	PS::F64vec gradvperp_x, gradvperp_y;
	PS::F64vec gradBperp_x, gradBperp_y;
	PS::F64vec gradBperp2;
	PS::F64vec gradVpara;
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

	PS::F64vec grav;
	PS::F64 pot;
	//Copy functions
	void copyFromForce(const RESULT::Dens& dens) {
		this->dens = dens.dens;
		this->smth = dens.smth;
	}
	void copyFromForce(const RESULT::Drvt& drvt) {
		this->gradPT = drvt.gradPT;
		this->gradBperp2 = drvt.gradBperp2;
		this->gradVpara = drvt.gradVpara;
		this->gradP = drvt.gradP;
		this->grad_dens = drvt.grad_dens;
		this->gradBperp_x = drvt.gradBperp_x;
		this->gradBperp_y = drvt.gradBperp_y;
		this->gradvperp_x = drvt.gradvperp_x;
		this->gradvperp_y = drvt.gradvperp_y;
		this->div_v = drvt.div_v;
		this->rot_v = drvt.rot_v;
		this->Bal = fabs(drvt.div_v) / (fabs(drvt.div_v) + sqrt(drvt.rot_v * drvt.rot_v) + 1.0e-4 * this->snds / this->smth); //Balsala switch
	}
	void copyFromForce(const RESULT::Hydro& force) {
		this->acc = force.acc;
		this->eng_dot = force.eng_dot;
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
		fprintf(fp, "%lld\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", this->id, this->mass, this->pos.x, this->pos.y, this->vel.x, this->vel.y, this->dens, this->eng, this->pres, this->MagneticB.x, this->MagneticB.y);
	}
	/*
	 void readAscii(FILE* fp){
	 fscanf(fp, "%ld\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->pos.z, &this->vel.x, &this->vel.y, &this->vel.z, &this->dens, &this->eng, &this->pres);
	 }
	 */
	void readAscii(FILE* fp) {
		fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t"
				"%lf\t%lf\t%lf", &this->id, &this->mass, &this->pos.x, &this->pos.y, &this->vel.x, &this->vel.y, &this->dens, &this->eng, &this->pres, &this->MagneticB.x, &this->MagneticB.y);
	}
	void setPressure() {
		//PS::F64 hcr = 1.4;//heat capacity ratio
		PS::F64 hcr = 5. / 3.; //heat capacity ratio
		pres = (hcr - 1.0) * dens * (eng - .5 * vel * vel - .5 * MagneticB * MagneticB / dens);
		snds = sqrt((hcr * pres ) / dens);
		if (pres < 0.0) {
			pres = 1e-4;
			std::cout << pres << std::endl;

		}
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
	PS::F64vec gradvperp_x, gradvperp_y;
	PS::F64vec gradBperp_x, gradBperp_y;

	PS::F64 dt;
	PS::F64vec MagneticB;
	PS::F64vec vel_half;
	PS::S64 id; ///DEBUG
	void copyFromFP(const RealPtcl& rp) {
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
		this->gradvperp_x = rp.gradvperp_x;
		this->gradvperp_y = rp.gradvperp_y;

		this->gradP = rp.gradP;
		this->grad_dens = rp.grad_dens;

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
	PS::F64vec getPos() const {
		return this->pos;
	}
	PS::F64 getEps2(void) const {
		return 1.0e-4;
	}
	void copyFromFP(const RealPtcl& rp) {
		pos = rp.pos;
		id = rp.id;
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
	PS::F64vec gradP;
	PS::F64vec gradPT;
	PS::F64vec gradBperp2;
	PS::F64vec gradvperp_x, gradvperp_y;
	PS::F64vec gradBperp_x, gradBperp_y;
	PS::F64vec gradVpara;
	PS::F64vec MagneticB;
	PS::F64vec vel_half;
	PS::S64 id; ///DEBUG
	void copyFromFP(const RealPtcl& rp) {
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
		this->gradvperp_x = rp.gradvperp_x;
		this->gradvperp_y = rp.gradvperp_y;
		this->gradP = rp.gradP;
		this->grad_dens = rp.grad_dens;

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
	}
};
}

namespace SPJ {

}

