//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#include<sys/stat.h>

#include "../CloudCloudCollision_Plane_Plane/header.h"

void makeOutputDirectory(const char * dir_name) {
	struct stat st;
	if (stat(dir_name, &st) != 0) {
		PS::S32 ret = -1;
		if (PS::Comm::getRank() == 0)
			ret = mkdir(dir_name, 0777);
		PS::Comm::broadcast(&ret, 1);
		if (ret == 0) {
			if (PS::Comm::getRank() == 0)
				fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
		} else {
			fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
			PS::Abort();
		}
	}
}

int main(int argc, char* argv[]) {
	//////////////////
	//Create vars.
	//////////////////
	PS::Initialize(argc, argv);
	makeOutputDirectory("result");
	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	PS::DomainInfo dinfo;

	boundary box;
	PS::F64 dt, end_time;
	SetupIC_BLASTWAVE(sph_system, &end_time,&box);
//	SetupIC_PolytropeTest_file(sph_system, &end_time);
	dinfo.initialize();

	//////////////////
	//Disp. Info
	//////////////////
	DisplayInfo();
	//////////////////
	//Setup Initial
	//////////////////
	dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);

	Initialize(sph_system);
	//Dom. info
	dinfo.setPosRootDomain(PS::F64vec(-box.x, -box.y, 0.0), PS::F64vec(box.x,box.y,box.z));

	dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
	dinfo.decomposeDomainAll(sph_system);
	sph_system.exchangeParticle(dinfo);

	//plant tree
	PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather dens_tree;
//	PS::TreeForForceShort<RESULT::Drvt, EPI::Drvt, EPJ::Drvt>::Gather drvt_tree;
	PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
	PS::TreeForForceLong <RESULT::Grav , EPI::Grav , EPJ::Grav >::Monopole grav_tree;

	dens_tree.initialize(sph_system.getNumberOfParticleGlobal());
//	drvt_tree.initialize(sph_system.getNumberOfParticleGlobal());
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal());
//	grav_tree.initialize(sph_system.getNumberOfParticleGlobal());

	for (int loop = 0; loop <= 5; ++loop) {
		dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
	}
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		sph_system[i].setPressure();
	}
//	drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
//	grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

	dt = getTimeStepGlobal(sph_system);

	PS::F64 time = 0.0;
	std::cout << std::scientific << std::setprecision(16) << "time = " << time << ", dt = " << dt << "phy_time = " << time * PARAM::ST / PARAM::yr / 1e6 << ", phy_dt = " << dt * PARAM::ST / PARAM::yr / 1e6 << std::endl;

	PS::S32 step = 0;
	for (; time < end_time; time += dt, ++step) {
		InitialKick(sph_system, dt);
		FullDrift(sph_system, dt);
		sph_system.adjustPositionIntoRootDomain(dinfo);
		Predict(sph_system, dt);
		dinfo.decomposeDomainAll(sph_system);
		sph_system.exchangeParticle(dinfo);
		for (int loop = 0; loop <= 5; ++loop) {
			dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
		}
		for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
			sph_system[i].setPressure();
		}
//		drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
		hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
//		grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);

		dt = getTimeStepGlobal(sph_system);

		FinalKick(sph_system, dt,time);
		if (step % PARAM::OUTPUT_INTERVAL == 0) {
			FileHeader header;
			header.time = time;
			header.Nbody = sph_system.getNumberOfParticleGlobal();
			char filename[256];
			sprintf(filename, "result/%04d.dat", step);
			sph_system.writeParticleAscii(filename, header);
			if (PS::Comm::getRank() == 0) {
				std::cout << "//================================" << std::endl;
				std::cout << "output " << filename << "." << std::endl;
				std::cout << "//================================" << std::endl;
			}
		}

		if (PS::Comm::getRank() == 0) {
			std::cout << "//================================" << std::endl;
			std::cout << std::scientific << std::setprecision(16) << "time = " << time << ", dt = " << dt << "phy_time = " << time * PARAM::ST / PARAM::yr / 1e6 << ", phy_dt = " << dt * PARAM::ST / PARAM::yr / 1e6 << std::endl;
			std::cout << "step = " << step << std::endl;
			std::cout << "//================================" << std::endl;
		}
		CheckConservativeVariables(sph_system);
	}

	PS::Finalize();
	return 0;
}

