//#define SANITY_CHECK_REALLOCATABLE_ARRAY
#include<sys/stat.h>
#include "header.h"

void makeOutputDirectory(const char * dir_name) {
	struct stat st;
	if (stat(dir_name, &st) != 0) {
		PS::S32 ret = -1;
		if (PS::Comm::getRank() == 0)
			ret = mkdir(dir_name, 0777);
		PS::Comm::broadcast(&ret, 1);
		if (ret == 0) {
			if (PS::Comm::getRank() == 0)
				fprintf(stderr, "Directory \"%s\" is successfully made.\n",
						dir_name);
		} else {
			fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
			PS::Abort();
		}
	}
}

int main(int argc, char* argv[]) {

//	for (int i = 0; i < 1e6; i++) {
//		double T = i;
//		double NUMDENS_CGS = 1e8;
//		double abundance_HI = 0.1;
//		double abundance_e = 0.1;
//		double abundance_H2 = 0.1;
//		double abundance_HII = 0.1;
//		double abundance_CII = 0.1;
//		double abundance_FeII = 0.1;
//		double abundance_SiII = 0.1;
//		double abundance_CI = 0.1;
//		double abundance_OI = 0.1;
//		double abundance_CO = 0.1;
//		double dust_to_gas_ratio = PARAM::dust_to_gas_ratio;
//		double gt = PARAM::Grain_T;
//		double gs = PARAM::grain_size;
//		double lambda = rate_Lambda(T, NUMDENS_CGS, abundance_HI, abundance_e,
//				abundance_H2, abundance_HII, abundance_CII, abundance_FeII,
//				abundance_SiII, abundance_CI, abundance_OI, abundance_CO,
//				dust_to_gas_ratio, gt, gs);
//		std::cout << lambda << std::endl;
//	}

	////////////////
//	Create vars.
	////////////////
	PS::Initialize(argc, argv);
	makeOutputDirectory("result");
	makeOutputDirectory("result_vtk");
	PS::ParticleSystem<RealPtcl> sph_system;
	sph_system.initialize();
	PS::DomainInfo dinfo;
	std::cout << "---" << std::endl;

	boundary box;
	PS::F64 dt, end_time;
//	SetupIC_BLASTWAVE(sph_system, &end_time,&box);
	SetupIC_HVCC_RESTART_file(sph_system, &end_time);
//	SetupIC_HVCC_file(sph_system, &end_time);
	dinfo.initialize();

	//////////////////
	//Disp. Info
	//////////////////
	DisplayInfo();
	//////////////////
	//Setup Initial
	//////////////////
	dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);

	Initialize(sph_system);
	//Dom. info
	dinfo.setDomain(PS::Comm::getNumberOfProc(), 1, 1);
	dinfo.decomposeDomainAll(sph_system);
	sph_system.exchangeParticle(dinfo);

	//plant tree
	PS::TreeForForceShort<RESULT::Dens, EPI::Dens, EPJ::Dens>::Gather dens_tree;
//	PS::TreeForForceShort<RESULT::Drvt, EPI::Drvt, EPJ::Drvt>::Gather drvt_tree;
	PS::TreeForForceShort<RESULT::Hydro, EPI::Hydro, EPJ::Hydro>::Symmetry hydr_tree;
	PS::TreeForForceLong<RESULT::Grav, EPI::Grav, EPJ::Grav>::Monopole grav_tree;

	dens_tree.initialize(sph_system.getNumberOfParticleGlobal());
//	drvt_tree.initialize(sph_system.getNumberOfParticleGlobal());
	hydr_tree.initialize(sph_system.getNumberOfParticleGlobal());
	grav_tree.initialize(sph_system.getNumberOfParticleGlobal());

	for (int loop = 0; loop <= 5; ++loop) {
		dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
	}
	for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {

		sph_system[i].setPressure();
	}
//	drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
	hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
	grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
	PS::F64 time = 0.0;

	dt = getTimeStepGlobal(sph_system,time);
//	dt = 1.3403086854443947e-05;

	std::cout << std::scientific << std::setprecision(16) << "time = " << time << ", dt = " << dt << "phy_time = " << time * PARAM::ST / PARAM::yr / 1e6 << ", phy_dt = " << dt * PARAM::ST / PARAM::yr / 1e6 << std::endl;

	PS::S32 step = 0;

	for (; time < end_time; time += dt, ++step) {



//		InitialKick(sph_system, dt);
//		FullDrift(sph_system, dt);
		sph_system.adjustPositionIntoRootDomain(dinfo);
		std::cout << "//====adjusted==================" << std::endl;

//		Predict(sph_system, dt);
		dinfo.decomposeDomainAll(sph_system);
		std::cout << "//====decomposed==================" << std::endl;

		sph_system.exchangeParticle(dinfo);
		for (int loop = 0; loop <= 5; ++loop) {
			dens_tree.calcForceAllAndWriteBack(CalcDensity(), sph_system, dinfo);
		}
		std::cout << "//====setdens==================" << std::endl;

		for (PS::S32 i = 0; i < sph_system.getNumberOfParticleLocal(); ++i) {
			sph_system[i].setPressure();
		}
		std::cout << "//====setpres==================" << std::endl;

//		drvt_tree.calcForceAllAndWriteBack(CalcDerivative(), sph_system, dinfo);
//		std::cout << "//====drvt==================" << std::endl;

		hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), sph_system, dinfo);
		std::cout << "//====hydra==================" << std::endl;

		grav_tree.calcForceAllAndWriteBack(CalcGravityForce<EPJ::Grav>(), CalcGravityForce<PS::SPJMonopole>(), sph_system, dinfo);
		std::cout << "//====sg==================" << std::endl;

		dt = getTimeStepGlobal(sph_system,time);
		std::cout << "//====dt================"<<dt << std::endl;

		FinalKick(sph_system, dt,time);
		std::cout << "//====kicked==================" << std::endl;

		if (step % PARAM::OUTPUT_INTERVAL == 0) {
			FileHeader header;
			header.time = time;
			header.Nbody = sph_system.getNumberOfParticleGlobal();
			char filename[256];
			sprintf(filename, "result/%04d.dat", step);
			sph_system.writeParticleAscii(filename, header);
#ifdef OUTPUT_VTK
			FileHeader header_vtk;
			header_vtk.phy_time = time * PARAM::ST / PARAM::yr / 1e6;
			header_vtk.Nbody = sph_system.getNumberOfParticleGlobal();
			char filename_vtk[256];
			sprintf(filename_vtk, "result_vtk/HVCC_CHEM_IMBH.vtk.%04d", step/PARAM::OUTPUT_INTERVAL);
			sph_system.writeParticle_VTK_Ascii(filename_vtk, header_vtk);

			if (PS::Comm::getRank() == 0) {
				std::cout << "//================================" << std::endl;
				std::cout << "output " << filename_vtk << "." << std::endl;
				std::cout << "//================================" << std::endl;
			}
#endif
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

