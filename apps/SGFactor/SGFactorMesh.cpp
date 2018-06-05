//============================================================================
// Name        : SGFactorMesh.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cmath>
#include "vector3.hpp"
using namespace std;
double PI = 3.141592654;
double dens(double radius, double h) {
	double rho = exp(-radius * radius / (h * h)) / pow(sqrt(PI) * h, 3);
	return rho;
}
int main() {

	double start1_R = 0;
	double start1_theta = 0;
	double start1_phi = 0;
	double finish1_R = 3.1415;
	double finish1_theta = PI;
	double finish1_phi = 2. * PI;
	double start2_R = 0;
	double start2_theta = 0;
	double start2_phi = 0;
	double finish2_R = 3.1415;
	double finish2_theta = PI;
	double finish2_phi = 2. * PI;

	int n1_R = 20;
	int n1_theta = 20;
	int n1_phi = 20;
	int n2_R = 20;
	int n2_theta = 20;
	int n2_phi = 20;

	double d1_R = (finish1_R - start1_R) / n1_R;
	double d1_theta = (finish1_theta - start1_theta) / n1_theta;
	double d1_phi = (finish1_phi - start1_phi) / n1_phi;
	double d2_R = (finish2_R - start2_R) / n2_R;
	double d2_theta = (finish2_theta - start2_theta) / n2_theta;
	double d2_phi = (finish2_phi - start2_phi) / n2_phi;
	double smth1 = 1;
	double smth2 = 1;
	double mass1 = 0.0;
	double mass2 = 0.0;

//	for (int i1 = 0; i1 < n1_R; i1++) {
//		double r1_i = start1_R + i1 * d1_R;
//		double rho1_ri = dens(r1_i, smth1);
//		for (int j1 = 0; j1 < n1_theta; j1++) {
//			double theta1_i = start1_theta + j1 * d1_theta;
//			for (int k1 = 0; k1 < n1_phi; k1++) {
//				double phi1_i = start1_phi + j1 * d1_phi;
//				Vector3<double> vec_r1i(r1_i * cos(theta1_i), r1_i * sin(theta1_i) * sin(phi1_i), r1_i * sin(theta1_i) * cos(phi1_i));
//				double r1i_2 = vec_r1i * vec_r1i;
//				double geom1_fac = r1i_2 * sin(theta1_i) * d1_R * d1_theta * d1_phi;
//				mass1 += geom1_fac * rho1_ri;
//			}
//		}
//	}
//	cout << mass1 << endl; // prints !!!Hello World!!!
#pragma omp parallel for
	for (int loop = 0; loop < 600; loop += 10) {
		Vector3<double> disp(0.00001 + loop * 0.01, 0, 0);
		Vector3<double> disp_norm_vec = disp / sqrt(disp * disp);
		Vector3<double> F;

		for (int i1 = 0; i1 < n1_R; i1++) {
			double r1_i = start1_R + i1 * d1_R;
			double rho1_ri = dens(r1_i, smth1);
			for (int j1 = 0; j1 < n1_theta; j1++) {
				double theta1_i = start1_theta + j1 * d1_theta;
				for (int k1 = 0; k1 < n1_phi; k1++) {
					double phi1_i = start1_phi + j1 * d1_phi;
					Vector3<double> vec_r1i(r1_i * cos(theta1_i), r1_i * sin(theta1_i) * sin(phi1_i), r1_i * sin(theta1_i) * cos(phi1_i));
					double r1i_2 = vec_r1i * vec_r1i;
					double geom1_fac = r1i_2 * sin(theta1_i) * d1_R * d1_theta * d1_phi;

					////////////////////////////////////////////////////////////////////
					for (int i2 = 0; i2 < n2_R; i2++) {
						double r2_i = start2_R + i2 * d2_R;
						double rho2_ri = dens(r2_i, smth2);
						for (int j2 = 0; j2 < n2_theta; j2++) {
							double theta2_i = start2_theta + j2 * d2_theta;
							for (int k2 = 0; k2 < n2_phi; k2++) {
								double phi2_i = start2_phi + j2 * d2_phi;
								Vector3<double> vec_r2i(r2_i * cos(theta2_i), r2_i * sin(theta2_i) * sin(phi2_i), r2_i * sin(theta2_i) * cos(phi2_i));
								double r2i_2 = vec_r2i * vec_r2i;
								Vector3<double> vec_r12 = vec_r1i - vec_r2i - disp;

								double r12 = sqrt(vec_r12 * vec_r12);
								if (r12 == 0) {
									continue;
								}
								double geom2_fac = r2i_2 * sin(theta2_i) * d2_R * d2_theta * d2_phi;
								F -= geom2_fac * geom1_fac * rho2_ri * rho1_ri * vec_r12 / pow(r12, 3);

							}
						}
					}

				}
			}
		}
		Vector3<double> F_newton = disp / pow(sqrt(disp * disp), 3);

		cout << sqrt(disp * disp) << " " << F * disp_norm_vec / (F_newton * disp_norm_vec) << endl; // prints !!!Hello World!!!
	}
	return 0;
}
