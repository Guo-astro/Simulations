/*
 * fittings.h
 *
 *  Created on: 2018/04/29
 *      Author: guo
 */

#ifndef APPS_SHOCKTUBE_GSPH_FITTINGS_H_
#define APPS_SHOCKTUBE_GSPH_FITTINGS_H_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

double fitLambda(double T);
double fitGamma(double T);
double fitTEMP(double numdens);
double fitAB_E(double numdens);
double fitAB_H2(double numdens);

double fitAB_CO(double numdens);

double fitAB_HI(double numdens);
#endif /* APPS_SHOCKTUBE_GSPH_FITTINGS_H_ */
