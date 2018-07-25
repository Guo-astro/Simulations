#pragma once
#include <particle_simulator.hpp>
PS::F64 rate_Lambda(PS::F64 T, PS::F64 NUMDENS_CGS, PS::F64 abundance_HI,
		PS::F64 abundance_e, PS::F64 abundance_H2, PS::F64 abundance_HII,
		PS::F64 abundance_CII, PS::F64 abundance_FeII, PS::F64 abundance_SiII,
		PS::F64 abundance_CI, PS::F64 abundance_OI, PS::F64 abundance_CO,
		PS::F64 dust_to_gas_ratio, PS::F64 T_gr, PS::F64 grain_size);
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "param.h"
#include "mathfunc.h"
#include "kernel.h"
#include "class.h"
#include "force.h"
#include "prototype.h"

