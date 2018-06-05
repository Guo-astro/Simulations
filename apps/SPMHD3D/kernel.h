#pragma once

struct kernel_t {
	kernel_t() {
	}

	PS::F64 W(const PS::F64vec dr, const PS::F64 h) const {
		PS::F64 value;
		const PS::F64 drnorm = sqrt(dr * dr);
		value = drnorm > 3.1 * h ? 0.0 : pow(h * sqrt(M_PI), -PARAM::Dim) * exp(-(dr * dr) / (h * h));
		return value;
	}
	//gradW
	PS::F64vec gradW(const PS::F64vec dr, const PS::F64 h) const {
		PS::F64vec vec_value;
		vec_value = dr * (pow(h * sqrt(M_PI), -PARAM::Dim) * (-2.0 / (h * h)) * exp(-(dr * dr) / (h * h)));
		return vec_value;
	}
	//gradW
	PS::F64 gradW_scoord(const PS::F64 dsij, const PS::F64 h) const {
		PS::F64 _value;
		//gradient of gaussian kernel
		_value = dsij * (pow(h * sqrt(M_PI), -PARAM::Dim) * (-2.0 / (h * h)) * exp(-(dsij * dsij) / (h * h)));
		return _value;
	}
	PS::F64 DW_h(const PS::F64vec dr, const PS::F64 h) const {
		PS::F64 value;
		//gradient of gaussian kernel
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
		value = (-2 * h * h + 2 * dr * dr) * exp(-(dr * dr) / (h * h)) / (pow(h, 5) * M_PI);
#else
		value = (-3 * h * h + 2 * dr * dr) * exp(-(dr * dr) / (h * h)) / (pow(h, 6.) * pow(M_PI, 1.5));
#endif
		return value;

	}
	static PS::F64 supportRadius() {
		return 3.0;
	}
};
