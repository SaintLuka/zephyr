#include <zephyr/phys/tests/swe/step.h>

namespace zephyr::phys::swe {

inline double r1(double h_r, double h_l, double u_l) {
	return u_l + 2.0 * (std::sqrt(GRAV * h_l) - std::sqrt(GRAV * h_r));
}

inline double spshock(double h_l, double u_l, double h_r) {
	return u_l + h_r * std::sqrt(GRAV * (h_l + h_r) / (2 * h_l * h_r));
}

Step::Step() {
	xdam = 0.0;

	// Case where: h_l > h_r + step_size and u_l = u_r = 0
	// WARNING: Since this solution required calculations we didn't implement in SWASHES for simplicity's sake,
	// the values of h_left, u_left, h_right and u_right must not be changed.
	// See the function compute for more details.
	h_left = 4.;
	h_right = 1.;
	u_left = 0.;
	u_right = 0.;
	step_size = 1.;

	//h_1 and h_2 are linked to h_L and h_R respectively by a rarefaction wave for h_1 and a shock wave for h_2.
	//Using the locus of those 2 waves this gives us a system of equation to find the points compatible
	//with the ideal steady step transition.
	//Those heights can be found using solving methods such as Newton-Raphson.
	//We didn't implement it in this code to keep it light.
	//Here is a sagemath code one can use to find the values of h_1 and h_2 in other cases:
	/*	reset()
		hr = 1
		hl = 4
		g = 9.81
		z = 1
		var('h1,h2')
		solve([ h1*(2*sqrt(g*hl)-2*sqrt(g*h1))-h2*sqrt(g/2*(h2-hr)^2*(1/h2+1/hr))==0, (2*sqrt(g*hl)-2*sqrt(g*h1))^2/2+g*h1-g/4*(h2-hr)^2*(1/h2+1/hr)-g*(h2+z)==0], h1,h2)*/
	h_1 = 3.0922845916198830698;
	h_2 = 1.8999144476020387087;

	//Since the state U_1 is in the rarefaction wave (in the first characteristic field) locus, we use it to find u_1.
	u_1 = r1(h_1, h_left, 0.0);

	//There is a constant discharge at the step: h_1*u_1=h_c*u_c
	u_2 = h_1 * u_1 / h_2;
}

double Step::bed(double x) const {
	return x < xdam ? 0.0 : step_size;
}

double Step::depth(double x, double t) const {
	if (x <= xdam) {
		// we check on which side of the dam we are
		// water height of the left side rarefaction wave:
		double h_wave = std::pow(2 * std::sqrt(GRAV * h_left) - (x - xdam) / t, 2) / (9 * GRAV);
		if (h_wave >= h_left) {
			return h_left;
		}
		else if (h_wave >= h_1) {
			return h_wave;
		}
		else {
			return h_1;
		}
	}
	else { //xex[i]>xdam
		if (x >= xdam + t * spshock(h_2, u_2, h_right)) {
			return h_right;
		}
		else {
			return h_2;
		}
	}
}

double Step::speed(double x, double t) const {
	if (x <= xdam) {
		// we check on which side of the dam we are
		double h_wave = std::pow(2 * std::sqrt(GRAV * h_left) - (x - xdam) / t, 2) / (9. * GRAV);
		if (h_wave >= h_left) {
			return 0.0;
		}
		else if (h_wave >= h_1) {
			return r1(h_wave, h_left, 0.0);
		}
		else {
			return u_1;
		}
	}
	else {
		if (x >= xdam + t * spshock(h_2, u_2, h_right)) {
			return 0.0;
		}
		else {
			return u_2;
		}
	}
}

IBed::Ptr Step::topography() const {
	return geom::StepBed::create(0.0, step_size, xdam);
}

} // zephyr::phys::swe