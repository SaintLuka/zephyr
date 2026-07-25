#include <cmath>
#include <iostream>
#include <zephyr/phys/tests/swe/dam_break.h>

namespace zephyr::phys {

constexpr double GRAV = 9.81;

DamBreak::DamBreak(double depth_L, double depth_R, double x_dam) {
	h_L = depth_L;
	h_R = depth_R;
	x0 = x_dam;

	compute();
}

/// @details Function to solve by dichotomy the equation
/// \f$ cm^6-9v_{right}^2cm^4+16v_{left}\, v_{right}^2cm^3-v_{right}^2(v_{right}^2+8v_{left}^2)cm^2+v_{right}^6=0\f$.
/// @return Value of \f$x^6-9v_{right}^2x^4+16v_{left}\, v_{right}^2x^3-v_{right}^2(v_{right}^2+8v_{left}^2)x^2+v_{right}^6\f$.
inline double function(double x, double v_left, double v_right) {
	return pow(x,6.)-9.*pow(v_right,2.)*pow(x,4.)+16.*v_left*pow(v_right,2.)*pow(x,3.)-pow(v_right,2.)*(pow(v_right,2.)+8.*pow(v_left,2.))*pow(x,2.)+pow(v_right,6.);
}


void DamBreak::compute() {	
	/// @details
	/// Computes the chosen dam break solution, see \cite Ritter92 \cite Stoker57.
	/// @par Modifies
	/// Solution#hex.

	c_L = std::sqrt(GRAV * h_L); // cl left wave velocity
	c_R = std::sqrt(GRAV * h_R); // cr right wave velocity

	double x_a{0.0}, x_b{0.0}; // variables for dichotomy
		
	// Parameters for the dichotomy
	double eps = 1.0e-6;
	int nmax = 1000;

	double func = function(c_L, c_L, c_R);
	std::cout << func << std::endl;
	if (func < 0.0){
		x_a = c_L; // func(cl)<0
	}else{
		x_b = c_L; // func(cl)>0
	}//end if

	func = function(c_R, c_L, c_R);
	std::cout << func << std::endl;
	if (func < 0.0){
		x_a = c_R; // func(cr)<0
	} else {
		x_b = c_R; // func(cr)>0
	}

	/* dichotomy in order to solve the equation in cm
	 * cm^6-9*cr^2*cm^4+16*cl*cr^2*cm^3-cr^2*(cr^2+8*cl^2)*cm^2+cr^6=0
	 * in order to get the water height hm (when hr is not null)
	 */
	int iter = 0;
	while(std::fabs(x_a - x_b) > eps && iter < nmax) {
		double mid = 0.5 * (x_a + x_b);
		std::cout << x_a << " " << mid << " " << x_b << "\n";

		func = function(mid, c_L, c_R);

		if (func < 0.0){
			x_a = mid;
		}else{
			x_b = mid;
		}
		++iter;
	}

	c_mid = 0.5 * (x_a + x_b); //cm

	if (std::abs(h_R) < 1.0e-6) {
		c_mid = 0.0;
	} // the dam break on dry soil

	h_mid = c_mid*c_mid/GRAV; //the water height hm
	u_mid = 2.0*(c_L-c_mid); //the velocity um
	v = h_mid*u_mid/(h_mid-h_R); //the velocity of the shock

	std::cout << h_L << " " << h_mid << " " << h_R << "\n";
	std::cout << u_mid << " " << v << "\n";
}

double DamBreak::depth(double x, double t) const {
	if (x <= x0 - c_L * t){
		return h_L;
	}
	if (x <= x0 + (2.0 * c_L - 3.0 * c_mid) * t) {
		double xi = (x - x0) / t;
		return (4.0/(9.0*GRAV)) * (c_L * (c_L - xi) + std::pow(0.5 * xi, 2));
	}
	if (x <= x0 + v * t){
		return h_mid;
	}
	return h_R;
}

double DamBreak::speed(double x, double t) const {
	if (x <= x0 - c_L * t) {
		return 0.0;
	}
	if (x <= x0 + (2.0 * c_L - 3.0 * c_mid) * t) {
		return (2.0 / 3.0) * ((x - x0) / t + c_L);
	}
	if (x <= x0 + v * t){
		return u_mid;
	}
	return 0.0;
}

} // namespace zephyr::phys