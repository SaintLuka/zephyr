#include <zephyr/phys/tests/swe/thacker_1d.h>

namespace zephyr::phys {

Thacker1D::Thacker1D()
	: Thacker1D(0.8, -0.3) {

}

Thacker1D::Thacker1D(double coeff, double bottom)
	: coeff(coeff), bottom(bottom) {
	h0 = -bottom;
	a = std::sqrt(h0 / coeff);
	omega = std::sqrt(2.0 * coeff * GRAV);

	lx = std::ceil(10.0 * (a + 0.5)) / 10.0;
}

double Thacker1D::max_time() const {
	return 2.0 * M_PI / omega;
}

double Thacker1D::bed(double x) const {
	return coeff * std::pow(x, 2) + bottom;
}

double Thacker1D::depth(double x, double t) const {
	double x1 = -0.5 * std::cos(omega * t) - a;
	double x2 = -0.5 * std::cos(omega * t) + a;
	if (x1 < x && x < x2) {
		return -coeff * (x + 0.5*cos(omega * t) - a)*(x + 0.5*std::cos(omega * t) + a);
	}else{
		return 0.0;
	}
}

double Thacker1D::speed(double x, double t) const {
	double x1 = -0.5 * std::cos(omega * t) - a;
	double x2 = -0.5 * std::cos(omega * t) + a;
	if (x1 < x && x < x2) {
		return 0.5 * omega * std::sin(omega * t);
	} else {
		return 0.0;
	}
}

}
