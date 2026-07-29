#include <zephyr/phys/tests/swe/thacker.h>

namespace zephyr::phys::swe {

// ------------------------------------------------------------------------------------------------
//                                         THACKER 1D
// ------------------------------------------------------------------------------------------------

Thacker1D::Thacker1D()
	: Thacker1D(1.0, -0.2) {
	lx = 1.0;
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

IBed::Ptr Thacker1D::topography() const {
	return geom::ParabolicBed::create(coeff, bottom);
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

// ------------------------------------------------------------------------------------------------
//                                         THACKER 2D
// ------------------------------------------------------------------------------------------------

Thacker2D::Thacker2D(int sol) : sol(sol) {
	bottom = -0.4;
	coeff  = 0.8;

	h0 = -bottom;
	a = std::sqrt(h0 / coeff);

	if (sol == 1) {
		eta = 0.2; // Не понял, что это
		omega = std::sqrt(2.0 * GRAV * h0) / a;
	}
	else {
		r0 = 0.4; // Пропорционально радиусу капли
		omega = std::sqrt(8.0 * GRAV * h0) / a;
	}
}

double Thacker2D::max_time() const {
	return 2.0 * M_PI / omega; // один период
}

double Thacker2D::bed(const Vector3d& v) const {
	return coeff * v.squaredNorm() + bottom;
}

IBed::Ptr Thacker2D::topography() const {
	return geom::ParabolicBed::create(coeff, bottom);
}

double Thacker2D::depth(const Vector3d& v, double t) const {
	double H;
	if (sol == 1) {
		H = ((eta*h0)/(a*a))*(2.0*v.x()*std::cos(omega * t) + 2.0*v.y()*std::sin(omega*t) - eta) - bed(v);
	}
	else {
		double A = (a*a - r0*r0)/(a*a + r0*r0);
		double B = std::sqrt(1.0 - A * A);
		double C = 1.0 - A * std::cos(omega * t);
		double xi = B / C;
		double r2 = v.squaredNorm();
		H = xi * (h0 - coeff * xi * r2);
	}
	return std::max(0.0, H);
}

Vector2d Thacker2D::speed(const Vector3d& v, double t) const {
	if (sol == 1) {
		return {
			-(eta * omega) * std::sin(omega * t),
			+(eta * omega) * std::cos(omega * t)
		};
	}
	else {
		double A = (a*a - r0*r0)/(a*a + r0*r0);
		double C = 1.0 - A * std::cos(omega * t);
		double mod = (0.5 * omega * A / C) * std::sin(omega * t);
		return mod * v.head<2>();
	}
}

} // zephyr::phys::swe
