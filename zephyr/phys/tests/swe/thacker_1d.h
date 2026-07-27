#pragma once
#include "solution.hpp"
#include <zephyr/phys/tests/swe/sw_test_1D.h>

namespace zephyr::phys {

/// @class Thacker1D
// * @brief Computes %Thacker1D solution
// * @details
// * Class that computes the solution for %Thacker1D parabola, see \cite Thacker81.
class Thacker1D : public SwTest1D {
public:
	/// @brief Тест по умолчанию
	Thacker1D();

	/// @brief Constructor
	Thacker1D(double coeff, double bottom);


	/// @brief Destructor
	~Thacker1D() override = default;

	double x_min() const { return -lx; }

	double x_max() const { return +lx; }

	double max_time() const;

	double bed(double x) const override;

	double depth(double x, double t) const override;

	double speed(double x, double t) const override;

	IBed::Ptr topography() const override {
		return math::ParabolicBed::create(coeff, bottom);
	}

private:
	double coeff, bottom;
	double omega;
	double a, h0;
	double lx;
};

}
