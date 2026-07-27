#pragma once
#include <zephyr/phys/tests/swe/sw_test_1D.h>

namespace zephyr::phys {

/// @details Computes the solutions for a dam break without friction
/// [1] A. Ritter. Die Fortpflanzung der Wasserwellen. Zeitschrift des Vereines
/// Deuscher Ingenieure, 36(33):947–954, 1892. 23, 24
/// [2] J. J. Stoker. Water Waves: The Mathematical Theory with Applications.
/// Pure and Applied Mathematics. Interscience Publishers, New York, USA, 1957. 23, 24.
class DamBreak : public SwTest1D {
public:
	/// @brief Constructor
	DamBreak(double depth_L, double depth_R, double x_dam = 0.0);

	/// @brief Destructor
	~DamBreak() override = default;

	double depth(double x, double t) const override;

	double speed(double x, double t) const override;

	IBed::Ptr topography() const override {
		return math::ConstBed::create(0.0);
	}

private:
	/// @brief Computes the solution
	void compute();

	double h_L, h_R; // hl (resp. hr) the water heights on the left resp. right) of the dam
	double x0; //the dam location

	double h_mid, u_mid, c_mid; //water height, velocity and wave velocity in the intermediate state (only for the dam break on wet soil: Stoker's solution)
	double v; //the shock velocity (only for the Stoker's solution)
	double c_L, c_R; // cl (resp. cr) left (resp. right) wave velocity
};

} // namespace zephyr::phys
