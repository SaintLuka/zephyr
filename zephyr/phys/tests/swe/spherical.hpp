#pragma once
#include "solution.hpp"

/// @brief Computes Static solutions in spherical geometry
/// @details Class that computes different static solutions
/// in spherical geometry, see \cite Williamson92.
class Spherical : public Solution {
public:
	/// @brief Constructor
	explicit Spherical(Parameters&);

	/// @brief Destructor
	~Spherical() override = default;

	/// @brief Computes the solution
	void compute() override;

	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double) const;

private:
	Table2D hex2D, vex2D, uex2D, rhoex;
	//here rhoex corresponds to the topography compared with the base sea level radius

	Table1D lambdaex; ///< lambdaex is the longitudinal angle
	Table1D thetaex;  ///< thetaex is the latitudinal angle

	double h0, u0, radius, omega, alpha;
	//Omega correponds to the pulsation of earth rotation 
	//alpha is the angle between the spherical pole and the earth axis

	int solu; // the number of the solution
};