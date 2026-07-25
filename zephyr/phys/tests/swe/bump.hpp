#pragma once
#include "solution.hpp"

/// @class Bump
/// @brief Computes bump solutions
/// @details
/// Class that computes the solutions with a bump for the topography, see \cite Delestre13 and \cite Goutal97.
class Bump : public Solution {
public:
	/// @brief Constructor
	explicit Bump(Parameters &);

	/// @brief Destructor
	virtual ~Bump();

	/// @brief Computes the solution
	void compute() override;
		
	/// @brief Coefficient p for Cardano method
	double p(double, double, double ) const;

	/// @brief Coefficient q for Cardano method
	double q(double, double, double, double ) const;

	/// @brief Determinant for Cardano method
	double determinant(double, double ) const;

	/// @brief Computation of the 3rd order polynomia roots
	double height(double, double, double, double, double ) const;

	/// @brief Defines a, b, c, d in order to solve \f$ ah^3+bh^2+ch+d \f$
	void abcd(double, double, double, double, double &, double &, double &, double &);

	/// @brief Steady state RH relation
	double RHJump(double, double, double ) const;
	
	/// @brief Writes the parameters of the solution
	void param(double, double) const;

private:
	
	double q_in, h_out, hmiddle; // discharge value, h at the outflow, the height at the top of the bump
	double a, b, c, d; // the coefficients of the 3rd order polynomia
	int solu; // the relative number of the bump solution
	double epsi, epsilon; //definition of epsilons
	double H_MAX; //the maximum water height
	double zmax; // max of the topography
};
