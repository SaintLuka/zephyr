#pragma once
#include "solution.hpp"

/// @class Inclined_plane
/// @brief Computes the solutions over an inclined plane
/// @details
/// Class that computes the solutions over an inclined plane, see \cite Delestre12.
class Inclined_plane : public Solution {
public:
	/// @brief Constructor
	explicit Inclined_plane(Parameters &);

	/// @brief Destructor
	~Inclined_plane() override = default;

	/// @brief Computes the solution
	void compute() override;
	
	/// @brief Coefficient p for Cardano method
	static double p(double, double, double ) ;
	
	/// @brief Coefficient q for Cardano method
	static double q(double, double, double, double );
	
	/// @brief Determinant for Cardano method
	static double determinant(double, double ) ;
	
	/// @brief Computation of the 3rd order polynomia roots
	double height(double, double, double, double, double ) const;
	
	/// @brief Defines a, b, c, d in order to solve \f$ ah^3+bh^2+ch+d \f$
	static void abcd(double, double, double, double, double &, double &, double &, double &);
	
	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double) const;

public:
	double q0, h0; // discharge value, h at the inflow
	double alpha, beta; // topography coefficients alpha x+beta
	double a, b, c, d; // the coefficients of the 3rd order polynomia
	double H_MAX;
};
