#pragma once
#include "solution.hpp"

/// @class Thacker2D
/// @brief Computes %Thacker1D solutions in 2D
/// @details
/// Class that computes the solutions for %Thacker1D paraboloid, see \cite Thacker81.
class Thacker2D : public Solution {
public:
	/// @brief Constructor
	explicit Thacker2D(Parameters &);

	/// @brief Destructor
	virtual ~Thacker2D();

	/// @brief Computes the solution
	void compute() override;
	
	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double) const;

private:
	double omega, eta, zmin;
	Table2D radius2, hex2D, zex2D, uex2D, vex2D;
	double a, h0, A, Aa, cAa, r0;
	int solu; // the number of the solution
};
