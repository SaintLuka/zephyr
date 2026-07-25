#pragma once
#include "solution.hpp"

/// @class Sampson
/// @brief Computes %Sampson solution
/// @details
/// Class that computes the solution for %Sampson parabola with friction, see \cite Sampson06 \cite Sampson08.

class Sampson : public Solution {
public:
	/// @brief Constructor
	explicit Sampson(Parameters &);

	/// @brief Destructor
	virtual ~Sampson();

	/// @brief Computes the solution
	void compute() override;
	
	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double ) const;

private:
	double x1, x2;
	double a, h0, B, p, tau, s;
};
