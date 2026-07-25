#pragma once
#include "solution.hpp"

/// @class MacDonaldB2
/// @brief Computes Mac Donald pseudo 2d solutions
/// @details
/// Class that computes Mac Donald pseudo 2d solutions with bottom B2, see \cite MacDonald96.

class MacDonaldB2 : public Solution {
public:
	/// @brief Constructor
	explicit MacDonaldB2(Parameters &);

	/// @brief Destructor
	virtual ~MacDonaldB2();

	/// @brief Computes the solution
	void compute() override;
	
	/// @brief Writes the parameters of the solution
	void param(double, double, double ) const;
	
	/// @brief Evaluation of the slope variation
	double Delta_topo(double, double, double, double, double, double, double, double, double ) const;

private:
	std::vector<double> hpex, b, bp;
	double Z;  // slope of the boundaries of the channel
	double Q;  // discharge
	double n;  // Manning friction coefficient
	double res;
	double expo_1, expo_2, h_r_bound;  // temporary values: exponents and boundary values of h
};
