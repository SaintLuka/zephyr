#pragma once
#include "solution.hpp"

/// @class DresslerDam
/// @brief Computes Dressler dam break solution
/// @details
/// Class that computes the solutions for a dam break with friction, see \cite Dressler52.
class DresslerDam : public Solution {
public:
	/// @brief Constructor
	explicit DresslerDam(Parameters &);

	/// @brief Destructor
	~DresslerDam() final = default;

	/// @brief Computes the solution
	void compute() override;
	
	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double ) const;

private:
	double h0; //water height behind the dam
	double xdam; //location of the dam
	double C; //Chezy friction coefficient C
	double dt, t; //time step and time for the algorithm
	double c0; //wave velocity linked to the water height h0
	double Cst; //constant g/C^2 linked to the Chezy friction coefficient
	double uTip0, uTip; //for the estimation of the velocity in the tip
	double xa, xb, xf; //several variables for coordinates
	double alpha1, alpha2;
	double Tg, c2, a, b; //variables for second order interpolation

	int mEnd, miTip, mEndTip; //variables for the location of the tip

	Table1D hexd; //array for the classical Dressler's solution
};
