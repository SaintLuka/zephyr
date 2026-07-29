#pragma once
#include "solution.hpp"

/// @class Rain
/// @brief Computes a solution with mobile rain
/// @details
/// Class that computes the solutions with mobile rain, with different velocities compared to the flow, see \cite DeLu25.

class Rain : public Solution {
public:
	/// @brief Constructor
	explicit Rain(Parameters&);

	/// @brief Destructor
	~Rain() override = default;

	/// @brief Computes the solution
	void compute() override;

	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double, double, double, double, double) const;

	/// @brief Computes the value of the integral of the rain*/
	double rainint(double, double, double);
	
	/// @brief Computes the solution at time t*/
	void computet(double);

private:
	double S0; // opposite of the slope of the domain
	double R0; // maximal rain intensity
	double vr; // rain velocity
	double C; // friction coefficient
	double h0, q0; // water height and discharge
	double x0, Lr; // characteristics of the rain distribution
	int raincase; // choice of the case
	
	double t0 = 0;
	double tl = 0;
	double t0l = 0;
	double tr = 1;
};
