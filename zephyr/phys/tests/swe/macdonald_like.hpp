#pragma once
#include "solution.hpp"

/// @class MacDonald_like
/// @brief Computes Mac Donald solutions
/// @details
/// Class that computes Mac Donald solutions in 1d,
/// see \cite MacDonald96, \cite MacDonald97, \cite Delestre13 and \cite Vo08.

class MacDonald_like : public Solution {
public:
	/// @brief Constructor
	explicit MacDonald_like(Parameters &);

	/// @brief Destructor
	~MacDonald_like() final = default;

	/// @brief Computes the solution
	void compute() override;

	/// @brief Evaluation of the slope variation for Manning friction law
	double Delta_topo_Manning(double, double, double, double, double ) const;

	/// @brief Evaluation of the slope variation for Darcy-Weisbach friction law
	double Delta_topo_Darcy_Weisbach(double, double, double, double, double ) const;
	
	/// @brief Writes the parameters of the solution
	void param(double, double ) const;

private:
	double varz, cf, R; //slope variation, friction coefficient, rain intensity
	int choice_fric; //variable for the choice of the friction law: 0=Manning, 1=Darcy-Weisbach
	Table1D dhex; //the water height variation (delta h)
	double h_l_bound, h_r_bound; // water height on the left or right bound of the domain (differs from the discrete values of h)
};
