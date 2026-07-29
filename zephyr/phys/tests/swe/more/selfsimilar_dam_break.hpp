#pragma once
#include "solution.hpp"

/// @class Selfsimilar_dam_break
/// @brief Computes self-similar dam break solutions
/// @details
/// Class that computes the self-similar solutions for dam break with friction,
/// see Self-similar_solutions.pdf in the doc folder or in the bibliography of sourcesup.
class Selfsimilar_dam_break : public Solution {
public:
	/// @brief Constructor
	explicit Selfsimilar_dam_break(Parameters &);

	/// @brief Destructor
	virtual ~Selfsimilar_dam_break();

	/// @brief Computes the solution
	void compute() override;
	
	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double ) const;

private:
	int choice;
	double xdam; // half-lenght of the dam
	double shift; // shift to get a domain between 0 and L.
	double xL, xR; //the dam location: the fluid is initially between xL and xR
	double hinit; //the initial height of the fluid
	double k1; // friction coefficient = -3*nu
	double k; // coefficient used in the approximate (kinematic or diffusive wave) equation
	double C1; // constant for the computation of the solution
	double Tm15; // Tm15 = T^{-1/5}
	double alpha, beta; // topography: zb = alpha x + beta
};
