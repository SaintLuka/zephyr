#pragma once
#include "solution.hpp"

/// @class Step
/// @brief Computes dam break with a step solutions
/// @details
/// Class that computes the solutions for a dam break with a step without friction, see \cite Han14.
class Step : public Solution {

public:

	/// @brief Constructor
	explicit Step(Parameters&);

	/// @brief Destructor
	virtual ~Step();

	/// @brief Computes the solution
	void compute() override;

	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double, double, double) const;

	/// @brief Computes the value of the rarefaction function from the left state*/
	double r1(double, double, double);

	/// @brief Computes the speed of the shock wave in the second characteristic field*/
	double spshock(double, double, double);


private:
	double xdam, step_size; //the dam location and size of the step
	double h_left, h_right; // hl (resp. hr) the water heights on the left (resp. right) of the dam
	double u_left, u_right;
	double h_1, h_2; //water heights before, at and after the step respectively
	double u_1, u_2; //speed before, at and after the step respectively

	int solu;
};
