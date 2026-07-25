#pragma once
#include "solution.hpp"

/// @class Sluice_gate
// * @brief Computes dam break with a sluice gate solutions
// * @details
// * Class that computes the solutions for a dam break with a sluice gate without friction, see \cite Cozzolino15.
// */
class Sluice_gate : public Solution {
public:
	/// @brief Constructor
	explicit Sluice_gate(Parameters&);

	/// @brief Destructor
	virtual ~Sluice_gate();

	/// @brief Computes the solution
	void compute() override;

	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double) const;

	/// @brief Computes the value of the free flow function
	double ff(double);

	/// @brief Computes the value of the rarefaction function from the left state*/
	double r1(double, double, double );
	
	/// @brief Computes the speed of the shock wave in the first characteristic field*/
	double spshock1(double, double, double);

	/// @brief Computes the speed of the shock wave in the second characteristic field*/
	double spshock2(double, double, double);

	/// @brief Computes the value of the shock function
	double s2(double, double, double );

	/// @brief Finds the intersection between two locus of admissible states to, the locus are chosen with the variable choice
	double dichotomie(int);

private:

	double Cc; //contraction coefficient

	double xdam, gate_size; //the dam location and size of the sluice gate
	double h_left, h_right; // hl (resp. hr) the water heights on the left (resp. right) of the dam
	double h_1, h_c, h_2; //water heights before, at and after the sluice gate respectively
	double u_1, u_c, u_2; //speed before, at and after the sluice gate respectively

	int solu;
};
