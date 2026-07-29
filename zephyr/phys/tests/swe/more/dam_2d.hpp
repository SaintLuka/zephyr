#pragma once
#include "solution.hpp"

/// @class Dam_2D
/// @brief Computes Static dam solutions in 2D
/// @details
/// Class that computes the solutions with a dam in 2d, see \cite Delestre13.
class Dam_2D : public Solution {
public:
	/// @brief Constructor
	explicit Dam_2D(Parameters&);

	/// @brief Destructor
	virtual ~Dam_2D();

	/// @brief Computes the solution
	void compute() override;

	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double) const;

	/// @brief Computes the norm of a vector*/
	double norm(double, double);

	/// @brief Computes the topography of the center ring for the second domain
	double ring(double, double, double, double);

	/// @brief Computes the topography of the cross for the second domain
	double cross(double, double, double, double);

private:
	Table2D hex2D, zex2D, uex2D, vex2D;
	double dam_d, dam_h, dam_w;

	int solu; // the number of the solution
};

