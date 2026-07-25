#pragma once
#include "solution.hpp"

/// @class Thacker
// * @brief Computes %Thacker solution
// * @details
// * Class that computes the solution for %Thacker parabola, see \cite Thacker81.
class Thacker : public Solution {
public:
	/// @brief Constructor
	explicit Thacker(Parameters &);

	/// @brief Destructor
	virtual ~Thacker();

	/// @brief Computes the solution
	void compute() override;
	
	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double) const;

private:
	double x1, x2;
	double omega;
	double a, h0, B;
};
