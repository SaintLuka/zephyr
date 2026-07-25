#pragma once
#include "solution.hpp"

/// @class Bedload
/// @brief Computes solutions with bedload
/// @details
/// Class that computes the solutions where the bed is moving with bedload, see \cite Berthon12.
class Bedload : public Solution {
public:
	/// @brief Constructor
	explicit Bedload(Parameters &);

	/// @brief Destructor
	~Bedload() override = default;

	/// @brief Computes the solution
	void compute() override;
	
	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double ) const;
	
	/// @brief Writes a warning about the the solution
	void paramwarning() const;

private:
	Table1D z0; // initial topography
	double uexl, hexl, z0l, zexl, uexr, hexr, z0r, zexr; // boundary values
	double alpha, beta, A, C, q, ucr2, p, ue2;
	double k, f, d, s, tcr, c1, c2; // for MPM

	/// @brief Copy constructor
	Bedload(const Bedload &);

	/// @brief operator=
	Bedload & operator=(const Bedload &);
};
