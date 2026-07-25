#pragma once
#include "solution.hpp"

/// @class MacDonald_like_diffus
/// @brief Computes Mac Donald solutions with diffusion
/// @details
/// Class that computes Mac Donald solutions in 1d with diffusion, see \cite Delestre10.

class MacDonald_like_diffus : public Solution {
public:
	/// @brief Constructor
	explicit MacDonald_like_diffus(Parameters &);

	/// @brief Destructor
	~MacDonald_like_diffus() override = default;

	/// @brief Computes the solution
	void compute() override;

	/// @brief Evaluation of the slope variation
	double Delta_topo_diffus(double, double, double, double, double, double, double, double ) const;

	/// @brief Writes the parameters of the solution
	void param(double, double ) const;

private:
	double varz, kt, kl, muv, muh; //slope variation, friction coefficients, viscosity
	Table1D dhex; //the water height variation (delta h)
	Table1D ddhex; //for the diffusion term (derivative of second order)
	double h_l_bound, h_r_bound; // water height on the left or right bound of the domain (differs from the discrete values of h)

	/// @brief Copy constructor
	MacDonald_like_diffus(const MacDonald_like_diffus &);

	/// @brief operator=
	MacDonald_like_diffus & operator=(const MacDonald_like_diffus &);
};
