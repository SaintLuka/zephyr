#pragma once
#include "solution.hpp"

/// @class Solute
/// @brief Computes solute solutions
/// @details
/// Class that computes the solutions for a solute problem, see \cite BZ24.
class Solute : public Solution {
public:
	/// @brief Constructor
	explicit Solute(Parameters&);

	/// @brief Destructor
	~Solute() override = default;

	/// @brief Computes the solution
	void compute() override;

	/// @brief Writes the parameters of the solution
	void param(double, double, double, double, double, double, double, double, double, double, double, double) const;
	
	/// @brief Compute the inital Gaussian distribution of Solute#phi
	double phi0(double, double, double);
	
	/// @brief Compute the inital zero distribution of Solute#psi
	double psi0(double);

private:
	double mu, sigma; // parameters of the initial Gaussian dissolved solute concentration
	double lambda; // degradation constant (in s^-1)
	double u; // water velocity (in m/s)
	double C; // sediment mass concentration in suspension (in kg/m^3)
	double kd; // equilibrium distribution coefficient (in m^3/kg)
	double Km1; // desorption rate (in s^-1)

	std::string choice; // text that characterize the solution
	
	Table1D phiex;   // dissolved solute concentration (in kg/m^3)
	Table1D psiex;   // adsorbed solute concentration (in kg/m^3)
	Table1D tabphi0; // initial dissolved solute concentration at xex values (in kg/m^3)
	Table1D tabpsi0; // initial adsorbed solute concentration at xex values (in kg/m^3)
	
	double phix0=0.0; // left boundary on the dissolved solute concentration (in kg/m^3)
	double psix0=0.0; // left boundary on the adsorbed solute concentration (in kg/m^3)
};
