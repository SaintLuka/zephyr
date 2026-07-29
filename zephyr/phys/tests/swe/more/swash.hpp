#pragma once
#include "solution.hpp"

/// @class Swash
// * @brief Computes the solutions of the swash over an inclined plane
// * @details
// * Class that computes the solutions of the swash over an inclined plane, see \cite Marche05, \cite CaGr58.
// */
class Swash : public Solution {
public:
	/// @brief Constructor
	explicit Swash(Parameters &);

	/// @brief Destructor
	virtual ~Swash();

	/// @brief Computes the solution
	void compute() override;

	/// @brief Computes the non-dimensional speed ua and the non-dimensional free surface eta
	void ua_eta(double, double,double,double,double,double*, int);

	/// @brief Writes the parameters of the solution
	void param(double, double, double, double) const;
	
	/// @brief Writes the left boundary condition (at each dt = 0.01s)
	void leftcondition(int, double, double) const;

	/// @brief Computes the kth Bessel function for k=0,1 or 2
	double J(int, double);

private:
	double e; // initial curvature of the wave
	double alpha; // topography coefficient
	double tf; // time value
	double dt; // time step for the backup of the left boundary condition
	int nt; // number of time iterations
	int sol; // 1 if transient, 2 if periodic
	int flag; // to know if h=0.
	double A; // amplitude of the periodic solution
	double a;
	double ta; // non-dimentional time
	double* xexa; //non-dimensional x
	double* zexa; //non-dimensional z
	double x0; //abscissa changement
	double* uexa; //non-dimensional speed
	double* hexa;//non-dimentional water height
	double eta0; // initial free surface
	double* etaa;//non-dimensional free surface
	std::string namefile;

	/// @brief Copy constructor
	Swash(const Swash &);

	/// @brief operator=
	Swash & operator=(const Swash &);
};
