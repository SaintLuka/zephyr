#include "solute.hpp"

Solute::Solute(Parameters& par) :Solution(par) {
	/// @details
	/// Defines the physical parameters, the final time and prints the header with the solute configuration.
	/// @param[in] par contains all the values from the parameters
	/// @par Modifies
	/// Solution#dx_ex, Solution#L, Solution#T, Solution#xex,
	/// to have the solute configuration.

	lx_ = 1000.; // m
	dx_ = lx_ / nx_;
	t_end_ = 500.; // s
	
	lambda = 0.0; // modified below if degradation
	u = 1.0; // m/s
	kd = 20; // m^3/kg
	C = 0.1; // kg/m^3
	Km1 = 0.002; // s^-1
	
	phiex.resize(nx_ + 1);
	psiex.resize(nx_ + 1);
	tabphi0.resize(nx_ + 1);
	tabpsi0.resize(nx_ + 1);

	for (int i = 0; i <= nx_; i++) {
		xs_[i] = (i - 0.5) * dx_;
	}
		
	if (par.get_choice()==1 || par.get_choice()==3){ // initial concentrations
		mu= 70; // m
		sigma=20; // m
		for (int i = 0; i <= nx_; i++) {
			tabphi0[i] = phi0(xs_[i], mu, sigma);
			tabpsi0[i] = psi0(xs_[i]) ;
		}
		choice = "Initial concentration";
	}
	else {
		phix0 = 0.001; // kg/m^3
		psix0 = 0.0; // kg/m^3
		choice = "Boundary concentration";
	}
		
	if (par.get_choice() >2) { // degradation
		lambda = 0.003; //s^-1
		choice = choice + " with degradation";
	}		
		
	head(par, "Solute solution", choice);
	param(lx_, phix0, psix0, lambda, mu, sigma, u, kd, C, Km1, dx_, t_end_);
}

void Solute::compute() {

	/**
/// @details
/// Computes the chosen solute solution.
/// @par Modifies
/// Solute#phiex, Solute#psiex.
	 */


	for (int i = 0; i <= nx_; i++) { //looks at all the different cases
		if (t_end_*u>xs_[i]){ // H(T-x/u) = 1
			phiex[i] = phix0/(C*kd+1)*exp(-lambda*xs_[i]/u)+ 1/(C*kd+1)*phi0(xs_[i]-u*t_end_, mu, sigma)*exp(-lambda*t_end_)+ C*kd*phix0/(C*kd+1)*exp(-(Km1*C*kd+Km1+lambda)*xs_[i]/u)+ C*kd/(C*kd+1)*phi0(xs_[i]-u*t_end_, mu, sigma)*exp(-(Km1*C*kd+Km1+lambda)*t_end_);
			psiex[i] = C*kd*phix0/(C*kd+1)*exp(-lambda*xs_[i]/u)-C*kd*phix0/(C*kd+1)*exp(-(Km1*C*kd+Km1+lambda)*xs_[i]/u) +C*kd/(C*kd+1)*phi0(xs_[i]-t_end_*u, mu, sigma)*exp(-lambda*t_end_)- C*kd/(C*kd+1)*phi0(xs_[i]-t_end_*u, mu, sigma)*exp(-(Km1*C*kd+Km1+lambda)*t_end_);
		}
		else{ // H(T-x/u) = 0
			phiex[i] = 1/(C*kd+1)*phi0(xs_[i]-u*t_end_, mu, sigma)*exp(-lambda*t_end_)+ C*kd/(C*kd+1)*phi0(xs_[i]-u*t_end_, mu, sigma)*exp(-(Km1*C*kd+Km1+lambda)*t_end_);
			psiex[i] = C*kd/(C*kd+1)*phi0(xs_[i]-t_end_*u, mu, sigma)*exp(-lambda*t_end_)- C*kd/(C*kd+1)*phi0(xs_[i]-t_end_*u, mu, sigma)*exp(-(Km1*C*kd+Km1+lambda)*t_end_);
		}
	}//end for
	
	
	save_final_concentrations(xs_, phiex, psiex, tabphi0, tabpsi0);

}

double Solute::phi0(double x, double mu, double sigma)
{
	/**
/// @details
/// Computes the initial Gaussian disctribution of Solute#phi at the point x (in kg/m^3)
/// @param[in] x coordinate of the point
/// @param[in] mu postion of the maximum of the initial dissolved solute concentration
/// @param[in] sigma standard deviation of the initial dissolved solute concentration
	 */
	return exp(-pow((x-mu),2)/(2*pow(sigma,2)))*0.001;
}

double Solute::psi0(double x)
{
	/**
/// @details
/// Computes the inital zero distribution of Solute#psi at the point x  (in kg/m^3)
/// @param[in] x coordinate of the point (unused)
	 */
	(void) x;
	return 0.0;
}


void Solute::param(double L, double phix0, double psix0, double lambda, double mu, double sigma, double u, double kd, double C, double Km1, double dx_ex, double T) const {

	/**
/// @details
/// @param[in] L length of the domain
/// @param[in] phix0 left boundary condition (input) on the dissolved solute concentration (in kg/m^3)
/// @param[in] psix0 left boundary condition (input) on the adsorbed solute concentration (in kg/m^3)
/// @param[in] lambda degradation constant
/// @param[in] mu postion of the maximum of the initial dissolved solute concentration
/// @param[in] sigma standard deviation of the initial dissolved solute concentration
/// @param[in] u water velocity
/// @param[in] kd equilibrium distribution coefficient
/// @param[in] C sediment mass concentration in suspension
/// @param[in] Km1 desorption rate
/// @param[in] dx_ex space step
/// @param[in] T final time
	 */

	cout << "# PARAMETERS OF THE SOLUTION" << endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters" << endl;
	cout << "# Space step: " << dx_ex << " meters" << endl;
	cout << "# Number of cells: " << nx_ << endl;
	cout << "# Degradation constant: lambda=" << lambda << " second^(-1)"<<endl;
	cout << "# Water velocity =" << u << " meters/second" << endl;
	cout << "# Sediment mass concentration in suspension: C=" << C << " kilograms/meter^3"<< endl;
	cout << "# Equilibrium distribution coefficient: kd=" << kd << " meters^3/kilogram"<<endl;
	cout << "# Desorption rate: K-1=" << Km1 <<" second^(-1)"<< endl;
	if (choice.at(0) == 'B'){ // if boundary conditions
		cout << "# Left phi value: phix0=" << phix0 << " kilograms/meter^3"<< endl;
		cout << "# Left psi value: psix0=" << psix0 << " kilograms/meter^3"<< endl;
	}
	else{ // if initial condition
		cout << "# Initial dissolved solute concentration: "<<endl;
		cout << "#    phii(x) = 0.001 exp(-(x-" << mu << ")^2/(2*"<<sigma <<"^2)) kilograms/meter^3"<<endl;
	}
	cout << "# Time value: " << T << " seconds" << endl;
	cout << "##############################################################################" << endl;
}

