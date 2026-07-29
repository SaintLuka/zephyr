#include "sampson.hpp"

Sampson::Sampson(Parameters & par):Solution(par){
	
	/** 
	 * @details 
	 * Defines the physical parameters, the final time and prints the header with the configuration.
	 * @param[in] par contains all the values from the parameters
	 * @par Modifies 
	 * Solution#dx_ex, Solution#L, Solution#T, Solution#xex, Solution#zex to have %Sampson configuration. 
	 */
	
	lx_ = 10000.;
	dx_ = lx_/nx_;
	a = 3000.;
	h0 = 10.;

	x1=0;
	x2=0;
	
	B=5.;
	tau=0.001;
	p = pow(8.*GRAV*h0/(a*a),1./2.);
	s = pow(p*p-tau*tau,1./2.)/2.;
	
	t_end_ = 6000.; //in seconds
	
	for (int i=0 ; i<=nx_ ; i++){
		xs_[i] = (i-0.5)*dx_;
		bed_[i] = h0*pow((xs_[i]-lx_/2.)/a,2.);
	}
	
	head(par, "Oscillations", "Planar surface in a parabola with a linear friction (Sampson's solution)");
	param(lx_, h0, a, B, tau, dx_, t_end_);
	
}

Sampson::~Sampson(){
}


void Sampson::compute(){
	
	/**
	 * @details 
	 * Computes %Sampson solution, see \cite Sampson06 \cite Sampson08.
	 * @par Modifies 
	 * Solution#hex, Solution#uex.
	 */
	
	x1 = a*a*exp(-tau*t_end_/2.)/(2.*GRAV*h0)*(-B*s*cos(s*t_end_)-tau*B*sin(s*t_end_)/2.)-a +lx_/2;
	x2 = a*a*exp(-tau*t_end_/2.)/(2.*GRAV*h0)*(-B*s*cos(s*t_end_)-tau*B*sin(s*t_end_)/2.)+a+lx_/2;
	
	for (int i=0 ; i<=nx_ ; i++){
		if (xs_[i]>x1 && xs_[i]<x2){
			speed_[i] = B*exp(-tau*t_end_/2.)*sin(s*t_end_);
			depth_[i] = h0+pow(a*B,2.)*exp(-tau*t_end_)/(8.*pow(GRAV,2.)*h0) *(-s*tau*sin(2.*s*t_end_)+(tau/2.-s)*(tau/2. +s)*cos(2.*s*t_end_))- pow(B,2.)*exp(-tau*t_end_)/(4.*GRAV)-exp(-tau*t_end_/2.)*(B*s*cos(s*t_end_)+tau*B*sin(s*t_end_)/2.)*(xs_[i]-lx_/2)/GRAV-bed_[i];
		}else{
			speed_[i] = 0.;
			depth_[i] = 0.;
		}
	}
	
	save_final_critical(xs_, depth_, speed_, bed_);

}


void Sampson::param(double L, double h0, double a, double B, double tau, double dx_ex, double T) const{
	
	/**
	 * @details
	 * @param[in] L length of the domain
	 * @param[in] h0 value of the topography in the center of the domain
	 * @param[in] a parameter of the topography
	 * @param[in] B constant for the initial condition
	 * @param[in] tau friction coefficient
	 * @param[in] dx_ex space step
	 * @param[in] T final time
	 */	
	
	cout << "# PARAMETERS OF THE SOLUTION"<< endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters"<<endl;
	cout << "# Space step: "<< dx_ex << " meters"<< endl;
	cout << "# Number of cells: " << nx_ << endl;
	cout << "# Time value: " << T << " seconds" << endl;
	cout << "# " << endl;
	cout << "# Topography: z(x) = h0 (x-L/2)^2/a^2, with h0="<< h0<< " meters and a=" << a <<" meters"<<endl;
	cout << "# Constant B for the initial condition B=" << B<< " m/s"<< endl;
	cout << "# Friction coefficient tau=" << tau<< " s^-1" << endl;
	cout << "##############################################################################"<<endl;
}	


