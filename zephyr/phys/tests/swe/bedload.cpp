#include "bedload.hpp"

Bedload::Bedload(Parameters & par) : Solution(par) {
	/// @details
	/// Defines the physical parameters, the final time and prints the header with the configuration.
	/// @param[in] par contains all the values from the parameters
	/// @warning Problem: allocation of z0 failed
	/// @par Modifies
	/// Solution#dx_ex, Solution#L, Solution#T, Solution#xex.
	/// @note If the vector z0 cannot be allocated, the code will exit with failure termination code.
	
	z0.resize(nx_+1);
	
	lx_ = 15.;
	dx_ = lx_/nx_;
	t_end_ = 7.;
	
	p = 1.5;
	alpha = 0.005;
	beta = 0.005;
	C = 1.0;
	q = 1.0;

	ue2 = 0;

	if (par.get_choice()==1){
		// bedload with Grass equation
		
		A = 0.005;
		ucr2 = 0.0;
		
		uexl = sqrt(pow((-alpha*0.5*dx_+beta)/A,1./p)+ucr2);
		hexl = q/uexl;
		z0l =- (pow(uexl,3.)+2.*GRAV*q)/(2.*GRAV*uexl)+ C;
		zexl =-alpha*t_end_+z0l;
		
		uexr = sqrt(pow((alpha*(lx_+0.5*dx_)+beta)/A,1./p));
		hexr = q/uexr;
		z0r =- (pow(uexr,3.)+2.*GRAV*q)/(2.*GRAV*uexr)+ C;
		zexr =-alpha*t_end_+z0r;
		
		head(par, "Bedload", "with Grass equation");
		param(lx_, dx_, t_end_, uexl, hexl, z0l, zexl, uexr, hexr, z0r, zexr, alpha, beta, A, q, C, p);
		
	}else{ // par.get_choice()==2
		// bedload with Meyer-Peter and Muler equation
		
		k=8;
		f=0.25;
		d=0.0005;
		s= 2600.0/1000.0; // rho_s/rho_water
		tcr=0.047;
		
		c1 = f/(8.*(s-1.)*GRAV*d);
		c2 = sqrt((s-1.)*GRAV*pow(d,3.));
		ucr2 = tcr/c1;
		A = k*pow(c1,p)*c2;
		
		uexl = sqrt(pow((-alpha*0.5*dx_+beta)/A,1./p)+ucr2);
		hexl = q/uexl;
		z0l =- (pow(uexl,3.)+2.*GRAV*q)/(2.*GRAV*uexl)+ C;
		zexl =-alpha*t_end_+z0l;
		
		uexr = sqrt(pow((alpha*(lx_+0.5*dx_)+beta)/A,1./p));
		hexr = q/uexr;
		z0r =- (pow(uexr,3.)+2.*GRAV*q)/(2.*GRAV*uexr)+ C;
		zexr =-alpha*t_end_+z0r;
		
		head(par, "Bedload", "with Meyer-Peter and Muler equation");
		param(lx_, dx_, t_end_, uexl, hexl, z0l, zexl, uexr, hexr, z0r, zexr, alpha, beta, A, q, C, p);
		cout << "# kappa="<<k<<", f="<<f<<", sedim. diam. d="<<d<<" m, sedim. density rs="<<s*1000.0<<" kg/m^3, Shields stress tau_cr="<<tcr << endl;
		
	}
	
	for (int i=0 ; i<=nx_ ; i++){
		xs_[i] = (i-0.5)*dx_;
		z0[i] = 0.0;
		bed_[i] = 0.0;
		speed_[i] = 0.0;
		depth_[i] = 0.0;
	}
}

void Bedload::compute() {
	/// @details
	/// Computes the chosen bedload solution, see \cite Berthon12.
	/// @par Modifies
	/// Solution#hex, Solution#uex, Solution#zex.

	for(int i=1; i<=nx_; i++) {
		ue2 = pow((alpha*xs_[i]+beta)/A,1./p);
		speed_[i] = sqrt(ue2+ucr2);
		depth_[i] = q/speed_[i];
		z0[i] =- (pow(speed_[i],3.)+2.*GRAV*q)/(2.*GRAV*speed_[i])+ C;
		bed_[i] =-alpha*t_end_+z0[i];
	}
	save_final_critical_init(xs_, depth_, speed_, bed_, z0);
}


void Bedload::param(double L , double dx_ex, double T, double uexl, double hexl, double z0l, double zexl, double uexr, double hexr, double z0r, double zexr, double alpha, double beta, double A,double q, double C, double p) const {
	/// @details
	/// @param[in] L length of the domain
	/// @param[in] dx_ex space step
	/// @param[in] T final time
	/// @param[in] uexl value of the velocity on the left boundary
	/// @param[in] hexl value of the water height on the left boundary
	/// @param[in] z0l value of the inital topography on the left boundary
	/// @param[in] zexl value of the final topography on the left boundary
	/// @param[in] uexr value of the velocity on the right boundary
	/// @param[in] hexr value of the water height on the right boundary
	/// @param[in] z0r value of the inital topography on the right boundary
	/// @param[in] zexr value of the final topography on the right boundary
	/// @param[in] alpha parameter for Exner equation
	/// @param[in] beta parameter for Exner equation
	/// @param[in] A parameter for Exner equation
	/// @param[in] q parameter for Exner equation
	/// @param[in] C parameter for Exner equation
	/// @param[in] p parameter for Exner equation
	
	cout << "# PARAMETERS OF THE SOLUTION"<< endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters"<<endl;
	cout << "# Space step: "<< dx_ex << " meters"<< endl;
	cout << "# Number of cells: " << nx_ << endl;
	cout << "# Time value: " << T << " seconds" << endl;
	cout << "# " << endl;
	cout << "# Initial conditions: h, u and z0 (see solution and documentation)" << endl;
	cout << "# Left boundary conditions x="<<-0.5*dx_ex <<" m:" << endl;
	cout << "#   u="<< uexl<< " m/s, h=" << hexl<< " m, z0=" << z0l << " m and z=" <<zexl<<" m" << endl;
	cout << "# Right boundary conditions x="<<L+0.5*dx_ex <<" m:" << endl;
	cout << "#   u="<< uexr<< " m/s, h=" << hexr<< " m, z0=" << z0r << " m and z=" <<zexr<<" m" << endl;
	cout << "# alpha="<<alpha<<" m/s, beta="<<beta<<" m^2/s, A="<<A<<" s^2/m, q="<<q<<" m^2/s, C="<<C<<" m, p="<<p<< endl; 
}
