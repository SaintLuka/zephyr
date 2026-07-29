#include "macdonaldb2.hpp"

MacDonaldB2::MacDonaldB2(Parameters & par):Solution(par){
	
	/** 
	 * @details 
	 * Defines the physical parameters and prints the header with the configuration.\n 
	 * The solution is saved at the steady state. 
	 * @param[in] par contains all the values from the parameters
	 * @par Modifies 
	 * Solution#dx_ex, Solution#L, Solution#xex, Solution#hex to have Mac Donald configuration. 
	 */	

	hpex.resize(nx_+1);
	b.resize(nx_+1);
	bp.resize(nx_+1);
	
	Q = 20. ;			//discharge
	Z=2.;				// slope of the boundaries of the channel
	lx_ = 400;			// length of the domain
	dx_ = lx_/nx_;	// space step
	n=0.03;				// Manning friction coefficient
	
	expo_1=4./3.;		// temporary values for the formula of res
	expo_2=10./3.;

	res=0;

	for (int i=1 ; i<=nx_ ; i++){
		xs_[i] = (i-0.5)*dx_;
		bed_[i] = 0.;
		b[i]=10.-5.*exp(-50.*pow((xs_[i]/400.)-(1./3.),2.))-5.*exp(-50.*pow((xs_[i]/400.)-(2./3.),2.));														// boundary
		bp[i]=(5./4.)*(((xs_[i]/400.)-(1./3.))*exp(-50.*pow((xs_[i]/400.)-(1./3.),2.))+((xs_[i]/400.)-(2./3.))*exp(-50.*pow((xs_[i]/400.)-(2./3.),2.)));	// derivative of the boundary
	}
	
	if (1==par.get_choice()){ // subcritical
		head(par, "MacDonald pseudo2D", "Trapezoidal long channel B2 with subcritical flow");
		param(lx_, dx_, n);
		
		for (int i=1 ; i<=nx_ ; i++){
			depth_[i]=0.9+0.3*exp(-40.*pow((xs_[i]/400.)-(1./3.),2.))+0.2*exp(-35.*pow((xs_[i]/400.)-(2./3),2.));
			hpex[i]=-0.06*((xs_[i]/400.)-(1./3.))*exp(-40.*pow((xs_[i]/400.)-(1./3),2.))-0.035*((xs_[i]/400.)-(2./3.))*exp(-35.*pow((xs_[i]/400.)-(2./3.),2.));
		}
		
		h_r_bound = 0.9+0.3*exp(-40.*pow((2./3),2.))+0.2*exp(-35.*pow((1./3),2.));
		
		cout << "# Initial conditions: h = max("<<h_r_bound << "- z(x), 0) m and q = 0 m^3/s" << endl; // note that z(L) = 0
		cout << "# Imposed discharge on the left boundary: " << Q << " m^3/s"<< endl;
		cout << "# Imposed water height on the right boundary: "<< h_r_bound << " m" << endl;
		cout << "############################################################################## " << endl;
	}
	else if (2==par.get_choice()){ // smooth transition and hydraulic jump
		head(par, "MacDonald pseudo2D", "Trapezoidal long channel B2 with smooth transition and hydraulic jump");
		param(lx_, dx_, n);
		
		for (int i=1 ; i<=nx_ ; i++){
			if(xs_[i]<=120.){
				depth_[i]=0.9+0.25*(expm1(-xs_[i]/40.))+0.25*exp(15.*((xs_[i]/400.)-(3./10.))); //For small magnitude values of x, expm1 may be more accurate than exp(x)-1.
				hpex[i]=0.009375*exp(15.*((xs_[i]/400.)-(3./10.)))-0.00625*exp(-xs_[i]/40.);
			}
			else{
				depth_[i]=1.5*exp(0.16*((xs_[i]/400.)-1.))-0.3*exp(2.*((xs_[i]/400.)-1.))+exp(-0.09*(xs_[i]-120.))*(-0.183691+1.519577*((xs_[i]-120.)/280.)-18.234429*pow((xs_[i]-120.)/280.,2.));
				hpex[i]=-0.0015*exp(2.*((xs_[i]/400.)-1.))+0.0006*exp(0.16*((xs_[i]/400.)-1.))+exp(-0.09*(xs_[i]-120.))*((1.519577/280.)-(18.234429/140.)*((xs_[i]-120.)/280.))-0.09*exp(-0.09*(xs_[i]-120.))*(-0.183691+1.519577*((xs_[i]-120.)/280.)-18.234429*pow((xs_[i]-120.)/280.,2.));
			}
		
		}
		
		h_r_bound = 1.5-0.3+exp(-0.09*(280.))*(-0.183691+1.519577-18.234429);
		
		cout << "# Initial conditions: h = max("<<h_r_bound << "- z(x), 0) m and q = 0 m^3/s" << endl; // note that z(L) = 0
		cout << "# Imposed discharge on the left boundary: " << Q << " m^3/s"<< endl;
		cout << "# Imposed water height on the right boundary: "<< h_r_bound << " m" << endl;
		cout << "############################################################################## " << endl;
	}
	
}


void MacDonaldB2::compute(){
	
	/**
	 * @details 
	 * Computes Mac Donald solutions with bottom B2, see \cite MacDonald96.
	 * @par Modifies 
	 * Solution#zex.
	 */
	
	// integration of the topography
	
	res = Delta_topo(depth_[nx_], hpex[nx_], b[nx_], bp[nx_], Q, n, Z, expo_1, expo_2);
	bed_[nx_] = - 0.5*dx_*res;
	
	for (int i=nx_;i>=1;i--){
		res = Delta_topo(depth_[i], hpex[i], b[i], bp[i], Q, n, Z, expo_1, expo_2);
		bed_[i-1] = bed_[i] - dx_*res;
	}
	
	save_final_mu(xs_, depth_, bed_);

}


void MacDonaldB2::param(double L, double dx_ex, double n) const{
	
	/**
	 * @details
	 * @param[in] L length of the domain
	 * @param[in] dx_ex space step
	 * @param[in] n friction coefficient 
	 */	
	
	cout << "# PARAMETERS OF THE SOLUTION"<< endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters"<<endl;
	cout << "# Space step in x: "<< dx_ex << " meters"<< endl;
	cout << "# Number of cells in x: " << nx_ << endl;
	cout << "# Topography: z(x) saved in the output"<<endl;
	cout << "# Solution at the steady state" << endl;
	cout << "# "<<endl;
	cout << "# Manning's friction coefficient: " << n << " m^-1/3 s" << endl;
	cout << "# " << endl;
	
}

double MacDonaldB2::Delta_topo(double h, double hp, double b, double bp, double Q, double n, double Z, double exp1, double exp2) const{
	
	/**
	 * @details
	 * @param[in] h water height
	 * @param[in] hp derivative of the water height 
	 * @param[in] b boundary function
	 * @param[in] bp derivative of the boundary function
	 * @param[in] Q discharge
	 * @param[in] n friction coefficient
	 * @param[in] Z slope
	 * @param[in] exp1 exponent, equal to 4/3
	 * @param[in] exp2 exponent, equal to 10/3
	 * @return Value of \f$\displaystyle hp\left(\frac{Q^2(b+2Zh)}{g(h(b+Zh))^3}-1\right)-Q^2 n^2 \frac{(b+2h\sqrt{1+Z^2})^{exp1}}{(h(b+Zh))^{exp2}}+ \frac{Q^2bp}{gh^2(b+Zh)^3}\f$.
	 */	
	
	return hp*(((pow(Q,2.)*(b+2.*Z*h))/(GRAV*pow(h*(b+Z*h),3.)))-1.)-pow(Q*n,2.)*(pow(b+2.*h*sqrt(1.+pow(Z,2.)),exp1)/pow(h*(b+Z*h),exp2))+(pow(Q,2.)*bp)/(GRAV*pow(h,2.)*pow(b+Z*h,3.));
}	

MacDonaldB2::~MacDonaldB2(){
	hpex.clear();
	b.clear();
	bp.clear();
}

