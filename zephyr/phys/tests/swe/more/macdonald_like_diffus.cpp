#include "macdonald_like_diffus.hpp"

MacDonald_like_diffus::MacDonald_like_diffus(Parameters & par) : Solution(par) {
	/** 
	 * @details 
	 * Defines the physical parameters and prints the header with the configuration.\n 
	 * The solution is saved at the steady state. 
	 * @param[in] par contains all the values from the parameters
	 * @warning Problem: allocation of dhex failed.
	 * @warning Problem: allocation of ddhex failed.
	 * @par Modifies 
	 * Solution#dx_ex, Solution#L, Solution#xex, Solution#hex, Solution#qex to have Mac Donald configuration. 
	 * @note If the vector dhex (or ddhex) cannot be allocated, the code will exit with failure termination code.
	 */	

	varz=0;

	dhex.resize(nx_ + 1);//array for the water height variations
	ddhex.resize(nx_ + 1);//array for the diffusion term (second order derivative)

	/***********************************************************************
	 * L=1000 m channels cases (linear and quadratic friction)
	 ***********************************************************************/
	lx_ = 1000.;
	dx_ = lx_/nx_; //space step
	
	for (int i=0 ; i<=nx_ ; i++){
		xs_[i] = (i-0.5)*dx_;
	}

	/***********************************************************************
	 * subcritical case
	 ***********************************************************************/
	if(par.get_choice()==1){
		
		for (int i=0 ; i<=nx_ ; i++){
			depth_[i] = (pow(4./GRAV,1./3.))*(1.+exp(-16.*pow(xs_[i]/1000.-1./2.,2.))/2.);
			dhex[i] = -2.*(pow(4./GRAV,1./3.))*(xs_[i]/1000.-1./2.)*exp(-16.*pow(xs_[i]/1000.-1./2.,2.))/125.;
			ddhex[i] = -pow(4./GRAV,1./3.)*(1.-32.*(xs_[i]/1000.-1./2.)*(xs_[i]/1000.-1./2.))*exp(-16.*(xs_[i]/1000.-1./2.)*(xs_[i]/1000.-1./2.))/62500.;
			flow_discharge_[i] = 1.5;
		}//end for
		
		kt = 0.01;
		kl = 0.001;
		muv = 0.01;
		muh = 0.001;
		
		h_r_bound = (pow(4./GRAV,1./3.))*(1.+exp(-16.*pow(lx_/1000.-1./2.,2.))/2.);
		
		head(par, "MacDonald", "long channel with subcritical flow and diffusion");
		param(lx_, dx_);
		cout << "# Values of the parameters: kt="<< kt<< " , kl="<< kl<< " , mu_v=" << muv<< " , mu_h=" << muh<< endl;
		cout << "# " << endl;
		cout << "# Initial conditions: h = 0 m and q = 0 m^2/s" << endl;
		cout << "# Imposed discharge on the left boundary: "<< flow_discharge_[1] <<" m^2/s"<< endl;
		cout << "# Imposed water height on the right boundary: "<< h_r_bound << " m" << endl;
		cout << "##############################################################################"<<endl;
		
		/***********************************************************************
		 * supercritical case
		 ***********************************************************************/
		}else{
			for (int i=0 ; i<=nx_ ; i++){
				depth_[i] = (pow(4./GRAV,1./3.))*(1.-exp(-36.*pow(xs_[i]/1000.-1./2.,2.))/5.);
				dhex[i] = (pow(4./GRAV,1./3.))*9.*exp(-36.*pow(xs_[i]/1000.-1./2.,2.))*((xs_[i]/1000.)-(1./2.))/625.;
				ddhex[i] = pow(4./GRAV,1./3.)*9.*(1.-72.*(xs_[i]/1000.-1./2.)*(xs_[i]/1000.-1./2.))*exp(-36.*(xs_[i]/1000.-1./2.)*(xs_[i]/1000.-1./2.))/625000.;
				flow_discharge_[i] = 2.5;
			}//end for
			
			kt = 0.005;
			kl = 0.001;
			muv = 0.01;
			muh = 0.1;
			
			h_l_bound = (pow(4./GRAV,1./3.))*(1.-exp(-36.*pow(-1./2.,2.))/5.);
			
			head(par, "MacDonald", "long channel with supercritical flow and diffusion");
			param(lx_, dx_);
			cout << "# Values of the parameters: kt="<< kt<< " , kl="<< kl<< " , mu_v=" << muv<< " , mu_h=" << muh<< endl; 
			cout << "# " << endl;
			cout << "# Initial conditions: h = 0 m and q = 0 m^2/s" << endl;
			cout << "# Imposed water height on the left boundary: "<< h_l_bound <<" m"<< endl;
			cout << "# Imposed discharge on the left boundary: "<< flow_discharge_[1] <<" m^2/s" << endl;
			cout << "##############################################################################"<<endl;
			
		} //end if
}

void MacDonald_like_diffus::compute(){
	
	/**
	 * @details 
	 * Computes Mac Donald solutions with diffusion, see \cite Delestre10.
	 * @par Modifies 
	 * Solution#zex.
	 */

	/***********************************************************************
	 * zex is the topography associated to the chosen steady flow
	 ***********************************************************************/
	
	varz = Delta_topo_diffus(flow_discharge_[nx_],depth_[nx_],dhex[nx_],ddhex[nx_],kt,kl,muv,muh);
	bed_[nx_] = 0.5*dx_*varz;
	// zex =0 on the boundary ("NX_EX + 1/2")
	
	for (int i=nx_ ; i>=1 ; i--){
			varz = Delta_topo_diffus(flow_discharge_[i],depth_[i],dhex[i],ddhex[i],kt,kl,muv,muh);
			bed_[i-1] = dx_*varz+bed_[i];
	}
	
	for (int i =0; i<=nx_; i++){
		if (abs(depth_[i])  > EPSILON) {
			speed_[i] = flow_discharge_[i] / depth_[i];
		}
		else {
			speed_[i] = 0.0;
		}
	}
	
	save_final_critical(xs_, depth_, speed_, bed_);
}


double MacDonald_like_diffus::Delta_topo_diffus(double q, double h, double dh, double ddh, double kt, double kl, double muv, double muh) const{
	
	/**
	 * @details
	 * @param[in] q discharge
	 * @param[in] h water height 
	 * @param[in] dh variation of the water height
	 * @param[in] ddh second order derivative of h
	 * @param[in] kt turbulent coefficient
	 * @param[in] kl laminar coefficient
	 * @param[in] muv vertical viscosity
	 * @param[in] muh horizontal viscosity
	 * @return Value of \f$\displaystyle \left(1-\frac{q^2}{gh^3}\right) dh + \frac{kl\,q}{gh^2(1+\frac{kl\, h}{3muv})}+ \frac{kt\,q^2}{gh^2(1+\frac{kl\, h}{3muv})^2} + 4muh \frac{q\, ddh - \frac{q\, dh^2}{h}}{gh^2}\f$.
	 */	
	
	return (1.-pow(q,2.)/(GRAV*pow(h,3.)))*dh+kl*q/(GRAV*pow(h,2.)*(1.+kl*h/(3.*muv)))+kt*pow(q,2.)/(GRAV*pow(h,2.)*pow(1.+kl*h/(3.*muv),2.))+4.*muh*(q*ddh-q*dh*dh/h)/(GRAV*h*h);
}

void MacDonald_like_diffus::param(double L, double dx_ex) const{
	
	/**
	 * @details
	 * @param[in] L length of the domain 
	 * @param[in] dx_ex space step 
	 */
	
	cout << "# PARAMETERS OF THE SOLUTION"<< endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters"<<endl;
	cout << "# Space step: "<< dx_ex << " meters"<< endl;
	cout << "# Number of cells: " << nx_ << endl;
	cout << "# Solution at the steady state" << endl;
	cout << "# "<< endl;
	
}

