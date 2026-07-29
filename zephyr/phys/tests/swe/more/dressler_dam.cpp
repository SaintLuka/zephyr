#include "dressler_dam.hpp"


DresslerDam::DresslerDam(Parameters & par):Solution(par) {
	/**
	 * @details 
	 * Defines the physical parameters, the final time and prints the header with the configuration.
	 * @param[in] par contains all the values from the parameters
	 * @warning Problem: allocation of hexd failed.
	 * @par Modifies 
	 * Solution#dx_ex, Solution#L, Solution#xex, Solution#zex to have Dressler dam break configuration. 
	 * @note If the vector hexd cannot be allocated, the code will exit with failure termination code.
	 */
	
	lx_ = 2000.;
	dx_ = lx_/nx_;
	t_end_ = 40.; //ttot
	h0 = 6.; //water height behind the dam
	xdam = lx_/2; //dam location
	C = 40.; //Chezy friction coefficient
	dt = 0.1; //time step for the tip location
	t = dt; //current time for the tip location

	uTip=0;
	xa=0;
	xb=0;
	alpha1=0;
	alpha2=0;
	Tg=0;
	c2=0;
	a=0;
	b=0;
	mEnd=0;
	miTip=0;
	mEndTip=0;
	
	hexd.resize(nx_ + 1);
	for (int i=0 ; i <= nx_ ; i++) {
		xs_[i] = (i-0.5)*dx_; //the domain definition
		bed_[i] = 0.; //flat topography
	}

	c0 = pow(GRAV*h0,0.5); //wave speed
	Cst = GRAV/(C*C); //constant linked to the friction
	uTip0 = 0.; //initialization of the velocity of the tip
	xf = xdam; //position of the wet/dry position
	
	head(par, "Dam break", "on a dry domain with friction (Dressler's solution)");
	param( lx_, xdam, C, dx_, t_end_);
}

void DresslerDam::compute(){
	
	/**
	 * @details 
	 * Computes Dressler solution, see \cite Dressler52.
	 * @par Modifies 
	 * Solution#hex, Solution#uex.
	 */

	for (int i=0 ; i<=nx_ ; i++){
		depth_[i] = 0.;
		speed_[i] = 0.;
	}


	while(t<t_end_){
		xa = xdam-c0*t; //location of the plateau h=h0 
		xb = xdam+2.*c0*t; //location of the plateau h=0 for the Ritter solution, the wet/dry transition is behind this location because of friction

		for (int i=0 ; i<=nx_ ; i++){
			if (xs_[i]<xa){
				speed_[i] = 0.;
				depth_[i] = h0; //the plateau
				hexd[i] = h0; //the plateau
			}else{
				if (xs_[i]>xb){
					speed_[i] = 0.;
					depth_[i] = 0.;
					hexd[i] = 0.;
				}else{
					alpha1 = 6./(5.*(2.-(xs_[i]-xdam)/(c0*t)))-2./3.+4.*pow(3.,0.5)/135.*pow(2.-(xs_[i]-xdam)/(c0*t),3./2.);
					alpha2 = 12./(2.-(xs_[i]-xdam)/(c0*t))-8./3.+8.*pow(3.,0.5)/189.*pow(2.-(xs_[i]-xdam)/(c0*t),3./2.)-108./(7.*pow(2.-(xs_[i]-xdam)/(c0*t),2.));

					speed_[i] = 2./3.*c0*(1.+(xs_[i]-xdam)/(c0*t))+Cst*GRAV*alpha2*t;
					depth_[i] = pow(c0/3.*(2.-(xs_[i]-xdam)/(c0*t))+Cst*GRAV*alpha1*t,2.)/GRAV;
					hexd[i] = depth_[i];

					mEnd = i;
				}//end if
			}//end if
		}//end for

		uTip = 0.; //intialization of the tip velocity

		/*
		 * Loop to get the tip velocity at time t
		 */
		for (int i=0 ; i<=nx_ ; i++){
			uTip = max(uTip,speed_[i]);//the velocity in the tip is maximum
		}

		xf = xf+dt*(uTip0+uTip)/2.; //location of the tip at time t by an iterative method (wet/dry transition)
		uTip0 = uTip; //velocity of the tip at time t

		t = t+dt; //incrementation of the time variable
	}//end while

	for (int i=0 ; i<=nx_ ; i++){
		if (abs(speed_[i]-uTip0)<EPSILON){
			miTip = i; //location (indice) of the begining of the tip thanks to the velocity: the velocity is constant in the tip
		}
	}

	for (int i=0 ; i<=nx_ ; i++){
		if (xf>xs_[i]){
			mEndTip = i; //location (indice) of the end (front) of the tip
		}
	}


	for (int i=miTip ; i<=mEndTip ; i++){
		speed_[i] = uTip0; //the velocity in the tip is constant
	}

	/*
	 * Second order interpolation in the tip to have an idea of the water height
	 * personnal communication with Valerio Caleffi
	 */

	Tg = (xs_[miTip]-xs_[miTip-1])/(depth_[miTip]-depth_[miTip-1]);

	c2 = xs_[mEndTip];
	a = (Tg*depth_[miTip]+c2-xs_[miTip])/pow(depth_[miTip],2.);
	b = Tg-2.*a*depth_[miTip];

	for (int i=miTip ; i<=mEndTip ; i++){
		depth_[i] = (-b-pow(pow(b,2.)-4.*a*(c2-xs_[i]),1./2.))/(2.*a);
	}


	for (int i=mEndTip+1 ; i<=mEnd ; i++){
		speed_[i] = 0.;
		depth_[i] = 0.;
		hexd[i] = 0.;
	}
	

	save_final_critical(xs_, depth_, speed_, bed_);
	
}

void DresslerDam::param(double L, double xdam, double C, double dx_ex, double T) const{
	
	/**
	 * @details
	 * @param[in] L length of the domain
	 * @param[in] xdam position of the dam
	 * @param[in] C Chezy friction coefficient
	 * @param[in] dx_ex space step
	 * @param[in] T final time
	 */	
	
	cout << "# PARAMETERS OF THE SOLUTION"<< endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters"<<endl;
	cout << "# Space step: "<< dx_ex << " meters"<< endl;
	cout << "# Number of cells: " << nx_ << endl;
	cout << "# Position of the dam: x=" << xdam << " meters" << endl;
	cout << "# Chezy friction coefficient: " << C << endl;
	cout << "# Time value: " << T << " seconds" << endl;
	cout << "##############################################################################"<<endl;
}	


