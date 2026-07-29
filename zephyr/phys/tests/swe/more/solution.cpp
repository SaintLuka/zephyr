#include "solution.hpp"

Solution::Solution(Parameters & par):nx_(par.get_nxex()),ny_(par.get_nyex()){
	
	/**
	 * @details	 
	 * @param[in] par contains all the values from the parameters file
	 */	
	
	allocation();

	t_end_=0;
	lx_=0;
	ly_=0;
	dx_=0;
	dy_=0;

}


void Solution::save_final_critical(const Table1D& xex, const Table1D& hex, const Table1D& uex, const Table1D& zex) const{
	
	/**
	 * @details	 
	 * Saves x (the position), h (the water height), u (the flow velocity), z (the topography), q (the flow discharge), 
	 * z+h (the free surface), Fr (the Froude number) and z+hc (the critical surface). 
	 * @param[in] xex abscissae
	 * @param[in] hex water height
	 * @param[in] uex flow velocity
	 * @param[in] zex topography
	 */	
	
	cout << "#(i-0.5)*dx " << "\t" << setw(9) << " h[i] " << "\t"<< setw(9) << " u[i] " << "\t"<< setw(9) << " topo[i] " << "\t"<< setw(9) << " q[i] " << "\t"<< setw(9) << " topo[i]+h[i] " <<"\t" <<setw(9) << "Fr[i]=Froude" << "\t"<< setw(9) << " topo[i]+hc[i] " << "\t"<< endl;
	for (int i=1;i<nx_+1;i++){
		if (hex[i]<EPSILON_H){
			cout << setprecision(7)<< setw(9)<< xex[i] << "\t"<< setw(9) << "0.0" << "\t"<< setw(9) <<"0.0" << "\t" << setw(9) <<zex[i] << "\t" << setw(9) <<"0.0"<< "\t" << setw(9) <<zex[i]<< "\t"<< setw(9) << "NaN" << "\t"<< setw(9) << pow(uex[i]*uex[i]*hex[i]*hex[i]/GRAV,1./3.)+zex[i]<< "\t"<< endl;
		}
		else{
			cout << setprecision(7)<< setw(9)<< xex[i] << "\t"<< setw(9) << hex[i] << "\t"<< setw(9) <<uex[i] << "\t" << setw(9) <<zex[i] << "\t" << setw(9) <<uex[i]*hex[i]<< "\t" << setw(9) <<hex[i]+zex[i]<< "\t"<< setw(9) << abs(uex[i])/sqrt(GRAV*hex[i]) << "\t"<< setw(9) << pow(uex[i]*uex[i]*hex[i]*hex[i]/GRAV,1./3.)+zex[i]<< "\t"<< endl;
		}
	}// end of i loop
}


void Solution::save_final_critical_init(const Table1D& xex, const Table1D& hex, const Table1D& uex, const Table1D& zex, std::span
                                     <const double> z0) const{
	
	/**
	 * @details	 
	 * Saves x (the position), h (the water height), u (the flow velocity), z (the topography), q (the flow discharge), 
	 * z+h (the free surface), Fr (the Froude number), z+hc (the critical surface),
	 * z0 (the initial topography) and z0+h (the initial surface).
	 * @param[in] xex abscissae
	 * @param[in] hex water height
	 * @param[in] uex flow velocity 
	 * @param[in] zex topography
	 * @param[in] z0 initial topography
	 */	
	
	cout << "#(i-0.5)*dx " << "\t" << setw(9) << " h[i] " << "\t"<< setw(9) << " u[i] " << "\t"<< setw(9) << " topo[i] " << "\t"<< setw(9) << " q[i] " << "\t"<< setw(9) << " topo[i]+h[i] " <<"\t" <<setw(9) << "Fr[i]=Froude" << "\t"<< setw(9) << " topo[i]+hc[i] " << "\t"<< setw(9) << "z0[i]=InitialTopo" << "\t"<< setw(9) << " z0[i]+h[i] " << "\t"<< endl;
	for (int i=1;i<nx_+1;i++){
		if (hex[i]<EPSILON_H){
			cout << setprecision(7)<< setw(9)<< xex[i] << "\t"<< setw(9) << "0.0" << "\t"<< setw(9) <<"0.0" << "\t" << setw(9) <<zex[i] << "\t" << setw(9) <<"0.0"<< "\t" << setw(9) <<zex[i]<< "\t"<< setw(9) << "NaN" << "\t"<< setw(9) << pow(uex[i]*uex[i]*hex[i]*hex[i]/GRAV,1./3.)+zex[i]<< "\t"<< setw(9) << z0[i] << "\t"<< setw(9) << z0[i] + hex[i] << "\t"<< endl;
		}
		else{
			cout << setprecision(7)<< setw(9)<< xex[i] << "\t"<< setw(9) << hex[i] << "\t"<< setw(9) <<uex[i] << "\t" << setw(9) <<zex[i] << "\t" << setw(9) <<uex[i]*hex[i]<< "\t" << setw(9) <<hex[i]+zex[i]<< "\t"<< setw(9) << abs(uex[i])/sqrt(GRAV*hex[i]) << "\t"<< setw(9) << pow(uex[i]*uex[i]*hex[i]*hex[i]/GRAV,1./3.)+zex[i]<< "\t"<< setw(9) << z0[i] <<"\t"<< setw(9) << z0[i]+hex[i] << "\t"<< endl;
		}
	}// end of i loop
}

void Solution::save_final_mu(const Table1D& xex, const Table1D& hex, const Table1D& zex) const{
	
	/**
	 * @details	 
	 * Saves x (the position), h (the water height),
	 * z (the topography) and z+h (the free surface). 
	 * @param[in] xex abscissae
	 * @param[in] hex water height
	 * @param[in] zex topography
	 */	
	
	cout << "#(i-0.5)*dx " << "\t" << setw(9) << " h[i] " << "\t"<< setw(9) <<" topo[i] " << "\t"<< setw(9) << " topo[i]+h[i] " << "\t"<< endl;
	for (int i=1;i<nx_+1;i++){
		cout << setprecision(7)<< setw(9) << xex[i] << "\t"<< setw(9) << hex[i] << "\t"<< setw(9) << zex[i]<< "\t" << setw(9) <<hex[i]+zex[i]<< "\t"<< endl;
	}// end of i loop
}

void Solution::save_final_2D(const Table1D& xex, const Table1D& yex, const Table2D& hex2D, const Table2D& uex2D, const Table2D& vex2D, const Table2D
                           & zex2D) const{
	
	/**
	 * @details	 
	 * Saves x and y (the position), h (the water height), u and v (the flow velocities in x and y),
	 * z+h (the free surface), z (the topography), U (the norm of the velocity), Fr (the Froude number), 
	 * qx, qy and q (the flow discharge in x, y and its norm).
	 * @param[in] xex abscissae
	 * @param[in] yex ordinates
	 * @param[in] hex2D water height
	 * @param[in] uex2D flow velocity in x
	 * @param[in] vex2D flow velocity in y
	 * @param[in] zex2D topography
	 */		
	
	cout << "#(i-0.5)*dx " << "\t" << setw(9) << "(j-0.5)*dy " << "\t" << setw(9) << " h[i][j] " << "\t"<< setw(9) << " u[i][j] " << "\t"<< setw(9) << " v[i][j] " << "\t"<< setw(9) << " topo[i][j]+h[i][j] "<< "\t"<< setw(9) << " topo[i][j] "<< "\t"<< setw(9) <<" ||U||[i][j] "<< "\t"<< setw(9) <<" Fr[i][j] "<< "\t"<< setw(9) <<" qx[i][j] "<< "\t"<< setw(9) <<" qy[i][j] " << "\t"<< setw(9) << "q[i][j]"<< endl;
	for (int i=1;i<nx_+1;i++){
		for (int j = 1 ; j< ny_+1;j++){
			if (hex2D[i][j]<EPSILON_H){
				cout << setprecision(7)<< setw(9)<< xex[i] << "\t"<< setw(9) << yex[j] << "\t"<< setw(9) << "0.0" << "\t"<< setw(9) <<"0.0" << "\t"<< setw(9) <<"0.0" << "\t" << setw(9) <<zex2D[i][j] << "\t" << setw(9) <<zex2D[i][j]<< "\t"<< setw(9) <<"0.0"<< "\t" << setw(9) <<"NaN" << "\t"<< setw(9) <<"0.0" << "\t" << setw(9) <<"0.0" << "\t"<< setw(9) << "0.0" << endl;
			}
			else{
				cout << setprecision(7)<< setw(9)<< xex[i] << "\t"<< setw(9) << yex[j] << "\t"<< setw(9) << hex2D[i][j] << "\t"<< setw(9) <<uex2D[i][j] << "\t"<< setw(9) <<vex2D[i][j] << "\t" << setw(9) <<hex2D[i][j]+zex2D[i][j] << "\t" << setw(9) <<zex2D[i][j]<< "\t"<< setw(9) << sqrt(uex2D[i][j]*uex2D[i][j] + vex2D[i][j]*vex2D[i][j])<< "\t" << setw(9) <<sqrt(pow(uex2D[i][j],2) + pow(vex2D[i][j],2))/sqrt(GRAV*hex2D[i][j]) << "\t"<< setw(9) <<hex2D[i][j]*uex2D[i][j] << "\t" << setw(9) <<hex2D[i][j]*vex2D[i][j] << "\t"<< setw(9) << hex2D[i][j]*sqrt(uex2D[i][j]*uex2D[i][j] + vex2D[i][j]*vex2D[i][j]) << endl;
			}
		}
		cout << endl;
		
	}// end of i loop
}

void Solution::save_final_spherical(const Table1D& lambdaex, const Table1D& thetaex, const Table2D& hex2D, const Table2D& uex2D, const Table2D& vex2D, const
                                  Table2D& rhoex) const{
	/**
	 * @details
	 * Saves x and y (the position), h (the water height), u and v (the flow velocities in x and y),
	 * z+h (the free surface), z (the topography), U (the norm of the velocity), Fr (the Froude number),
	 * qx, qy and q (the flow discharge in x, y and its norm).
	 * @param[in] lambdaex longitudinal angle
	 * @param[in] thetaex latitudinal angle
	 * @param[in] hex2D water height
	 * @param[in] uex2D longitudinale flow velocity component
	 * @param[in] vex2D latitudinale flow velocity component
	 * @param[in] rhoex topography
	 */

	cout << "#lambda " << "\t" << setw(9) << "theta " << "\t" << setw(9) << " h[i][j] " << "\t" << setw(9) << " u[i][j] " << "\t" << setw(9) << " v[i][j] " << "\t" << setw(9) << " topo[i][j]+h[i][j] " << "\t" << setw(9) << " topo[i][j] " << "\t" << setw(9) << " ||U||[i][j] " << "\t" << setw(9) << " Fr[i][j] " << "\t" << setw(9) << " qx[i][j] " << "\t" << setw(9) << " qy[i][j] " << "\t" << setw(9) << "q[i][j]" << endl;
	for (int i = 0; i < nx_+1; i++) {
		for (int j = 0; j < ny_ ; j++) {
			if (hex2D[i][j] < EPSILON_H) {
				cout << setprecision(7)<< setw(9)<< lambdaex[i] << "\t" << setw(9) << thetaex[j] << "\t" << setw(9) << "0.0" << "\t" << setw(9) << "0.0" << "\t" << setw(9) << "0.0" << "\t" << setw(9) << rhoex[i][j] << "\t" << setw(9) << rhoex[i][j] << "\t" << setw(9) << "0.0" << "\t" << setw(9) << "NaN" << "\t" << setw(9) << "0.0" << "\t" << setw(9) << "0.0" << "\t" << setw(9) << "0.0" << endl;
			}
			else {
				cout << setprecision(7)<< setw(9)<< lambdaex[i] << "\t" << setw(9) << thetaex[j] << "\t" << setw(9) << hex2D[i][j] << "\t" << setw(9) << uex2D[i][j] << "\t" << setw(9) << vex2D[i][j] << "\t" << setw(9) << hex2D[i][j] + rhoex[i][j] << "\t" << setw(9) << rhoex[i][j] << "\t" << setw(9) << sqrt(uex2D[i][j] * uex2D[i][j] + vex2D[i][j] * vex2D[i][j]) << "\t" << setw(9) << sqrt(pow(uex2D[i][j], 2) + pow(vex2D[i][j], 2)) / sqrt(GRAV * hex2D[i][j]) << "\t" << setw(9) << hex2D[i][j] * uex2D[i][j] << "\t" << setw(9) << hex2D[i][j] * vex2D[i][j] << "\t" << setw(9) << hex2D[i][j] * sqrt(uex2D[i][j] * uex2D[i][j] + vex2D[i][j] * vex2D[i][j]) << endl;
			}
		}
		cout << endl;

	}// end of i loop
}

void Solution::save_final_concentrations(const Table1D& xex, const Table1D& phiex, const Table1D& psiex, const Table1D& phi0, const Table1D& psi0) const{
	/**
	 * @details
	 * Saves x (the position), phi (the dissolved solute concentration), psi (the adsorbed solute concentration),
	 * and phi0 (the initial dissolved concentration), psib (the initial adsorbed solute concentration).
	 * @param[in] xex abscissae
	 * @param[in] phiex dissolved solute concentration
	 * @param[in] psiex adsorbed solute concentration
	 * @param[in] phi0 initial dissolved solute concentration
	 * @param[in] psi0 initial adsorbed solute concentration
	 */

	cout << "# (i-0.5)*dx " << "\t" << setw(9) << " phi[i] " << "\t" << setw(9) << " psi[i] " << "\t" << setw(9) << " phi0[i]"<< "\t" << setw(9) << " psi0[i]" << "\t"<< endl;
	for (int i=1;i<nx_+1;i++){
		cout << setprecision(7)<< setw(9)<< xex[i] << "\t"<< setw(9) << phiex[i] << "\t"<< setw(9) << psiex[i]<< "\t" << setw(9) << phi0[i]<< "\t" << setw(9) << psi0[i]<< "\t"<< endl;
	}// end of i loop
}


void Solution::head(const Parameters & par, const string & solutiontype, const string & solutionchoice) {
	/**
	 * @details	 
	 * @param[in] par parameter, contains all the values from the parameters file
	 * @param[in] solutiontype name of the type of the solution
	 * @param[in] solutionchoice name of the solution
	 */	
	
	cout << "##############################################################################" << endl;
	cout << "# Generated by " << VERSION <<endl;
	cout << "##############################################################################" << endl;
	cout << "# Dimension: " << par.get_choicedim() << endl;
	cout << "# Type: "<< par.get_choicetype() <<" (="<< solutiontype << ")"<< endl;
	cout << "# Domain: " << par.get_choicedomain() << endl;
	cout << "# Choice: "<< par.get_choice() <<" (="<< solutionchoice<<")" << endl;
	cout << "##############################################################################" << endl;
}

void Solution::allocation() {
	xs_.resize(nx_ + 1);
	ys_.resize(ny_ + 1);
	depth_.resize(nx_ + 1);
	speed_.resize(nx_ + 1);
	flow_discharge_.resize(nx_ + 1);
	bed_.resize(nx_ + 1);
}
