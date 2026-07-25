/**
 * @file selfsimilar_dam_break.cpp
 * @author Serge Bodjona (2013)
 * @author Carine Lucas <carine.lucas@univ-orleans.fr> (2014-2015)
 * @version 1.03.00
 * @date 2015-10-28
 *
 * @brief Computes self-similar dam break solutions
 * @details 
 * Analytic solution: self-similar solution for dam break with friction, \n
 * see Self-similar_solutions.pdf in the doc folder or in the bibliography of sourcesup.
 *
 * @copyright License Cecill-V2 \n
 * <http://www.cecill.info/licences/Licence_CeCILL_V2-en.html>
 *
 * (c) CNRS - Universite d'Orleans - INRA (France)
 */
/*
 * This file is part of SWASHES software.
 * <https://sourcesup.renater.fr/projects/swashes/>
 *
 * SWASHES = Shallow-Water Analytic Solutions for Hydraulic and
 * Environmental Studies.
 * This software is a computer program whose purpose is to compute analytic
 * solutions for Shallow-Water equations.
 *
 * LICENSE
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * <http://www.cecill.info>.
 *
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading, using, modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean that it is complicated to manipulate, and that also
 * therefore means that it is reserved for developers and experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and, more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 ******************************************************************************/

#include "selfsimilar_dam_break.hpp"

Selfsimilar_dam_break::Selfsimilar_dam_break(Parameters & par):Solution(par){
	
	/** 
	 * @details 
	 * Defines the physical parameters, the final time and prints the header with the configuration.
	 * @param[in] par contains all the values from the parameters
	 * @par Modifies 
	 * Solution#dx_ex, Solution#L, Solution#T, Solution#xex, Solution#zex to have the self-similar dam break configuration.
	 */

	lx_ = 20.;
	dx_ = lx_/nx_;
	k1 = -3*0.1;
	
	if (1==par.get_choice()){
		// self-similar dam break on a flat bottom
		choice = 1;
		
		t_end_ = 30.;
		shift = lx_/2;
		Tm15 = pow(t_end_, -1./5.);
		hinit = 0.4;
		xdam = 2.5;
		xL = shift-xdam;
		xR = shift+xdam;
		
		for (int i=1 ; i<=nx_ ; i++){
			xs_[i] = (i-0.5)*dx_;
			bed_[i] = 0.0;
			/*// To save the initial condition as the standard error output (using 2> name_of_file)
			if (xex[i]<xL||xex[i]>xR){
				cerr << xex[i] <<  "\t"<< setw(9) << 0. <<  "\t"<< setw(9) <<0. << endl;
			}
			else{
				cerr << xex[i] <<  "\t"<< setw(9) << hinit <<  "\t"<< setw(9) <<0. << endl;
			}*/
		}
		k = GRAV/k1;
		C1 = 0.811774 * pow(-5*k/3, 2./5.);
		
		head(par, "Self-similar dam break", "on a flat bottom");
		param(lx_, xL, xR, hinit, k1, dx_, t_end_);
		cout << "##############################################################################" << endl;
		
	}else{ // 2==par.get_choice()
		// self-similar dam break on an inclined plane
		choice = 2;
		
		t_end_ = 100.;
		alpha = -0.1;
		beta = 3;
		hinit = 0.1;
		shift = 2;
		xL = shift;
		xR = 10+shift;
		
		for (int i=1 ; i<=nx_ ; i++){
			xs_[i] = (i-0.5)*dx_;
			bed_[i] = alpha *xs_[i] + beta;
			/*// To save the initial condition as the standard error output (using 2> name_of_file)
			if (xex[i]<xL||xex[i]>xR){
				cerr << xex[i] <<  "\t"<< setw(9) << 0. <<  "\t"<< setw(9) <<0. << endl;
			}
			else{
				cerr << xex[i] <<  "\t"<< setw(9) << hinit <<  "\t"<< setw(9) <<0. << endl;
			}*/
		}
		k = alpha *GRAV* 3/k1;
		C1 = 1; // value of A, the initial mass of the fluid
		
		head(par, "Self-similar dam break", "on an inclined plane");
		param(lx_, xL, xR, hinit, k1, dx_, t_end_);
		cout << "# Topography: zb = " << alpha << " x + " << beta << endl;
		cout << "##############################################################################" << endl;
	}

}


Selfsimilar_dam_break::~Selfsimilar_dam_break(){

}


void Selfsimilar_dam_break::compute(){
	
	/**
	 * @details
	 * Computes the chosen self-similar dam break solution,
	 * see Self-similar_solutions.pdf in the doc folder or in the bibliography of sourcesup.
	 * @par Modifies
	 * Solution#hex.
	 */

	
	if (1 == choice){
		// self-similar dam break on a flat bottom
		for (int i=1 ; i<=nx_ ; i++){
			if ((xs_[i]-shift)*Tm15<-sqrt(2*C1) || (xs_[i]-shift)*Tm15>sqrt(2*C1)){
				depth_[i] = 0.;
			}
			else{
				depth_[i] = Tm15*pow(-3/(5*k)*(C1-0.5*(xs_[i]-shift)*(xs_[i]-shift)*Tm15*Tm15),1./3.);
			}
		}
	}
	else{// 2 == choice
		// self-similar dam break on an inclined plane
		for (int i=1 ; i<=nx_ ; i++){
			if ((xs_[i]-shift)<0 || (xs_[i]-shift)>pow(9*C1*C1*k*t_end_/4, 1./3.)){
				depth_[i] = 0.;
			}
			else{
				depth_[i] = sqrt((xs_[i]-shift)/(k*t_end_));
			}
		}
	}
	
	save_final_mu(xs_, depth_, bed_);
	
}

void Selfsimilar_dam_break::param(double L, double xL, double xR, double hinit, double k1, double dx_ex, double T) const{
	
	/**
	 * @details 
	 * @param[in] L length of the domain
	 * @param[in] xL left bound for the water
	 * @param[in] xR right bound for the water
	 * @param[in] hinit initial height of the fluid
	 * @param[in] k1 inverse of the friction coefficient in SW equations
	 * @param[in] dx_ex space step 
	 * @param[in] T final time
	 */	
	
	cout << "# PARAMETERS OF THE SOLUTION"<< endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters"<<endl;
	cout << "# Space step: "<< dx_ex << " meters"<< endl;
	cout << "# Number of cells: " << nx_ << endl;
	cout << "# Initial position of the fluid: between " << xL << " meters and " << xR << " meters"<< endl;
	cout << "# Initial height of the fluid: " << hinit << " meters"<< endl;
	cout << "# Laminar friction: " << k1 << " q / h^2" << endl;
	cout << "# Time value: " << T << " seconds" << endl;
}	
