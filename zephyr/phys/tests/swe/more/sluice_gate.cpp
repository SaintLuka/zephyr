#include "sluice_gate.hpp"

Sluice_gate::Sluice_gate(Parameters& par) :Solution(par) {

	/**
	 * @details
	 * Defines the physical parameters, the final time and prints the header with the configuration.
	 * @param[in] par contains all the values from the parameters
	 * @par Modifies
	 * Solution#dx_ex, Solution#L, Solution#T, Solution#xex, Solution#zex, Sluice_gate#h_left, Sluice_gate#h_right, Sluice_gate#gate_size,
	 * Sluice_gate#solu, Sluice_gate#Cc, Sluice_gate#xdam
	 * to have the sluice gate opening configuration.
	 */
	Cc = 0.611;
	lx_ = 10.;
	xdam = lx_ / 2.;
	dx_ = lx_ / nx_;
	t_end_ = 6.e-8;


	//we also initialize a few variable that will be usefull later:
	h_1 = 0.;
	h_2 = 0.;
	u_1 = 0.;
	u_2 = 0.;
	h_c = 0.;
	u_c = 0.;

	for (int i = 0; i <= nx_; i++) {
		xs_[i] = (i - 0.5) * dx_;
		bed_[i] = 0.0;
	}

	if (par.get_choice() == 1) {
		// sluice gate opening on dry domain
		h_left = 0.005;
		h_right = 0.;
		gate_size = 0.2 * h_left; //gate_size must be kept <4/9*h_left for the solution to be correct.
		solu = 1;
		head(par, "Sluice Gate", "on a dry domain without friction");
		param(lx_, xdam, gate_size, h_left, h_right, dx_, t_end_);
	}
	else if (par.get_choice() == 2) {
		// sluice gate opening on wet domain
		h_left = 0.005;
		h_right = 0.002*h_left; //for this solution to hold h_right must be < Cc*a
		gate_size = 0.2 * h_left; //gate_size must be kept <4/9*h_left for the solution to be correct.
		solu = 2;
		head(par, "Sluice Gate", "on a slightly wet domain without friction");
		param(lx_, xdam, gate_size, h_left, h_right, dx_, t_end_);
	}
	else {//par.get_choice ==3
		// sluice gate opening on wet domain
		h_left = 0.005;
		h_right = 0.2 * h_left; //for this solution to hold h_right must be < Cc*a
		gate_size = 0.2 * h_left; //gate_size must be kept <4/9*h_left for the solution to be correct.
		solu = 3;
		head(par, "Sluice Gate", "on a dry domain without friction");
		param(lx_, xdam, gate_size, h_left, h_right, dx_, t_end_);
	}
}


Sluice_gate::~Sluice_gate() {

}


void Sluice_gate::compute() {

	/**
	 * @details
	 * Computes the chosen sluice gate solution.
	 * @par Modifies
	 * Solution#hex.
	 */

	// For the following solutions we will start each time by computing all the intermediate states: h_1,u_1,h_c,u_c as well as h_2,u_2 if necessary
	// Most intermediate states are found at the intersection of two locus of admissible states.
	// We will use a dichotomie function to find them.
	if (solu == 1) {
		
		//we find the first intermediate state at the intersection between the locus corresponding to the rarefaction wave 
		//and the locus of admissible discharge at the sluice gate.
		h_1 = dichotomie(1);
		//Since this state is in the rarefaction wave (in th efirst characteristic field) locus we use it to find u_1.
		u_1 = r1(h_1,h_left,0.0);
		//h_c is given by the definition of Cc
		h_c = Cc * gate_size;
		//There is a constant discharge at the sluice gate: h_1*u_1=h_c*u_c
		u_c = h_1 * u_1 / h_c;

		for (int i = 0; i <= nx_; i++) { //looks at all the different cases
			if (xs_[i] <= xdam) {// we check on which side of the dam we are

				//water height of the left side rarefaction wave:
				double h_wave = pow(2 * sqrt(GRAV * h_left) - (xs_[i] - xdam) / t_end_, 2.) / (9. * GRAV);
				if (h_wave >= h_left) {
					depth_[i] = h_left;
					speed_[i] = 0.;
				}
				else if(h_wave >= h_1) {
					depth_[i] = h_wave;
					speed_[i] = r1(h_wave,h_left, 0.0);
				}
				else {
					depth_[i] = h_1;
					speed_[i] = u_1;
				}
			}
			else { //xex[i]>xdam

				//water height of the right side rarefaction wave:
				double h_wave = pow(u_c + 2 * sqrt(GRAV * h_c) - (xs_[i] - xdam) / t_end_, 2.) / (9. * GRAV);

				if (xs_[i] >= xdam + t_end_ * (u_c + 2 * sqrt(GRAV * h_c))) {
					depth_[i] = h_right;
					speed_[i] = 0.;
				}
				else if (h_wave >= h_c) {
					depth_[i] = h_c;
					speed_[i] = u_c;
				}
				else {
					depth_[i] = h_wave;
					speed_[i] = r1(h_wave,h_c,u_c);
				}
			}

		}//end for
	}

	else if(solu==2) {
		//we find the first intermediate state at the intersection between the locus corresponding to the rarefaction wave 
		//and the locus of admissible discharge at the sluice gate.
		h_1 = dichotomie(1);
		//Since this state is in the rarefaction wave (in th efirst characteristic field) locus we use it to find u_1.
		u_1 = r1(h_1, h_left, 0.0);
		//h_c is given by the definition of Cc.
		h_c = Cc * gate_size;
		//There is a constant discharge at the sluice gate: h_1*u_1=h_c*u_c.
		u_c = h_1 * u_1 / h_c;
		//h_2 is found at the intersection between the locus of the rarefaction wave from the U_c and the shock with the state U_R.
		h_2 = dichotomie(2);
		//same as with u_1.
		u_2 = r1(h_2,h_c,u_c);

		for (int i = 0; i <= nx_; i++) { //looks at all the different cases
			if (xs_[i] <= xdam) {// we check on which side of the dam we are

				//water height of the left side rarefaction wave:
				double h_wave = pow(2 * sqrt(GRAV * h_left) - (xs_[i] - xdam) / t_end_, 2.) / (9. * GRAV);

				if (h_wave >= h_left) {
					depth_[i] = h_left;
					speed_[i] = 0.;
				}
				else if (h_wave >= h_1) {
					depth_[i] = h_wave;
					speed_[i] = r1(h_wave, h_left,0.0);
				}
				else {
					depth_[i] = h_1;
					speed_[i] = u_1;
				}
			}
			else { //xex[i]>xdam

				//water height of the right side rarefaction wave:
				double h_wave = pow(u_c + 2 * sqrt(GRAV * h_c) - (xs_[i] - xdam) / t_end_, 2.) / (9. * GRAV);

				if (xs_[i] >= xdam + t_end_ * spshock2(h_2,u_2,h_right)) {
					depth_[i] = h_right;
					speed_[i] = 0.;
				}
				else if (h_wave >= h_c) {
					depth_[i] = h_c;
					speed_[i] = u_c;
				}
				else if (h_wave >= h_2) {
					depth_[i] = h_wave;
					speed_[i] = r1(h_wave,h_c,u_c);
				}
				else {
					depth_[i] = h_2;
					speed_[i] = u_2;
				}
			}
		}//end for
	}

	else {//solu =3
		//we find the first intermediate state at the intersection between the locus corresponding to the rarefaction wave 
		//and the locus of admissible discharge at the sluice gate.
		h_1 = dichotomie(1);
		//Since this state is in the rarefaction wave (in th efirst characteristic field) locus we use it to find u_1.
		u_1 = r1(h_1, h_left, 0.0);
		//h_c is given by the definition of Cc
		h_c = Cc * gate_size;
		//There is a constant discharge at the sluice gate: h_1*u_1=h_c*u_c
		u_c = h_1 * u_1 / h_c;
		//h_2 is found at the intersection between the locus of the shock wave from the U_c and the shock with the state U_R.
		h_2 = dichotomie(3);
		//same principle as with u_1.
		u_2 = s2(h_c, u_c, h_2);
	
		for (int i = 0; i <= nx_; i++) { //looks at all the different cases
			if (xs_[i] <= xdam) {// we check on which side of the dam we are

				//water height of the left side rarefaction wave:
				double h_wave = pow(2 * sqrt(GRAV * h_left) - (xs_[i] - xdam) / t_end_, 2.) / (9. * GRAV);

				if (h_wave >= h_left) {
					depth_[i] = h_left;
					speed_[i] = 0.;
				}
				else if (h_wave >= h_1) {
					depth_[i] = h_wave;
					speed_[i] = r1(h_wave, h_left, 0.0);
				}
				else {
					depth_[i] = h_1;
					speed_[i] = u_1;
				}
			}
			else { //xex[i]>xdam

				if (xs_[i] >= xdam + t_end_ * spshock2(h_2, u_2, h_right)) {
					depth_[i] = h_right;
					speed_[i] = 0.;
				}
				else if (xs_[i] <= xdam + t_end_ * spshock1(h_c, u_c, h_2)) {
					depth_[i] = h_c;
					speed_[i] = u_c;
				}
				else {
					depth_[i] = h_2;
					speed_[i] = u_2;
				}
			}
		}//end for
	}
	save_final_critical(xs_, depth_, speed_, bed_);

}



void Sluice_gate::param(double L, double xdam, double gate_size, double h_left, double h_right, double dx_ex, double T) const {

	/**
	 * @details
	 * @param[in] L length of the domain
	 * @param[in] xdam position of the dam
	 * @param[in] gate_size size of the opening of the sluice gate
	 * @param[in] h_left height of the left side water
	 * @param[in] h_right height of the right side water
	 * @param[in] dx_ex space step
	 * @param[in] T final time
	 */

	cout << "# PARAMETERS OF THE SOLUTION" << endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters" << endl;
	cout << "# Space step: " << dx_ex << " meters" << endl;
	cout << "# Number of cells: " << nx_ << endl;
	cout << "# Position of the dam: x=" << xdam << " meters" << endl;
	cout << "# Size of the gate: gate_size=" << gate_size << " meters" << endl;
	cout << "# height of water on the left side=" << h_left << " meters" << endl;
	cout << "# height of water on the right side=" << h_right << " meters" << endl;
	cout << "# Time value: " << T << " seconds" << endl;
	cout << "##############################################################################" << endl;
}

double Sluice_gate::ff(double h){
	/**
	 * @details
	 * Computes the value of the free flow function
	 * @param[in] h water height, corresponds here to the where we compute the function
	 */
	h_c = Cc * gate_size;
	return h_c * sqrt(2 * GRAV * h) / (h * sqrt(1 + h_c / h));
}

double Sluice_gate::r1(double h_r, double h_l, double u_l){
	/**
	 * @details
	 * Computes the value of the speed of the rarefaction function in the first characteristic field
	 * @param[in] h_r water height of the right side state
	 * @param[in] h_l water height of the left side state, origin of the rarefacion wave
	 * @param[in] u_l water speed of the left side state, origin of the rarefacion wave
	 */
	return u_l + 2. * (sqrt(GRAV * h_l) - sqrt(GRAV * h_r));
}

double Sluice_gate::spshock1(double h_l, double u_l, double h_r) {
	/**
	 * @details
	 * Computes the speed of the shock wave in the first characteristic field
	 * @param[in] h_l water height of the left side state, origin of the rarefacion wave
	 * @param[in] u_l water speed of the left side state, origin of the rarefacion wave
	 * @param[in] h_r water height of the right side state
	 */
	return u_l - h_r * sqrt(GRAV * (h_l + h_r) / (2 * h_l * h_r));
}

double Sluice_gate::spshock2(double h_l, double u_l, double h_r) {
	/**
	 * @details
	 * Computes the speed of the shock wave in the second characteristic field
	 * @param[in] h_l water height of the left side state, origin of the rarefacion wave
	 * @param[in] u_l water speed of the left side state, origin of the rarefacion wave
	 * @param[in] h_r water height of the right side state
	 */
	return u_l + h_r * sqrt(GRAV * (h_l + h_r) / (2 * h_l * h_r));
}

double Sluice_gate::s2(double h_l, double u_l, double h_r){
	/**
	 * @details
	 * Computes the value of second characteristic field shock function
	 * @param[in] h_l water height of the left side state, origin of the rarefacion wave
	 * @param[in] u_l water speed of the left side state, origin of the rarefacion wave
	 * @param[in] h_r water height of the right side state
	 */
	return u_l - (h_r-h_l)*sqrt((GRAV/2)*(1/h_l+1/h_r));
}

double Sluice_gate::dichotomie(int choice){
	/**
	 * @details
	 * Finds the intersection between two locus of admissible states to, the locus are chosen with the variable choice
	 * @param[in] choice choice of which dichotomy to apply
	 */
	// Parameters for the dichotomy
	double eps = 0.0000001;
	int nmax = 1000;
	int iter = 0.;
	double h_a = eps;
	double h_b = h_left;
	double mid=0.0;
	double func=10,funca=10;

	switch (choice) {
	case 1: //finds the intersection between the free flow function and the rarefaction function

		while (fabs(func) > eps && iter < nmax) {
			iter = iter + 1;
			mid = (h_b + h_a) * 0.5;
			func = (ff(mid) - r1(mid, h_left, 0.0));
			funca = (ff(h_a) - r1(h_a, h_left, 0.0));
			if (func * funca >= 0.) {
				h_a = mid;
			}
			else {
				h_b = mid;
			}
		} //end while
		return mid;

	case 2: //finds the intersection between the rarefaction function and the shock function
		func = eps + 1;
		while (fabs(func) > eps && iter < nmax) {
			iter = iter + 1;
			mid = (h_b + h_a) * 0.5;
			func = (s2(mid, 0.0, h_right) - r1(mid, h_c, u_c));
			funca = (s2(h_a, 0.0, h_right) - r1(h_a, h_c, u_c));
			if (func * funca >= 0.) {
				h_a = mid;
			}
			else {
				h_b = mid;
			}
		} //end while
		return mid;
	case 3: //finds the intersection between two shock functions
		func = eps + 1;
		while (fabs(func) > eps && iter < nmax) {
			iter = iter + 1;
			mid = (h_b + h_a) * 0.5;
			func = (s2(mid, 0.0, h_right) - s2(h_c, u_c, mid));
			funca = (s2(h_a, 0.0, h_right) - s2(h_c, u_c, h_a));
			if (func * funca >= 0.) {
				h_a = mid;
			}
			else {
				h_b = mid;
			}
		} //end while
		return mid;
	default:
	break;
	}
	return 0.0;
}


