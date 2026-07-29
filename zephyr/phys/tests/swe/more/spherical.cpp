#include "spherical.hpp"

Spherical::Spherical(Parameters& par) :Solution(par) {
	/// @details
	/// Defines the physical parameters.
	/// @param[in] par contains all the values from the parameters
	/// @warning Problem: allocation of lambdaex failed
	/// @warning Problem: allocation of thetaex failed
	/// @par Modifies
	/// Solution#dx_ex, Solution#dy_ex, Spherical#rhoex, Spherical#uex2D, Spherical#vex2D, Spherical#lambdaex,
	/// Spherical#thetaex, Spherical#alpha, Spherical#omega, Spherical#radius,
	/// to have a Spherical configuration and the domain paramters.


	//To correctly represent the values graphically the longitudinal angle (lambda) must loop from 0 to 2PI.
	//Since the points where lambda=2PI don't add new informations compared to the ones where lambda=0,
	//we decided to add those points without counting them in the number of points demanded by the user.
	//So overall there will be (NX_EX+1)*NY_EX calculations. 
	//But only NX_EX*(NY_EX-2)+2 points in the actual mesh (considering the poles).			

	lambdaex.resize(nx_ + 1); //lambdaex is the longitudinal angle
	thetaex.resize(ny_); //thetaex is the latitudinal angle

	rhoex.resize(nx_ + 1); // i : 0->NX_EX
	uex2D.resize(nx_ + 1); // i : 0->NX_EX
	vex2D.resize(nx_ + 1); // i : 0->NX_EX
	hex2D.resize(nx_ + 1); // i : 0->NX_EX

	for (int i = 0; i < nx_ + 1; i++) {
		rhoex[i].resize(ny_); // j : 0->NY_EX - 1
		vex2D[i].resize(ny_); // j : 0->NY_EX - 1
		uex2D[i].resize(ny_); // j : 0->NY_EX - 1
		hex2D[i].resize(ny_); // j : 0->NY_EX - 1

		for (int j = 0; j < ny_; j++) {
			//The topography equals 0 for all the domains of this class
			rhoex[i][j] = 0.;
		}
	}

	if (par.get_choicedomain() == 1) {
		alpha = 0.; //alpha is the angle between the spherical pole and the earth axis
		if (par.get_choice() == 1) { //for now there's only one choice in this domain
			solu = 1;
		}
	}
	else { //choicedomain = 2
		alpha = 0.406; //alpha is the angle between the spherical pole and the earth axis
		if (par.get_choice() == 1) { //for now there's only one choice in this domain
			solu = 2;
		}
	}
	dx_ = 2*PI / (nx_); //In spherical geometry, dx and dy correspond to incremants in the angle lambda and theta
	dy_ = PI / (ny_-1);

	h0 = 3000.;
	radius = 6.37122 * pow(10, 6);
	u0 = 2*PI*radius/(12*24*3600);
	omega = 7.292*pow(10,-5); //Omega correponds to the pulsation of earth rotation 

	lambdaex[0] = 0.;
	thetaex[0] = -PI/2;
	for (int i = 1; i < nx_ + 1; i++) {
		lambdaex[i] = lambdaex[i-1]+dx_; //lambdaex is the longitudinal angle
	}

	for (int j = 1; j < ny_; j++) {
		thetaex[j] = thetaex[j - 1] + dy_; //thetaex is the latitudinal angle
	}

	for (int i = 0; i < nx_ + 1; i++) {
		for (int j = 0; j < ny_; j++) {
			uex2D[i][j] = u0*(cos(thetaex[j]) * cos(alpha) + cos(lambdaex[i]) * sin(thetaex[j]) * sin(alpha));
			vex2D[i][j] = -u0*sin(lambdaex[i])*sin(alpha);
		}
	}

	head(par, "SPHERICAL", "global steady state");
	param(radius, alpha, omega, dx_, dy_);
	cout << "##############################################################################" << endl;
}

void Spherical::compute() {

	/**
/// @details
/// Computes the chosen Spherical solution.
/// @par Modifies
/// Spherical#hex2D.
	 */

	if (solu == 1) {

		for (int i = 0; i < nx_ + 1; i++) {
			for (int j = 0; j < ny_; j++) {
				hex2D[i][j] = h0 - (radius * omega * u0 + pow(u0, 2) / 2) * pow(-cos(lambdaex[i]) * cos(thetaex[j]) * sin(alpha) + cos(alpha) * sin(thetaex[j]), 2) / GRAV;			
			}
		}
	}
	else { //solu=2

		for (int i = 0; i < nx_ + 1; i++) {
			for (int j = 0; j < ny_; j++) {
				hex2D[i][j] = h0 - (radius * omega * u0 + pow(u0, 2) / 2) * pow(-cos(lambdaex[i]) * cos(thetaex[j]) * sin(alpha) + cos(alpha) * sin(thetaex[j]), 2) / GRAV;
			}
		}
	}
	save_final_spherical(lambdaex, thetaex, hex2D, uex2D, vex2D, rhoex);
}


void Spherical::param(double radius, double alpha, double omega, double dx_ex, double dy_ex) const {

	/**
/// @details
/// @param[in] radius radius of the sphere
/// @param[in] omega pulsation of the rotation of the sphere
/// @param[in] alpha angle between the sphere's rotation axis and the polar axis
/// @param[in] dx_ex angle step in lambda
/// @param[in] dy_ex angle step in theta
	 */

	cout << "# PARAMETERS OF THE SOLUTION" << endl;
	cout << "# " << endl;
	cout << "# Radius of the sphere: " << radius << " meters" << endl;
	cout << "# Pulsation of the rotation of the sphere: " << omega << " s^-1" << endl;
	cout << "# Angle between the sphere's rotation axis and the polar axis: " << alpha << " radiants" << endl;
	cout << "# Longitudinale angle step: " << dx_ex << " radiants" << endl;
	cout << "# Latitudinale angle step: " << dy_ex << " radiants" << endl;
	cout << "# Number of cells in lambda: " << nx_ + 1 << endl;
	cout << "# Number of cells in theta: " << ny_ << endl;
	cout << "# Topography: z(x) = 0" << endl;

}
