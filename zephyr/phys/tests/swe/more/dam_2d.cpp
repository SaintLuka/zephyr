#include "dam_2d.hpp"

Dam_2D::Dam_2D(Parameters& par) :Solution(par) {

	/**
	 * @details
	 * Defines the physical parameters.
	 * @param[in] par contains all the values from the parameters
	 * @par Modifies
	 * Solution#dx_ex, Solution#L, Solution#l, Solution#xex, Solution#yex, Dam_2D#zex2D, Dam_2D#uex2D, Dam_2D#vex2D
	 * to have Dam_2D configuration.
	 */

	zex2D.resize(nx_ + 1); // i : 0->NX_EX
	uex2D.resize(nx_ + 1); // i : 0->NX_EX
	vex2D.resize(nx_ + 1); // i : 0->NX_EX
	hex2D.resize(nx_ + 1); // i : 0->NX_EX

	for (int i = 0; i <= nx_; i++) {

		zex2D[i].resize(ny_ + 1); // j : 0->NY_EX
		uex2D[i].resize(ny_ + 1); // j : 0->NY_EX
		vex2D[i].resize(ny_ + 1); // j : 0->NY_EX
		hex2D[i].resize(ny_ + 1); // j : 0->NY_EX


		for (int j = 0; j <= ny_; j++) {
			//The speed equals 0 for all the domains of this class
			uex2D[i][j] = 0;
			vex2D[i][j] = 0;
		}
	}
	if (par.get_choicedomain() == 1) {

		lx_ = 25.;
		ly_ = 10.;
		dx_ = lx_ / nx_;
		dy_ = ly_ / ny_;

		dam_d = 10.;
		dam_h = 0.5;
		dam_w = 1.;

		//parameter that garanties the width of the dam as a whole is 4*dam_w:
		float alpha = (exp(-pow(dam_w / 2., 2.)) - exp(-pow(dam_w * 2., 2.))) / dam_h; 

		//parameter of the topography that garanties the width of the flat surface of the dam is dam_w meter:
		float beta = exp(-pow(dam_w / 2, 2.)) / alpha - dam_h; 

		for (int j = 0; j <= ny_; j++) {
			ys_[j] = (j - 0.5) * dy_; //We look at values taken in the center of each cell of size dx*dy
		}

		for (int i = 0; i <= nx_; i++) {
			xs_[i] = (i - 0.5) * dx_;
			for (int j = 0; j <= ny_; j++) {
				zex2D[i][j] = MIN(dam_h, MAX(0., exp(-pow(xs_[i] - dam_d - pow(ys_[j] - 5., 2.) / 25., 2.)) / alpha - beta));
			}
		}

		if (par.get_choice() == 1) { //for now there's only one choice in this domain
			solu = 1;
			head(par, "DAM_2D", "Curved DAM");
			param(lx_, ly_, dam_d, dam_h, dam_w, dx_, dy_);
			cout << "##############################################################################" << endl;

		}
	}
	else { //choicedomain()==2

		lx_ = 10.;
		ly_ = 10.;
		dx_ = lx_ / nx_;
		dy_ = ly_ / ny_;

		dam_d = 2.5;
		dam_h = 1.;
		dam_w = 1.;

		float alpha_ring = (exp(-pow(dam_w / 2., 2.)) - exp(-pow(dam_w * 2., 2.))) / dam_h; //parameter that garanties the width of the dam is 4*dam_w
		float beta_ring = exp(-pow(dam_w / 2, 2.)) / alpha_ring - dam_h; //parameter of the topography that garanties the width of the flat surface of the dam is dam_w meter

		float alpha_cross = (exp(-pow(dam_w / 2., 2.)) - exp(-pow(dam_w * 2., 2.))) / (dam_h/2); //same parameter than the ring but the cross is 2 times lower
		float beta_cross = exp(-pow(dam_w / 2, 2.)) / alpha_cross - (dam_h/2);

		for (int j = 0; j <= ny_; j++) {
			ys_[j] = (j - 0.5) * dy_; //We look at values taken in the center of each cell of size dx*dy
		}

		for (int i = 0; i <= nx_; i++) {
			xs_[i] = (i - 0.5) * dx_;
			for (int j = 0; j <= ny_; j++) {
				double croix = cross(xs_[i], ys_[j], alpha_cross, beta_cross);
				zex2D[i][j] = MAX(ring(xs_[i], ys_[j], alpha_ring, beta_ring), croix );
			}
		}
		if (par.get_choice() == 1) { //for now there's only one choice in this domain
			solu = 2;
			head(par, "DAM_2D", "Symmetrical cross Dam");
			param(lx_, ly_, dam_d, dam_h, dam_w, dx_, dy_);
			cout << "##############################################################################" << endl;
		}
	}
}

Dam_2D::~Dam_2D() {
	for (int i = 0; i <= nx_; i++) {
		zex2D[i].clear();
		uex2D[i].clear();
		vex2D[i].clear();
		hex2D[i].clear();
	}
	zex2D.clear();
	uex2D.clear();
	vex2D.clear();
	hex2D.clear();
}


void Dam_2D::compute() {

	/**
	 * @details
	 * Computes the chosen Dam_2D solution.
	 * @par Modifies
	 * Dam_2D#hex2D.
	 */

	if (solu == 1) {
		
		for (int j = 0; j <= ny_; j++) {
			for (int i = 0; i <= nx_; i++) {
				if (xs_[i] - dam_d - pow(ys_[j] - ly_/2, 2.) / 25.<0.) {//
					hex2D[i][j] = dam_h- zex2D[i][j];
				}
				else {
					hex2D[i][j] = 0;
				}
			}
		}
	}
	else{ //solu == 2
		for (int j = 0; j <= ny_; j++) {
			for (int i = 0; i <= nx_; i++) {
				if (norm(xs_[i] - lx_ / 2, ys_[j] - ly_ / 2) < dam_d) {
					hex2D[i][j] = dam_h - zex2D[i][j];
				}
				else{
					hex2D[i][j] = 0;
				}
			}
		}
	}

	save_final_2D(xs_, ys_, hex2D, uex2D, vex2D, zex2D);

}


void Dam_2D::param(double L, double l, double dam_d, double dam_h, double dam_w, double dx_ex, double dy_ex) const {

	/**
	 * @details
	 * @param[in] L length of the domain in x
	 * @param[in] l length of the domain in y
	 * @param[in] dam_h height of the dam
	 * @param[in] dam_d distance of the center of the dam to the upstream border
	 * @param[in] dam_w width of the flat section at the top of the dam
	 * @param[in] dx_ex space step in x
	 * @param[in] dy_ex space step in y
	 */

	cout << "# PARAMETERS OF THE SOLUTION" << endl;
	cout << "# " << endl;
	cout << "# Length of the domain: " << L << " meters" << endl;
	cout << "# Width of the domain: " << l << " meters" << endl;
	cout << "# Space step in x: " << dx_ex << " meters" << endl;
	cout << "# Space step in y: " << dy_ex << " meters" << endl;
	cout << "# Number of cells in x: " << nx_ << endl;
	cout << "# Number of cells in y: " << ny_ << endl;
	cout << "# Height of the dam: " << dam_h << " meters" << endl;
	cout << "# Distance between the upstream limit and the center of the dam: " << dam_d << " meters" << endl;
	cout << "# width of the dam: " << dam_w << " meters" << endl;
	cout << "# Topography of the form: z(x,y) = min(dam_h, max(0, exp( ( g(x,y) ^2 )/alpha - beta));" << endl;

}

double Dam_2D::norm(double x, double y)
{
	/**
	 * @details
	 * computes the norm of the point (x,y)
	 * @param[in] x first coordinate of the point
	 * @param[in] y second coordinate of the point
	 */
	return sqrt(pow(x, 2) + pow(y, 2));
}

double Dam_2D::ring(double x, double y, double alpha, double beta) {
	/**
	 * @details
	 * computes the height of the center ring at the point (x,y) with a topography of the form 
	 * z(x,y) = min(dam_h, max(0, exp( ( g(x,y) ^2 )/alpha - beta))
	 * @param[in] x first coordinate of the point
	 * @param[in] y second coordinate of the point
	 * @param[in] alpha parameter of the dam shape
	 * @param[in] beta parameter of the dam shape
	 */
	return (MIN(dam_h, MAX(0., exp(-pow(norm(x-lx_/2,y-ly_/2) - dam_d, 2.)) / alpha - beta)));
}

double Dam_2D::cross(double x, double y, double alpha, double beta) {
	/**
	 * @details
	 * computes the height of the cross shaped dam at the point (x,y) with a topography of the form 
	 * z(x,y) = min(dam_h, max(0, exp( ( g(x,y) ^2 )/alpha - beta))
	 * @param[in] x first coordinate of the point
	 * @param[in] y second coordinate of the point
	 * @param[in] alpha parameter of the dam shape
	 * @param[in] beta parameter of the dam shape
	 */
	double x_center = x - lx_ / 2, y_center = y - ly_ / 2;
	if (norm(x_center,y_center) < dam_d) {
		return 0;
	}
	else {
		double proj_orth = 0.0;
		//we check whether we should project on x=y or x=-y then compute the norm of the projection
		if (abs(x_center+y_center)>abs(x_center-y_center)) { 
			proj_orth = norm((x_center - y_center) / 2,-(x_center-y_center)/2);
		}
		else {
			proj_orth = norm((x_center + y_center) / 2, (x_center + y_center) / 2);
		}
		
		return (MIN(dam_h/2, MAX(0., exp(-pow(proj_orth, 2.)) / alpha - beta)));
	}
}
