#include "choice_solution.hpp"

int main(int argc, char ** argv) {
	
	/**
	 * @details
	 * Main function of SWASHES, see \cite Delestre13.
	 * @param[in] argc number of the arguments.
	 * @param[in] argv value of the arguments.
	 */
	
	Parameters par(argc,argv);
	Choice_solution * solution;
	solution = new Choice_solution(par);
	cout << "# " << endl;
	solution->compute();
	delete solution;
	exit (EXIT_SUCCESS);

}
