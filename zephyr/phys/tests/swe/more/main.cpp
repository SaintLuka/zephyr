#include "choice_solution.hpp"

void run_main(int argc, char ** argv) {
	Parameters par(argc,argv);
	Choice_solution* solution;
	solution = new Choice_solution(par);
	std::cout << "# " << "\n";
	solution->compute();
	delete solution;
}
