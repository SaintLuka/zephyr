#include <zephyr/phys/tests/swe/inclined_plane.hpp>


#include <zephyr/utils/pyplot.h>
#include <zephyr/utils/numpy.h>

#include <zephyr/phys/tests/swe/dam_break.h>

#include "zephyr/phys/tests/swe/rain.hpp"
#include <zephyr/phys/tests/swe/sluice_gate.hpp>
#include <zephyr/phys/tests/swe/thacker_1d.h>

using namespace zephyr;
using namespace zephyr::utils;
using namespace zephyr::phys;

int main() {

	Parameters params;
	params.nx_ex = 1000;
	params.ny_ex = 1;

	params.choice = 3;

	//Sluice_gate test(params);

	Thacker1D test;

	pyplot plt;

	auto xs = np::linspace(-2.0, 2.0, 1000);

	auto bed = np::zeros_like(xs);
	auto depth = np::zeros_like(xs);
	auto level = np::zeros_like(xs);

	double t = 0.0;
	for (int i = 0; i < xs.size(); ++i) {
		bed[i] = test.bed(xs[i]);
		depth[i] = test.depth(xs[i], t);
		level[i] = bed[i] + depth[i]; //test.level(xs[i], t);
	}

	plt.plot(xs, bed, {.color="black"});
	plt.plot(xs, level, {.color="blue"});
	//plt.plot(xs, depth, {.color="green"});
	plt.show();

	return 0;
}