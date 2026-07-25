#include <zephyr/phys/tests/swe/inclined_plane.hpp>


#include <zephyr/utils/pyplot.h>
#include <zephyr/utils/numpy.h>

#include <zephyr/phys/tests/swe/dam_break.h>

#include "zephyr/phys/tests/swe/rain.hpp"

using namespace zephyr;
using namespace zephyr::utils;
using namespace zephyr::phys;

int main() {

	Parameters params;
	params.nx_ex = 1000;
	params.ny_ex = 1;

	params.choice = 3;
	params.choicedomain = 1;

	Rain test(params);
	test.compute();

	//DamBreak test(1.0, 0.0, 0.0);

	pyplot plt;

	auto xs = test.xs_; //np::linspace(-10.0, 10.0, 1000);

	auto bed = test.bed_; //np::zeros_like(xs);
	auto depth = test.depth_; //np::zeros_like(xs);
	auto level = np::zeros_like(xs);

	double t = 0.5;
	for (int i = 0; i < xs.size(); ++i) {
		//bed[i] = bed(xs[i], t);
		//depth[i] = test.depth(xs[i], t);
		level[i] = bed[i] + depth[i]; //test.level(xs[i], t);
	}

	plt.plot(xs, bed, {.color="black"});
	plt.plot(xs, level, {.color="blue"});
	//plt.plot(xs, depth, {.color="green"});
	plt.show();

	return 0;
}