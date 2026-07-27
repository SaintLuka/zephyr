#pragma once
#include "solution.hpp"
#include "dam_break.h"
#include "selfsimilar_dam_break.hpp"
#include "dressler_dam.hpp"
#include "sluice_gate.hpp"
#include "rain.hpp"
#include "step.hpp"
#include "inclined_plane.hpp"
#include "bump.hpp"
#include "solute.hpp"
#include "dam_2d.hpp"
#include "spherical.hpp"
#include "macdonald_like.hpp"
#include "macdonald_like_diffus.hpp"
#include "thacker_1d.h"
#include "bedload.hpp"
#include "thacker2d.hpp"
#include "macdonaldb1.hpp"
#include "macdonaldb2.hpp"
#include "sampson.hpp"
#include "swash.hpp"

/// @class Choice_solution
/// @brief Choice of the solution
/// @details
/// Class that calls the chosen solution.
class Choice_solution {
private:
	Solution* sol;
	int dim2;

	/// @brief Copy constructor
	Choice_solution(const Choice_solution &);

	/// @brief operator=
	Choice_solution & operator=(const Choice_solution &);

public:
	/// @brief Constructor
	explicit Choice_solution(Parameters &);

	/// @brief Computes the solution
	void compute();

	/// @brief Destructor
	virtual ~Choice_solution();
};
