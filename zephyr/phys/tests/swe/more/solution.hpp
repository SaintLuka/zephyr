#pragma once
#include <span>
#include "parameters.hpp"

/// @class Solution
/// @brief Analytic solution
/// @details
/// Class that contains all the common declarations for the solutions.
class Solution {
public:
	double t_end_; ///< Final time

	const int nx_; ///< Number of cells in x
	const int ny_; ///< Number of cells in y

	double lx_; ///< Dimensions of the domain in x
	double ly_; ///< Dimensions of the domain in y
	double dx_; ///< Space step in x
	double dy_; ///< Space step in y

	Table1D xs_; ///< Array for the first coordinate
	Table1D ys_; ///< Array for the second coordinate

	Table1D depth_; ///< Array for the water height
	Table1D speed_; ///< Array for the flow velocity
	Table1D flow_discharge_; ///< Array for the flow discharge
	Table1D bed_; ///< Array for the topography

public:
	/// @brief Constructor
	explicit Solution(Parameters &);

	/// @brief Destructor
	virtual ~Solution() = default;

	/// @brief Allocations of the tables
	void allocation();

	/// @brief Function to be specified in case
	virtual void compute() =0;
		
	/// @brief Saves the analytic solution at the final time with the critical height
	void save_final_critical(const Table1D&, const Table1D&, const Table1D&, const Table1D&) const;
	
	/// @brief Saves the analytic solution at the final time with the critical height and the initial topography
	void save_final_critical_init(const Table1D&, const Table1D&, const Table1D&, const Table1D&, std::span<const double>) const;
	
	/// @brief Saves the analytic solution at the final time without u
	void save_final_mu(const Table1D&, const Table1D&, const Table1D&) const;
	
	/// @brief Saves the analytic solution at the final time in 2D
	void save_final_2D(const Table1D&, const Table1D&, const Table2D&, const Table2D&, const Table2D&, const Table2D&) const;

	/// @brief Saves the analytic solution at the final time in a spherical geometry
	void save_final_spherical(const Table1D&, const Table1D&, const Table2D&, const Table2D&, const Table2D&, const Table2D&) const;
	
	/// @brief Saves the analytic solution at the final time when written in concentrations
	void save_final_concentrations(const Table1D&, const Table1D&, const Table1D&, const Table1D&, const Table1D&) const;
		
	/// @brief Writes the version of the software and the choice of the solution
	static void head(const Parameters &, const std::string &, const std::string &);
};
