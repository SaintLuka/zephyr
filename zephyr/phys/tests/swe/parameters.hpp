#pragma once
#include "misc.hpp"

/// @class Parameters
/// @brief Gets parameters
/// @details
/// Class that reads the parameters, checks their values and
/// contains all the common declarations to get the values of the parameters.
class Parameters{
public:
	/// Number of cells in x.
	int nx_ex; 
	/// Number of cells in y.
	int ny_ex;
	/// Value corresponding to the dimension of the solution.
	double choicedim;
	/// Value corresponding to the type of the solution.
	int choicetype;
	/// Value corresponding to the chosen solution.
	int choice;
	/// Value corresponding to the domain of the solution.
	int choicedomain;
	
public:
	Parameters() = default;

	/// @brief Constructor
	Parameters(int, char**);

	/// @brief Destructor
	virtual ~Parameters();
	
	/// @brief Prints help
	void help() const;

	/// @brief Gives the number of cells in x
	int get_nxex() const;
	
	/// @brief Gives the number of cells in y
	int get_nyex() const;
	
	/// @brief Gives the dimension
	double get_choicedim() const;
	
	/// @brief Gives the type
	int get_choicetype() const;

	/// @brief Gives the chosen solution
	int get_choice() const;
	
	/// @brief Gives the domain
	int get_choicedomain() const;
};
