#pragma once

#include <vector>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <cstdlib>

#define MAX(a,b) (a>=b?a:b)
#define MIN(a,b) (a<=b?a:b)

#define GRAV_DEM 4.905
#define PI 3.14159265
#define EPSILON_H 1.e-12
#define EPSILON 1.e-12

#define VERSION "SWASHES version 1.05.01, 2026-04-09"

using std::pow;
using std::sqrt;
using std::exp;
using std::complex;

using namespace std;
inline constexpr double GRAV = 9.81;

using Table1D = std::vector<double>;
using Table2D = std::vector<std::vector<double>>;

