/*
  Copyright (C) 2017-2018  Semyon Glazyrin

  Released under the GNU GPL, see the GPL.txt file
  in the source distribution for its full text.
*/

#include <vector>
#include <math.h>
#include <algorithm>

#include <zephyr/math/calc/interpolation.h>


inline double sign(double x) { return (x > 0.0) ? 1.0 : ( (x < 0.0) ? -1.0 : 0.0 ); };
/// \brief Find index of the last element that less or equal than a given value
/// \param[in] x sorted vector for searching
/// \param[in] x0 value
/// \return first index i that x0 in the range [x[i], x[i+1]), -1 if x0 < x[0], x.size()-1 if x.back() < x0
inline int findArrayPos(std::vector<double>& x, double x0) { return std::upper_bound(x.begin(), x.end(), x0) - 1 - x.begin(); };
/// \brief Find index i that a given value in the range [x[i], x[i+1]) or this range is closest
/// \param[in] x sorted vector for searching
/// \param[in] x0 value
/// \return first index i that x0 in the range [x[i], x[i+1]), 0 if x0 < x[0], x.size()-2 if x.back() < x0
inline int findArrayPosInside(std::vector<double>& x, double x0) {
    //algorithm returns iter. in the range [x.begin()+1, x.end()-1], x0 < *iter; decrement and return vector index
    return std::upper_bound(x.begin()+1, x.end()-1, x0) - 1 - x.begin();
};


Interpolation1D::Interpolation1D(std::vector<double>& x, std::vector<double>& f)
{
  Nx = x.size();
  this->x = x;
  this->f = f;
}

void Interpolation1D::changeFunction(std::vector<double>& x, std::vector<double>& f)
{
  Nx = x.size();
  this->x = x;
  this->f = f;
}

// ---------------------------------------------------------------
//
//                     InterpLinear1D
//
// ---------------------------------------------------------------

double InterpLinear1D::calc(double x0)
{
  int i = findArrayPosInside(this->x, x0);
  return f[i]+(f[i+1]-f[i])*(x0-x[i])/(x[i+1]-x[i]);
}

double InterpLinear1D::diffX(double x0)
{
  int i = findArrayPosInside(this->x, x0);
  return (f[i+1]-f[i])/(x[i+1]-x[i]);
}

double InterpLinear1D::inverseX(double f0)
{
  int i = findArrayPosInside(this->f, f0);
  return x[i]+(x[i+1]-x[i])*(f0-f[i])/(f[i+1]-f[i]);
}

// ---------------------------------------------------------------
//
//                     InterpLog1D
//
// ---------------------------------------------------------------

double InterpLog1D::calc(double x0)
{
  double log_x0 = log(x0);
  double log_f = InterpLinear1D::calc(log_x0);
  return exp(log_f);
}

double InterpLog1D::diffX(double x0)
{
  double log_x0 = log(x0);
  double log_f = InterpLinear1D::calc(log_x0);
  double dlogfdlogx = InterpLinear1D::diffX(log_x0);
  return dlogfdlogx*exp(log_f)/x0;
}

double InterpLog1D::inverseX(double f0)
{
  double log_f0 = log(f0);
  double log_x = InterpLinear1D::inverseX(log_f0);
  return exp(log_x);
}

// ---------------------------------------------------------------
//
//                     InterpLogx1D
//
// ---------------------------------------------------------------

double InterpLogx1D::calc(double x0)
{
  double log_x0 = log(x0);
  return InterpLinear1D::calc(log_x0);
}


// ---------------------------------------------------------------
//
//                     Interpolation2D
//
// ---------------------------------------------------------------

Interpolation2D::Interpolation2D(std::vector<double>& x, std::vector<double>& y, double **f)
{
  Nx = x.size();
  Ny = y.size();
  this->x = x;
  this->y = y;
  this->f.resize(Nx);
  for(unsigned int i=0; i<Nx; i++)
  {
    this->f[i].resize(Ny);
    for(unsigned int j=0; j<Ny; j++)
      this->f[i][j] = f[i][j];
  }
}

Interpolation2D::Interpolation2D(std::vector<double>& x, std::vector<double>& y, std::vector< std::vector<double> >& f)
{
  Nx = x.size();
  Ny = y.size();
  this->x = x;
  this->y = y;
  this->f = f;
}

// ---------------------------------------------------------------
//
//                     InterpBilinear2D
//
// ---------------------------------------------------------------

void lim_value(double& v, double v_min, double v_max){
  if(v < v_min)
    v = v_min;
  else if(v > v_max)
    v = v_max;
}

double InterpBilinear2D::calc(double x0, double y0)
{
  int i = findArrayPosInside(this->x, x0);
  int j = findArrayPosInside(this->y, y0);
  if(this->extrapolate == 0) // no extrapolation
  {
    lim_value(x0, x.front(), x.back());
    lim_value(y0, y.front(), y.back());
  }
  double fracx = (x0-x[i])/(x[i+1]-x[i]);
  double fracy = (y0-y[j])/(y[j+1]-y[j]);
  return f[i][j]+(f[i+1][j]-f[i][j])*fracx+(f[i][j+1]-f[i][j])*fracy+
    (f[i+1][j+1]-f[i][j+1]-f[i+1][j]+f[i][j])*fracx*fracy;
}

double InterpBilinear2D::diffX(double x0, double y0)
{
  int i = findArrayPosInside(this->x, x0);
  int j = findArrayPosInside(this->y, y0);
  if(this->extrapolate == 0) // no extrapolation
  {
    lim_value(x0, x.front(), x.back());
    lim_value(y0, y.front(), y.back());
  }
  return (f[i+1][j]-f[i][j])/(x[i+1]-x[i])+
    (f[i+1][j+1]-f[i][j+1]-f[i+1][j]+f[i][j])/(x[i+1]-x[i])*(y0-y[j])/(y[j+1]-y[j]);
}

double InterpBilinear2D::diffY(double x0, double y0)
{
  int i = findArrayPosInside(this->x, x0);
  int j = findArrayPosInside(this->y, y0);
  if(this->extrapolate == 0) // no extrapolation
  {
    lim_value(x0, x.front(), x.back());
    lim_value(y0, y.front(), y.back());
  }
  return (f[i][j+1]-f[i][j])/(y[j+1]-y[j])+
    (f[i+1][j+1]-f[i][j+1]-f[i+1][j]+f[i][j])/(y[j+1]-y[j])*(x0-x[i])/(x[i+1]-x[i]);
}

double InterpBilinear2D::inverseX(double f0, double y0)
{
  int j = findArrayPosInside(this->y, y0);
  if(this->extrapolate == 0)
    lim_value(y0, y.front(), y.back());

  double fracy = (y0-y[j])/(y[j+1]-y[j]);

//   std::vector<double>& fY(Nx);
//   for(unsigned int i=0; i<Nx; i++)
//     fY[i] = f[i][j]+(f[i][j+1]-f[i][j])*fracy;
//   int i = findArrayPosInside(fY, f0);

  //interpolate f(x_i, y0) from f(x_i, y_j) and f(x_i, y_{j+1}), then use it for comparison in binary search
  auto comp = [&j, &fracy](const double& val, const std::vector<double>& fi) { return val < fi[j]+(fi[j+1]-fi[j])*fracy; };
  //same as findArrayPosInside implementation with array values calculated on the fly
  int i = std::upper_bound(f.begin()+1, f.end()-1, f0, comp) - 1 - f.begin();

  double f1 = f[i][j]+(f[i][j+1]-f[i][j])*fracy;
  double f2 = f[i+1][j]+(f[i+1][j+1]-f[i+1][j])*fracy;
  if(this->extrapolate == 0)
    lim_value(f0, f1, f2);

  return x[i]+(x[i+1]-x[i])*(f0-f1)/(f2-f1);
}

double InterpBilinear2D::inverseY(double f0, double x0)
{
  int i = findArrayPosInside(this->x, x0);
  if(this->extrapolate == 0)
    lim_value(x0, x.front(), x.back());

  double fracx = (x0-x[i])/(x[i+1]-x[i]);

//   std::vector<double>& fX(Ny);
//   for(unsigned int j=0; j<Ny; j++)
//     fX[j] = f[i][j]+(f[i+1][j]-f[i][j])*fracx;
//   int j = findArrayPosInside(fX, f0);

  //more complicated implementation compared to inverseX(), since the second index is searched
  auto start = f[i].data();
  auto comp = [&i, &start, &fracx, this](const double& val, const double& fij) { return val < fij+(this->f[i+1][&fij-start]-fij)*fracx; };
  //same as findArrayPosInside implementation with array values calculated on the fly
  int j = std::upper_bound(f[i].begin()+1, f[i].end()-1, f0, comp) - 1 - f[i].begin();

  double f1 = f[i][j]+(f[i+1][j]-f[i][j])*fracx;
  double f2 = f[i][j+1]+(f[i+1][j+1]-f[i][j+1])*fracx;
  if(this->extrapolate == 0)
    lim_value(f0, f1, f2);

  return y[j]+(y[j+1]-y[j])*(f0-f1)/(f2-f1);
}

// ---------------------------------------------------------------
//
//                     InterpBilinearLog2D
//
// ---------------------------------------------------------------

InterpBilinearLog2D::InterpBilinearLog2D(std::vector<double>& x, std::vector<double>& y, double **f) : InterpBilinear2D(x, y, f)
{
  for(unsigned int i=0; i<x.size(); i++)
    this->x[i] = log(this->x[i]);
  for(unsigned int i=0; i<y.size(); i++)
    this->y[i] = log(this->y[i]);
  for(unsigned int i=0; i<Nx; i++)
    for(unsigned int j=0; j<Ny; j++)
      this->f[i][j] = log(this->f[i][j]);
}

InterpBilinearLog2D::InterpBilinearLog2D(std::vector<double>& x, std::vector<double>& y, std::vector< std::vector<double> >& f) : InterpBilinear2D(x, y, f)
{
  for(unsigned int i=0; i<x.size(); i++)
    this->x[i] = log(this->x[i]);
  for(unsigned int i=0; i<y.size(); i++)
    this->y[i] = log(this->y[i]);
  for(unsigned int i=0; i<Nx; i++)
    for(unsigned int j=0; j<Ny; j++)
      this->f[i][j] = log(this->f[i][j]);
}

double InterpBilinearLog2D::calc(double x0, double y0)
{
  double log_x0 = log(x0);
  double log_y0 = log(y0);
  double log_f = InterpBilinear2D::calc(log_x0, log_y0);
  return exp(log_f);
}

double InterpBilinearLog2D::diffX(double x0, double y0)
{
  double log_x0 = log(x0);
  double log_y0 = log(y0);
  double log_f = InterpBilinear2D::calc(log_x0, log_y0);
  double dlogfdlogx = InterpBilinear2D::diffX(log_x0, log_y0);
  return dlogfdlogx*exp(log_f)/x0;
}

double InterpBilinearLog2D::diffY(double x0, double y0)
{
  double log_x0 = log(x0);
  double log_y0 = log(y0);
  double log_f = InterpBilinear2D::calc(log_x0, log_y0);
  double dlogfdlogy = InterpBilinear2D::diffY(log_x0, log_y0);
  return dlogfdlogy*exp(log_f)/y0;
}

double InterpBilinearLog2D::inverseX(double f0, double y0)
{
  double log_f0 = log(f0);
  double log_y0 = log(y0);
  double log_x = InterpBilinear2D::inverseX(log_f0, log_y0);
  return exp(log_x);
}

double InterpBilinearLog2D::inverseY(double f0, double x0)
{
  double log_f0 = log(f0);
  double log_x0 = log(x0);
  double log_y = InterpBilinear2D::inverseY(log_f0, log_x0);
  return exp(log_y);
}

// ---------------------------------------------------------------
//
//                     InterpBilinearLinLog2D
//
// ---------------------------------------------------------------

double InterpBilinearLinLog2D::calc(double x0, double y0)
{
  int i = findArrayPosInside(this->x, x0);
  int j = findArrayPosInside(this->y, y0);
  if( (f[i][j]*f[i+1][j]>=0.0) && (f[i][j]*f[i][j+1]>=0.0) && (f[i][j]*f[i+1][j+1]>=0.0) )
  {
    // Same signes: log interpolation
    double log_x2_x1 = log(x[i+1]/x[i]);
    double log_x0_x1 = log(x0/x[i]);
    double log_y2_y1 = log(y[j+1]/y[j]);
    double log_y0_y1 = log(y0/y[j]);
    double log_f21_f11 = log(f[i+1][j]/f[i][j]);
    double log_f12_f11 = log(f[i][j+1]/f[i][j]);
    double log_f22_f12 = log(f[i+1][j+1]/f[i][j+1]);

    double log_f_f11 = log_f21_f11*log_x0_x1/log_x2_x1
      + log_f12_f11*log_y0_y1/log_y2_y1
      + (log_f22_f12 - log_f21_f11)*(log_x0_x1*log_y0_y1)/(log_x2_x1*log_y2_y1);
    return f[i][j]*exp( log_f_f11 );
  }
  else
  {
    // Different signes: linear interpolation
    return f[i][j]+(f[i+1][j]-f[i][j])*(x0-x[i])/(x[i+1]-x[i])+(f[i][j+1]-f[i][j])*(y0-y[j])/(y[j+1]-y[j])+
      (f[i+1][j+1]-f[i][j+1]-f[i+1][j]+f[i][j])*(x0-x[i])/(x[i+1]-x[i])*(y0-y[j])/(y[j+1]-y[j]);
  }
}

// ---------------------------------------------------------------
//
//                     Interpolation3D
//
// ---------------------------------------------------------------

Interpolation3D::Interpolation3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double ***f)
{
  Nx = x.size();
  Ny = y.size();
  Nz = z.size();
  this->x = x;
  this->y = y;
  this->z = z;

  this->f.resize(Nx);
  for(unsigned int i=0; i<Nx; i++)
  {
    this->f[i].resize(Ny);
    for(unsigned int j=0; j<Ny; j++)
    {
      this->f[i][j].resize(Nz);
      for(unsigned int k=0; k<Nz; k++)
	this->f[i][j][k] = f[i][j][k];
    }
  }
}

Interpolation3D::Interpolation3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector< std::vector< std::vector<double> > >& f)
{
  Nx = x.size();
  Ny = y.size();
  Nz = z.size();
  this->x = x;
  this->y = y;
  this->z = z;
  this->f = f;
}

// ---------------------------------------------------------------
//
//                     InterpTrilinearLog3D
//
// ---------------------------------------------------------------

InterpTrilinearLog3D::InterpTrilinearLog3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double ***f)
{
  Nx = x.size();
  Ny = y.size();
  Nz = z.size();
  this->x = x;
  for(unsigned int i=0; i<x.size(); i++)
    this->x[i] = log(this->x[i]);
  this->y = y;
  for(unsigned int i=0; i<y.size(); i++)
    this->y[i] = log(this->y[i]);
  this->z = z;
  for(unsigned int i=0; i<z.size(); i++)
    this->z[i] = log(this->z[i]);

  this->f.resize(Nx);
  for(unsigned int i=0; i<Nx; i++)
  {
    this->f[i].resize(Ny);
    for(unsigned int j=0; j<Ny; j++)
    {
      this->f[i][j].resize(Nz);
      for(unsigned int k=0; k<Nz; k++)
	this->f[i][j][k] = log(f[i][j][k]);
    }
  }
}

InterpTrilinearLog3D::InterpTrilinearLog3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector< std::vector< std::vector<double> > >& f)
{
  Nx = x.size();
  Ny = y.size();
  Nz = z.size();
  this->x = x;
  for(unsigned int i=0; i<x.size(); i++)
    this->x[i] = log(this->x[i]);
  this->y = y;
  for(unsigned int i=0; i<y.size(); i++)
    this->y[i] = log(this->y[i]);
  this->z = z;
  for(unsigned int i=0; i<z.size(); i++)
    this->z[i] = log(this->z[i]);
  this->f = f;
  for(unsigned int i=0; i<Nx; i++)
    for(unsigned int j=0; j<Ny; j++)
      for(unsigned int k=0; k<Nz; k++)
	this->f[i][j][k] = log(f[i][j][k]);
}

double InterpTrilinearLog3D::calc(double x0, double y0, double z0)
{
  double log_x0 = log(x0);
  double log_y0 = log(y0);
  double log_z0 = log(z0);
  double log_f = InterpTrilinear3D::calc(log_x0, log_y0, log_z0);
  return exp(log_f);
}

// ---------------------------------------------------------------
//
//                     InterpTrilinear3D
//
// ---------------------------------------------------------------

double InterpTrilinear3D::calc(double x0, double y0, double z0)
{
  int i = findArrayPosInside(this->x, x0);
  int j = findArrayPosInside(this->y, y0);
  int k = findArrayPosInside(this->z, z0);

  double xd = (x0-x[i])/(x[i+1]-x[i]);
  double yd = (y0-y[j])/(y[j+1]-y[j]);
  double zd = (z0-z[k])/(z[k+1]-z[k]);
/*
  double c00 = f[i][j][k]*(1.0-xd) + f[i+1][j][k]*xd;
  double c01 = f[i][j][k+1]*(1.0-xd) + f[i+1][j][k+1]*xd;
  double c10 = f[i][j+1][k]*(1.0-xd) + f[i+1][j+1][k]*xd;
  double c11 = f[i][j+1][k+1]*(1.0-xd) + f[i+1][j+1][k+1]*xd;

  double c0 = c00*(1.0-yd) + c10*yd;
  double c1 = c01*(1.0-yd) + c11*yd;
  */

  int t = i;
  double tl1 = f[t][j][k]+(f[t][j+1][k]-f[t][j][k])*yd+(f[t][j][k+1]-f[t][j][k])*zd+(f[t][j+1][k+1]-f[t][j][k+1]-f[t][j+1][k]+f[t][j][k])*yd*zd;
  t = i+1;
  double tl2 = f[t][j][k]+(f[t][j+1][k]-f[t][j][k])*yd+(f[t][j][k+1]-f[t][j][k])*zd+(f[t][j+1][k+1]-f[t][j][k+1]-f[t][j+1][k]+f[t][j][k])*yd*zd;

  return tl1*(1.0-xd)+tl2*xd;

//  return c0*(1.0-zd) + c1*zd;
}

// ---------------------------------------------------------------
//
//
//
// ---------------------------------------------------------------


