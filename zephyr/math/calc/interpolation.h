/*
  Copyright (C) 2017-2018  Semyon Glazyrin

  Released under the GNU GPL, see the GPL.txt file
  in the source distribution for its full text.
*/

#ifndef INTERPOLATIONS_H_FILE

#define INTERP_BILINEAR2D     0
#define INTERP_BILINEARLOG2D  1

#include <cmath>
#include <iostream>
#include <vector>

/// \addtogroup Interp Interpolations
/// \ingroup Tools
/// @{

class Interpolation
{
protected:
  int extrapolate; // 1 - do extrapolation, 0 - no extrapolation
public:
  Interpolation()
    {
      this->extrapolate = 1;
    }
  virtual ~Interpolation() {};
  void setExtrapolation(int extrapolate)
    {
      this->extrapolate = extrapolate;
    }
};

/// \brief Base class for interpolation of the function \f$ y=f(x) \f$ or the inverse function \f$ x=f^{-1}(y) \f$.
///
/// All interpolation use binary search, so argument values must be sorted in ascending order for the function interpolation, same as
/// function values for the inverse function interpolation.
class Interpolation1D : public Interpolation
{
protected:
  unsigned int Nx;
  std::vector<double> x;
  std::vector<double> f;
public:
  Interpolation1D() {}
  Interpolation1D(std::vector<double>& x, std::vector<double>& f);

  /// \brief Calculate an interpolated function value for a given argument.
  /// Initial argument values must be in ascending order.
  virtual double calc(double x)     = 0;
  /// Calculate an interpolated function derivative value for given argument
  virtual double diffX(double x)    = 0;
  /// \brief Calculate an interpolated argument value for a given function value.
  /// Initial function values must be in ascending order.
  virtual double inverseX(double f) = 0;

  virtual void changeFunction(std::vector<double>& x, std::vector<double>& f);
};

/// \brief Linear 1-dimensional interpolation.
class InterpLinear1D : public Interpolation1D
{
public:
  InterpLinear1D() {}
  InterpLinear1D(std::vector<double>& x, std::vector<double>& f) : Interpolation1D(x, f) {}
  double calc(double x)     override;
  double diffX(double x)    override;
  double inverseX(double f) override;
};

/// \brief Linear 1-dimensional interpolation in a logarithmic scale of both function and argument.
/// Appropriate for a power-law functions.
class InterpLog1D : public InterpLinear1D
{
public:
  InterpLog1D() {}
  InterpLog1D(std::vector<double>& x, std::vector<double>& f) : InterpLinear1D(x, f) {
    for(unsigned int i=0; i<x.size(); i++)
    {
      this->x[i] = log(this->x[i]);
      this->f[i] = log(this->f[i]);
    }
  }
  double calc(double x)     override;
  double diffX(double x)    override;
  double inverseX(double f) override;
};

class InterpLogx1D : public InterpLinear1D
{
public:
  InterpLogx1D() {}
  InterpLogx1D(std::vector<double>& x, std::vector<double>& f) : InterpLinear1D(x, f) {
    for(unsigned int i=0; i<x.size(); i++)
      this->x[i] = log(this->x[i]);
  }
  ~InterpLogx1D() {}
  virtual double calc(double x);
};


/// \brief Base class for interpolation of the function \f$ z=f(x,y) \f$ or the inverse functions \f$ x=f^{-1}(y,z) \f$ or \f$ y=f^{-1}(x,z) \f$.
///
/// All interpolation use binary search, so argument values must be sorted in ascending order for the function interpolation, same as
/// function values for the inverse function interpolation.
class Interpolation2D : public Interpolation
{
protected:
  unsigned int Nx, Ny;
  std::vector<double> x, y;
  std::vector< std::vector<double> > f;

  Interpolation2D() {};
  Interpolation2D(std::vector<double>& x, std::vector<double>& y, double **f);
  Interpolation2D(std::vector<double>& x, std::vector<double>& y, std::vector< std::vector<double> >& f);
public:
  /// \brief Calculate an interpolated function value for given arguments.
  /// Initial argument values must be in ascending order.
  virtual double calc(double x, double y)     = 0;
  /// Calculate an interpolated function derivative with respect to the first variable for given argument values.
  virtual double diffX(double x, double y)    = 0;
  /// Calculate an interpolated function derivative with respect to the second variable for given argument values.
  virtual double diffY(double x, double y)    = 0;
  /// \brief Calculate a first argument value for a given function value and a second argument value.
  /// Initial function values must be in ascending order.
  virtual double inverseX(double f, double y) = 0;
  /// \brief Calculate a second argument value for a given function value and a first argument value.
  /// Initial function values must be in ascending order.
  virtual double inverseY(double f, double x) = 0;
};

class Interpolation3D : public Interpolation
{
protected:
  unsigned int Nx, Ny, Nz;
  std::vector<double> x, y, z;
  std::vector< std::vector< std::vector<double> > > f;
public:
  Interpolation3D() {};
  Interpolation3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double ***f);
  Interpolation3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector< std::vector< std::vector<double> > >& f);
  ~Interpolation3D() {};
  virtual double calc(double x, double y, double z)      { std::cerr << "Default Interpolation3D::calc()\n"; exit(1); };
  virtual double diffX(double x, double y, double z)     { std::cerr << "Default Interpolation3D::diffX()\n"; exit(1); };
  virtual double diffY(double x, double y, double z)     { std::cerr << "Default Interpolation3D::diffY()\n"; exit(1); };
  virtual double diffZ(double x, double y, double z)     { std::cerr << "Default Interpolation3D::diffZ()\n"; exit(1); };
  virtual double inverseX(double f, double y, double z)  { std::cerr << "Default Interpolation3D::inverseX()\n"; exit(1); };
  virtual double inverseY(double f, double x, double z)  { std::cerr << "Default Interpolation3D::inverseY()\n"; exit(1); };
  virtual double inverseZ(double f, double x, double y)  { std::cerr << "Default Interpolation3D::inverseZ()\n"; exit(1); };
};


/// \brief Linear 2-dimensional interpolation.
class InterpBilinear2D : public Interpolation2D
{
public:
  InterpBilinear2D() {};
  InterpBilinear2D(std::vector<double>& x, std::vector<double>& y, double **f)
    : Interpolation2D(x, y, f) {}
  InterpBilinear2D(std::vector<double>& x, std::vector<double>& y, std::vector< std::vector<double> >& f)
    : Interpolation2D(x, y, f) {}
  double calc(double x, double y)     override;
  double diffX(double x, double y)    override;
  double diffY(double x, double y)    override;
  double inverseX(double f, double y) override;
  double inverseY(double f, double x) override;
};

/// \brief Linear 2-dimensional interpolation in a logarithmic scale of both function and arguments.
class InterpBilinearLog2D : public InterpBilinear2D
{
public:
  InterpBilinearLog2D(std::vector<double>& x, std::vector<double>& y, double **f);
  InterpBilinearLog2D(std::vector<double>& x, std::vector<double>& y, std::vector< std::vector<double> >& f);
  double calc(double x, double y)     override;
  double diffX(double x, double y)    override;
  double diffY(double x, double y)    override;
  double inverseX(double f, double y) override;
  double inverseY(double f, double x) override;
};


/// \brief 2-dimensional interpolation that choose between linear (for different signes) and log (same signes) scales
class InterpBilinearLinLog2D : public InterpBilinear2D
{
public:
  InterpBilinearLinLog2D(std::vector<double>& x, std::vector<double>& y, double **f) : InterpBilinear2D(x, y, f) {};
  InterpBilinearLinLog2D(std::vector<double>& x, std::vector<double>& y, std::vector< std::vector<double> >& f) : InterpBilinear2D(x, y, f) {};
  double calc(double x, double y);
  // double diffX(double x, double y);
  // double diffY(double x, double y);
  // double inverseX(double f, double y);
  // double inverseY(double f, double x);
};


class InterpTrilinear3D : public Interpolation3D
{
public:
  InterpTrilinear3D() {};
  InterpTrilinear3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double ***f)
    : Interpolation3D(x, y, z, f) {}
  InterpTrilinear3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector< std::vector< std::vector<double> > >& f)
    : Interpolation3D(x, y, z, f) {}
  virtual double calc(double x, double y, double z);
};

class InterpTrilinearLog3D : public InterpTrilinear3D
{
public:
  InterpTrilinearLog3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, double ***f);
  InterpTrilinearLog3D(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z, std::vector< std::vector< std::vector<double> > >& f);
  double calc(double x, double y, double z);
};

//
//     Interpolations on Field
//


/// @}

#endif
#define INTERPOLATIONS_H_FILE
