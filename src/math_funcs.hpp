#ifndef MATH_HPP
#define MATH_HPP
#include <cmath>      // sqrt, fabs
#include <iostream>   // std::endl
#include <sstream>    // std::stringstream
#include <stdexcept>  // std::runtime_error
#include <vector>

// small & simple functions
template<typename T>
inline T _sqr(T x) { return x * x; }

template<typename T>
inline T _cub(T x) { return x * x * x; }

template<typename T>
inline T _max(T a, T b) { return a > b ? a : b; }

template<typename T>
inline T _min(T a, T b) { return a < b ? a : b; }

template<typename T>
inline T _min(T a, T b, T c) { return _min(_min(a, b), c); }

template<typename T>
inline int _sign(T x) { return x < 0. ? -1 : 1; }

template<typename T>
inline T _quadsol1(T a, T b, T c) {return (-b + sqrt(b*b - 4.*a*c))/(2.*a);}

template<typename T>
inline T _quadsol2(T a, T b, T c) {return (-b - sqrt(b*b - 4.*a*c))/(2.*a);}

// binary search for the index of the value in an ordered array
template<typename T>
int _locate(
        T const *xx,
        T x,
        int n,
        int d = 1
        )
{
    if (n <= 0) return -1;
    int ju, jm, jl, j;
    int monotonicity = 0;

    jl = 0;
    ju = n + 1;

    if (n == 1)
        monotonicity = 1;
    else if (xx[d] > xx[0])
        monotonicity = 1;

    xx -= d;

    while (ju - jl > 1) {
        jm = (ju + jl) >> 1;
        if ((x >= xx[jm * d] && monotonicity == 1) || (x <= xx[jm * d] && monotonicity == 0))
            jl = jm;
        else
            ju = jm;
    }

    if (x == xx[d]) {
        j = 1;
    } else {
        if (x == xx[n * d]) j = n - 1;
        else j = jl;
    }

    j -= 1;

    return j;
}

// multidimensional linear interpolation
// method = 0 linear-linear
// method = 1 log-log
template<typename T>
void _interpn(
        T *val,
        T const *coor,
        T const *data,
        T const *axis,
        int const *len,
        int ndim,
        int nval = 1,
        int method = 0
        )
{
    int i1, i2;

    i1 = _locate(axis, *coor, *len);

    // if the interpolation value is out of bound
    // use the closest value
    if (i1 == -1) {
        i1 = 0; i2 = 0;
    } else if (i1 == *len - 1) {
        i1 = *len - 1; i2 = *len - 1;
    } else i2 = i1 + 1;

    T *v1 = new T [nval],
      *v2 = new T [nval];

    T x1 = axis[i1],
      x2 = axis[i2];

    if (ndim == 1) {
        for (int j = 0; j < nval; ++j) {
            v1[j] = data[i1 * nval + j];
            v2[j] = data[i2 * nval + j];
        }
    } else {
        int s = nval;
        for (int j = 1; j < ndim; ++j) s *= len[j];
        _interpn(v1, coor + 1, data + i1 * s, axis + *len, len + 1, ndim - 1, nval, method);
        _interpn(v2, coor + 1, data + i2 * s, axis + *len, len + 1, ndim - 1, nval, method);
    }

    if (method == 0) {
      if (x2 != x1)
          for (int j = 0; j < nval; ++j)
              val[j] = ((*coor - x1) * v2[j] + (x2 - *coor) * v1[j]) / (x2 - x1);
      else
          for (int j = 0; j < nval; ++j)
              val[j] = (v1[j] + v2[j]) / 2.;
    } else if (method == 1) {
      if (x2 != x1)
          for (int j = 0; j < nval; ++j)
              val[j] = exp((log(*coor/x1)*log(v2[j]) + log(x2 / *coor)*log(v1[j]))/log(x2/x1));
      else
          for (int j = 0; j < nval; ++j)
              val[j] = sqrt(v1[j]*v2[j]);
    }

    delete[] v1;
    delete[] v2;
}

#undef  MAX_IT
#define MAX_IT 100

#undef  UNLIKELY_VAL
#define UNLIKELY_VAL -1.11111e+30

template<typename T>
T _interp1(T x, T const *data, T const *axis, int len, int method = 0)
{
  T value;
  _interpn(&value, &x, data, axis, &len, 1, 1, method);
  return value;
}

template<typename R, typename T>
int _root(
        R x1,
        R x2,
        R xacc,
        R *x_root,
        T func
        )
{
  int
    iter,
    compare;
  R 
    fh,fl,fm,fnew,
    s,xh,xl,xm,xnew;
  static char
    name[]="_root";

  std::stringstream msg;

  fl = func(x1);
  fh = func(x2);
  if ((fl > 0. && fh < 0.) || (fl < 0. && fh > 0.)) {
    xl      = x1;
    xh      = x2;
    /* Set *x_root to an unlikely value: */
    *x_root = UNLIKELY_VAL;

    for (iter = 0; iter < MAX_IT; iter++) {
      xm = 0.5*(xl+xh);
      fm = func(xm);
      s  = sqrt(fm*fm-fl*fh);
      if (s == 0.) {
        return 0;
      }
      xnew = xm+(xm-xl)*((fl > fh ? 1. : -1.)*fm/s);

      if (fabs(xnew-*x_root) <= xacc) {
        return 0;
      }
      *x_root = xnew;

      fnew    = func(*x_root);
      if (fnew == 0.) {
        return 0;
      }

      if ((fnew > 0. ? fabs(fm) : -fabs(fm)) != fm) {
        xl = xm;
        fl = fm;
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fl) : -fabs(fl)) != fl) {
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fh) : -fabs(fh)) != fh) {
        xl = *x_root;
        fl = fnew;
      }
      else {
        msg << "### FATAL ERROR in _root function: should never get here" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      if (fabs(xh-xl) <= xacc) {
        return 0;
      }
    }
    msg << "### FATAL ERROR in _root function: exceeded MAX_IT = ";
    msg << MAX_IT;
    msg << "current root calc = ";
    msg << *x_root;
    msg << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  else {
    if (fl == 0.) {
      *x_root = x1;
      return 0;
    }
    if (fh == 0.) {
      *x_root = x2;
      return 0;
    }

    compare = fabs(fl) >= fabs(fh) ? 1 : -1;
    if (compare < 0) {
      return -1;
    }
    else {
      return 1;
    }
  }

  /* Should never get here. */
  msg << "### FATAL ERROR in _root function: should never get here" << std::endl;
  throw std::runtime_error(msg.str().c_str());
}

#undef  MAX_IT
#undef  UNLIKELY_VAL

template<typename T>
class CubicSpline {
public:
	// constructors
	CubicSpline() {}

 	CubicSpline(T *x0, T *y0, int n) :
		x(n), y(n), dx(n), dy(n), a(n), k(0)
	{
  	for (int i = 0; i < n; i++) {
    	x[i] = x0[i];
      y[i] = y0[i];
   	}
   	for (int i = 0; i < n-1; i++) {
      dx[i] = x[i+1] - x[i];
      dy[i] = (y[i+1] - y[i]) / dx[i];
		}
		T *t = new T[n];
		t[0] = 0.0;
   	a[0] = 0.0;
   	a[n-1] = 0.0;
   
   	for (int i = 1; i < n-1; i++) {
    	T piv = 2.0 * (dx[i-1] + dx[i]) - dx[i-1]*t[i-1];
      t[i] = dx[i] / piv;
      a[i] = (dy[i] - dy[i-1] - dx[i-1]*a[i-1]) / piv;
   	}
   	for (int i = n-2; i > 0; i--) {
      a[i] = a[i] - t[i]*a[i+1];
   	}
		delete[] t;
	}

	// destructor
  ~CubicSpline() {}

	// evaluate function
 	T operator()(T xx) {
		int N = a.size();
  	while (xx < x[k]) {
    	k--;
      if (k == -1) {
      	k = 0;
        return y[0] + (dy[0] - dx[0]*a[1]) * (xx - x[0]);
      } 
   	}
		while (xx >= x[k+1]) {
    	k++;
      if (k == N-1) {
      	k = N-2;
        return y[N-1] + (dy[k] + dx[k]*a[k]) * (xx - x[N-1]);
      }
		 }
		 T w = (xx - x[k]) / dx[k];
		 T v = 1.0 - w;
		 return w*y[k+1] + v*y[k] + dx[k]*dx[k]*((w*w - 1.0)*w*a[k+1] + (v*v - 1.0)*v*a[k]);
	}

protected:
	 std::vector<T> x, y, dx, dy, a;
   int k;
};

#endif
