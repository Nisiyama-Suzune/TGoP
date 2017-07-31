#include <bits/stdc++.h>

#ifndef __GEOMETRY_BASIC
#define __GEOMETRY_BASIC

namespace geometry {

	/*	Constants & basic functions :
			EPS : fixes the possible error of data.
				i.e. x == y iff |x - y| < EPS.
			PI : the value of PI.
			int sgn (const double &x) : returns the sign of x.
			int cmp (const double &x, const double &y) : returns the sign of x - y.
			double sqr (const double &x) : returns x * x.
	*/

	const double EPS = 1E-8;
	const double PI = acos (-1);

	int sgn (const double &x) { return x < -EPS ? -1 : x > EPS; }
	int cmp (const double &x, const double &y) { return sgn (x - y); }
	double sqr (const double &x) { return x * x; }

}

#endif

