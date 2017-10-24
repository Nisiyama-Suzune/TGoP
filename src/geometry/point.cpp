#include "basic.cpp"

#ifndef GEOMETRY_POINT
#define GEOMETRY_POINT

namespace geometry {

	/*	struct point : defines a point and its various utility.
			point (const double &x, const double &y) gives a point at (x, y).
			It also represents a vector on a 2D plane.
			point unit () const : returns the unit vector of (x, y).
			point rot90 () const :
				returns a point rotated 90 degrees counter-clockwise with respect to the origin.
			point _rot () const : same as above except clockwise.
			point rotate (const double &t) const : returns a point rotated t radian(s) counter-clockwise.
			Operators are mostly vector operations. i.e. vector +, -, *, / and dot/det product.
	*/

	struct point {
		double x, y;
		explicit point (const double &x = 0, const double &y = 0) : x (x), y (y) {}
		double norm () const { return sqrt (x * x + y * y); }
		double norm2 () const { return x * x + y * y; }
		point unit () const { double l = norm (); return point (x / l, y / l); }
		point rot90 () const {return point (-y, x); }
		point _rot90 () const {return point (y, -x); }
		point rot (const double &t) const {
			double c = cos (t), s = sin (t);
			return point (x * c - y * s, x * s + y * c);
		}
	};

	bool operator == (const point &a, const point &b) { return cmp (a.x, b.x) == 0 && cmp (a.y, b.y) == 0; }
	bool operator != (const point &a, const point &b) { return !(a == b); }
	bool operator < (const point &a, const point &b) {
		return (cmp (a.x, b.x) == 0) ? cmp (a.y, b.y) < 0 : cmp (a.x, b.x) < 0;
	}
	point operator - (const point &a) { return point (-a.x, -a.y); }
	point operator + (const point &a, const point &b) { return point (a.x + b.x, a.y + b.y); }
	point operator - (const point &a, const point &b) { return point (a.x - b.x, a.y - b.y); }
	point operator * (const point &a, const double &b) { return point (a.x * b, a.y * b); }
	point operator / (const point &a, const double &b) { return point (a.x / b, a.y / b); }
	double dot (const point &a, const point &b) { return a.x * b.x + a.y * b.y; }
	double det (const point &a, const point &b) { return a.x * b.y - a.y * b.x; }
	double dis (const point &a, const point &b) { return sqrt (sqr (a.x - b.x) + sqr (a.y - b.y)); }

}

#endif

