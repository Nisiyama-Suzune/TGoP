#include "point.cpp"

#ifndef __GEOMETRY_CENTER
#define __GEOMETRY_CENTER

namespace geometry {

	/*	Centers of a triangle :
			returns various centers of a triangle with vertices (a, b, c).
	*/

	point incenter (const point &a, const point &b, const point &c) {
		double p = dis (a, b) + dis (b, c) + dis (c, a);
		return (a * dis (b, c) + b * dis (c, a) + c * dis (a, b)) / p;
	}

	point circumcenter (const point &a, const point &b, const point &c) {
		point p = b - a, q = c - a, s (dot (p, p) / 2, dot (q, q) / 2);
		double d = det (p, q);
		return a + point (det (s, point (p.y, q.y)), det (point (p.x, q.x), s)) / d;
	}

	point orthocenter (const point &a, const point &b, const point &c) {
		return a + b + c - circumcenter (a, b, c) * 2.0;
	}

}

#endif

