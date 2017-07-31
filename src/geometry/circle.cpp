#include "line.cpp"

#ifndef __GEOMETRY_CIRCLE
#define __GEOMETRY_CIRCLE

namespace geometry {

	/*	struct circle : defines a circle.
		circle (point c, double r) gives a circle with center c and radius r.
	*/

	struct circle {
		point c;
		double r;
		explicit circle (point c = point (), double r = 0) : c (c), r (r) {}
	};

	bool operator == (const circle &a, const circle &b) {
		return a.c == b.c && cmp (a.r, b.r) == 0;
	}

	bool operator != (const circle &a, const circle &b) {
		return ! (a == b);
	}

	/*	Circle interaction :
			bool in_circle (const point &a, const circle &b) : checks if a is in or on b.
			circle make_circle (const point &a, const point &b) :
				generates a circle with diameter ab.
			circle make_circle (const point &a, const point &b, const point &c) :
				generates a circle passing a, b and c.
			std::pair <point, point> line_circle_intersect (const line &a, const circle &b) :
				returns the intersections of a and b.
				Fails if a and b do not intersect.
			std::pair <point, point> circle_intersect (const circle &a, const circle &b):
				returns the intersections of a and b.
				Fails if a and b do not intersect.
			std::pair <line, line> tangent (const point &a, const circle &b) :
				returns the tangent lines of b passing through a.
				Fails if a is in b.
	*/

	bool in_circle (const point &a, const circle &b) {
		return cmp (dis (a, b.c), b.r) <= 0;
	}


	circle make_circle (const point &a, const point &b) {
		return circle ((a + b) / 2, dis (a, b) / 2);
	}

	circle make_circle (const point &a, const point &b, const point &c) {
		point p = circumcenter (a, b, c);
		return circle (p, dis (p, a));
	}

	std::pair <point, point> line_circle_intersect (const line &a, const circle &b) {
		double x = sqrt (sqr (b.r) - sqr (point_to_line (b.c, a)));
		return std::make_pair (project_to_line (b.c, a) + (a.s - a.t).unit () * x,
		                       project_to_line (b.c, a) - (a.s - a.t).unit () * x);
	}

	point __circle_intersect (const circle &a, const circle &b) {
		point r = (b.c - a.c).unit ();
		double d = dis (a.c, b.c);
		double x = .5 * ((sqr (a.r) - sqr (b.r)) / d + d);
		double h = sqrt (sqr (a.r) - sqr (x));
		return a.c + r * x + r.rot90 () * h;
	}

	std::pair <point, point> circle_intersect (const circle &a, const circle &b) {
		return std::make_pair (__circle_intersect (a, b), __circle_intersect (b, a));
	}

	std::pair <line, line> tangent (const point &a, const circle &b) {
		circle p = make_circle (a, b.c);
		auto d = circle_intersect (p, b);
		return std::make_pair (line (d.first, a), line (d.second, a));
	}

}

#endif

