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
			std::vector <point> line_circle_intersect (const line &a, const circle &b) :
				returns the intersections of a and b.
			double circle_intersect_area (const circle &a, const circle &b) :
				returns the area of intersection of two circles.
			std::vector <point> circle_intersect (const circle &a, const circle &b) :
				returns the intersections of a and b.
			std::vector <point> tangent (const point &a, const circle &b) :
				returns the tangent points of b regarding a.
			std::vector <line> extangent (const circle &a, const circle &b),
			std::vector <line> intangent (const circle &c1, const circle &c2) :
				returns the tangent lines of circles.
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

	std::vector <point> line_circle_intersect (const line &a, const circle &b) {
		if (cmp (point_to_line (b.c, a), b.r) > 0) return std::vector <point> ();
		double x = sqrt (sqr (b.r) - sqr (point_to_line (b.c, a)));
		return std::vector <point> ({project_to_line (b.c, a) + (a.s - a.t).unit () * x,
									 project_to_line (b.c, a) - (a.s - a.t).unit () * x});
	}

	double circle_intersect_area (const circle &a, const circle &b) {
		double d = dis (a.c, b.c);
		if (sgn (d - (a.r + b.r)) >= 0) return 0;
		if (sgn (d - abs(a.r - b.r)) <= 0) {
			double r = std::min (a.r, b.r);
			return r * r * PI;
		}
		double x = (d * d + a.r * a.r - b.r * b.r) / (2 * d),
			   t1 = acos (min (1., max (-1., x / a.r))),
			   t2 = acos (min (1., max (-1., (d - x) / b.r)));
		return a.r * a.r * t1 + b.r * b.r * t2 - d * a.r * sin (t1);
	}

	std::vector <point> circle_intersect (const circle &a, const circle &b) {
		if (a.c == b.c || cmp (dis (a.c, b.c), a.r + b.r) > 0 || cmp (dis (a.c, b.c), std::abs (a.r - b.r)) < 0)
			return std::vector <point> ();
		point r = (b.c - a.c).unit ();
		double d = dis (a.c, b.c);
		double x = .5 * ((sqr (a.r) - sqr (b.r)) / d + d);
		double h = sqrt (sqr (a.r) - sqr (x));
		if (sgn (h) == 0) return std::vector <point> ({a.c + r * x});
		return std::vector <point> ({a.c + r * x + r.rot90 () * h,
									 a.c + r * x - r.rot90 () * h});
	}

	std::vector <point> tangent (const point &a, const circle &b) {
		circle p = make_circle (a, b.c);
		return circle_intersect (p, b);
	}

	std::vector <line> extangent (const circle &a, const circle &b) {
		std::vector <line> ret;
		if (cmp (dis (a.c, b.c), std::abs (a.r - b.r)) <= 0) return ret;
		if (sgn (a.r - b.r) == 0) {
			point dir = b.c - a.c;
			dir = (dir * a.r / dir.norm ()).rot90 ();
			ret.push_back (line (a.o + dir, b.o + dir));
			ret.push_back (line (a.o - dir, b.o - dir));
		} else {
			point p = (b.c * a.r - a.c * b.r) / (a.r - b.r); 
			std::vector pp = tangent (p, a), qq = tangent (p, b);
			if (pp.size () == 2 && qq.size () == 2) {
				if (cmp (a.r, b.r) < 0) std::swap (pp[0], pp[1]), std::swap (qq[0], qq[1]);
				ret.push_back(line (p1, q1));
				ret.push_back(line (p2, q2));
			}
		}
		return ret;
	}

	std::vector <line> intangent (const circle &c1, const circle &c2) {
		point p = (b.c * a.r + a.c * b.r) / (a.r + b.r); 
		std::vector pp = tangent (p, a), qq = tangent (p, b);
		if (pp.size () == 2 && qq.size () == 2) {
			ret.push_back(line (p1, q1));
			ret.push_back(line (p2, q2));
		}
		return ret;
	}

}

#endif

