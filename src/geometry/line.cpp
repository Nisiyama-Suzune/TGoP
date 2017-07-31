#include "point.cpp"

#ifndef __GEOMETRY_LINE
#define __GEOMETRY_LINE

namespace geometry {

	/*	struct line : defines a line (segment) based on two points, s and t.
			line (const point &s, const point &t) gives a basic line from s to t.
			double length () const : returns the length of the segment.
	*/

	struct line {
		point s, t;
		explicit line (const point &s = point (), const point &t = point ()) : s (s), t (t) {}
		double length () const { return dis (s, t); }
	};

	/*	Point & line interactions :
			bool point_on_segment (const point &a, const line &b) : checks if a is on b.
			bool intersect_judgement (const line &a, const line &b) : checks if segment a and b intersect.
			point line_intersect (const line &a, const line &b) : returns the intersection of a and b.
				Fails on colinear or parallel situations.
			double point_to_line (const point &a, const line &b) : returns the distance from a to b.
			double point_to_segment (const point &a, const lint &b) : returns the distance from a to b.
				i.e. the minimized length from a to segment b.
			bool in_polygon (const point &p, const std::vector <point> &po) :
				checks if a is in a polygon with vetices po (clockwise or counter-clockwise order).
			double polygon_area (const std::vector <point> &a) :
				returns the signed area of polygon a (positive for counter-clockwise order, and vise-versa).
			point project_to_line (const point &a, const line &b) :
				returns the projection of a on b,
	*/

	bool point_on_segment (const point &a, const line &b) {
		return sgn (det (a - b.s, b.t - b.s)) == 0 && sgn (dot (b.s - a, b.t - a)) <= 0;
	}

	bool two_side (const point &a, const point &b, const line &c) {
		return sgn (det (a - c.s, c.t - c.s)) * sgn (det (b - c.s, c.t - c.s)) < 0;
	}

	bool intersect_judgment (const line &a, const line &b) {
		if (point_on_segment (b.s, a) || point_on_segment (b.t, a)) return true;
		if (point_on_segment (a.s, b) || point_on_segment (a.t, b)) return true;
		return two_side (a.s, a.t, b) && two_side (b.s, b.t, a);
	}

	point line_intersect (const line &a, const line &b) {
		double s1 = det (a.t - a.s, b.s - a.s);
		double s2 = det (a.t - a.s, b.t - a.s);
		return (b.s * s2 - b.t * s1) / (s2 - s1);
	}

	double point_to_line (const point &a, const line &b) {
		return fabs (det (b.t - b.s, a - b.s)) / dis (b.s, b.t);
	}

	point project_to_line (const point &a, const line &b) {
		return b.s + (b.t - b.s) * (dot (a - b.s, b.t - b.s) / (b.t - b.s).norm2 ());
	}

	double point_to_segment (const point &a, const line &b) {
		if (sgn (dot (b.s - a, b.t - b.s) * dot (b.t - a, b.t - b.s)) <= 0)
			return fabs (det (b.t - b.s, a - b.s)) / dis (b.s, b.t);
		return std::min (dis (a, b.s), dis (a, b.t));
	}

	bool in_polygon (const point &p, const std::vector <point> & po) {
		int n = (int) po.size ();
		int counter = 0;
		for (int i = 0; i < n; ++i) {
			point a = po[i], b = po[ (i + 1) % n];
			/*			The following statement checks is p is on the border of the polygon.
						The boolean returned may be changed if necessary.
							i.e. the algorithm may check if p is strictly in the polygon.
			*/
			if (point_on_segment (p, line (a, b))) return true;
			int x = sgn (det (p - a, b - a)), y = sgn (a.y - p.y), z = sgn (b.y - p.y);
			if (x > 0 && y <= 0 && z > 0) counter++;
			if (x < 0 && z <= 0 && y > 0) counter--;
		}
		return counter != 0;
	}

	double polygon_area (const std::vector <point> &a) {
		double ans = 0.0;
		for (int i = 0; i < (int) a.size (); ++i)
			ans += det (a[i], a[ (i + 1) % a.size ()]) / 2.0;
		return ans;
	}

}

#endif

