#include "line.cpp"

#ifndef __GEOMETRY_HALF_PLANE_INTERSECTION
#define __GEOMETRY_HALF_PLANE_INTERSECTION

namespace geometry {

	/*	Online half plane intersection (complexity = O(c.size ())) :
			std::vector <point> cut (const std::vector<point> &c, line p) :
				returns the convex polygon cutting convex polygon c with half plane p.
					(left hand with respect to vector p)
				If such polygon does not exist, returns an empty set.
				e.g.
					static const double BOUND = 1e5;
					convex.clear ();
					convex.push_back (point (-BOUND, -BOUND));
					convex.push_back (point (BOUND, -BOUND));
					convex.push_back (point (BOUND, BOUND));
					convex.push_back (point (-BOUND, BOUND));
					convex = cut (convex, line(point, point));
					if (convex.empty ()) { ... }
	*/

	std::vector <point> cut (const std::vector<point> &c, line p) {
		std::vector <point> ret;
		if (c.empty ()) return ret;
		for (int i = 0; i < (int) c.size (); ++i) {
			int j = (i + 1) % (int) c.size ();
			if (turn_left (p.s, p.t, c[i])) ret.push_back (c[i]);
			if (two_side (c[i], c[j], p))
				ret.push_back (line_intersect (p, line (c[i], c[j])));
		}
		return ret;
	}


	/*	Offline half plane intersection (complexity = O(nlogn), n = h.size ()) :
			std::vector <point> half_plane_intersect (std::vector <line> h) :
				returns the intersection of half planes h.
					(left hand with respect to the vector)
				If such polygon does not exist, returns an empty set.
	*/

	bool turn_left (const line &l, const point &p) {
		return turn_left (l.s, l.t, p);
	}

	std::vector <point> half_plane_intersect (std::vector <line> h) {
		typedef std::pair <double, line> polar;
		std::vector <polar> g;
		g.resize (h.size ());
		for (int i = 0; i < (int) h.size (); ++i) {
			point v = h[i].t - h[i].s;
			g[i] = std::make_pair (atan2 (v.y, v.x), h[i]);
		}
		sort (g.begin (), g.end (), [] (const polar &a, const polar &b) {
			if (cmp (a.first, b.first) == 0)
				return sgn (det (a.second.t - a.second.s, b.second.t - a.second.s)) < 0;
			else
				return cmp (a.first, b.first) < 0;
		});
		h.resize (std::unique (g.begin (), g.end (), [] (const polar &a, const polar &b) {
			return cmp (a.first, b.first) == 0;
		}) - g.begin ());
		for (int i = 0; i < (int) h.size (); ++i)
			h[i] = g[i].second;
		int fore = 0, rear = -1;
		std::vector <line> ret;
		for (int i = 0; i < (int) h.size (); ++i) {
			while (fore < rear && !turn_left (h[i], line_intersect (ret[rear - 1], ret[rear]))) {
				--rear;
				ret.pop_back ();
			}
			while (fore < rear && !turn_left (h[i], line_intersect (ret[fore], ret[fore + 1])))
				++fore;
			++rear;
			ret.push_back (h[i]);
		}
		while (rear - fore > 1 && !turn_left (ret[fore], line_intersect (ret[rear - 1], ret[rear]))) {
			--rear;
			ret.pop_back ();
		}
		while (rear - fore > 1 && !turn_left (ret[rear], line_intersect (ret[fore], ret[fore + 1])))
			++fore;
		if (rear - fore < 2) return std::vector <point> ();
		std::vector <point> ans;
		ans.resize (ret.size ());
		for (int i = 0; i < (int) ret.size (); ++i)
			ans[i] = line_intersect (ret[i], ret[ (i + 1) % ret.size ()]);
		return ans;
	}

}

#endif

