#include "point.cpp"

#ifndef __GEOMETRY_MINIMUM_CIRCLE
#define __GEOMETRY_MINIMUM_CIRCLE

namespace geometry {

	/*	Minimum circle of a point set :
		circle minimum_circle (std::vector <point> p) : returns the minimum circle of point set p.
	*/

	circle minimum_circle (std::vector <point> p) {
		circle ret;
		std::random_shuffle (p.begin (), p.end ());
		for (int i = 0; i < (int) p.size (); ++i)
			if (!in_circle (p[i], ret)) {
				ret = circle (p[i], 0);
				for (int j = 0; j < i; ++j)
					if (!in_circle (p[j], ret)) {
						ret = make_circle (p[j], p[i]);
						for (int k = 0; k < j; ++k)
							if (!in_circle (p[k], ret)) ret = make_circle (p[i], p[j], p[k]);
					}
			}
		return ret;
	}

}

#endif

