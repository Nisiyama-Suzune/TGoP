#include "basic.cpp"

#ifndef __NUMBER_BSGS
#define __NUMBER_BSGS

namespace number {

	/*	Baby step giant step algorithm :
			Solves a^x = b (mod c) in O (sqrt (c) * log (sqrt (c))).
			int bsgs::solve (int a, int b, int c) :
				returns -1 when no solution.
	*/

	struct bsgs {

		int solve (int a, int b, int c) {
			std::map <int, int> bs;
			int m = (int) sqrt ((double) c) + 1, res = 1;
			for (int i = 0; i < m; ++i) {
				if (bs.find (res) == bs.end ()) bs[res] = i;
				res = int (1LL * res * a % c);
			}
			int mul = 1, inv = (int) inverse (a, c);
			for (int i = 0; i < m; ++i)
				mul = int (1LL * mul * inv % c);
			res = b % c;
			for (int i = 0; i < m; ++i) {
				if (bs.find (res) != bs.end ())
					return i * m + bs[res];
				res = int (1LL * res * mul % c);
			}
			return -1;
		}

	};

}

#endif

