#include "basic.cpp"

#ifndef __NUMBER_MILLER_RABIN
#define __NUMBER_MILLER_RABIN

namespace number {

	/*	Miller Rabin :
			bool miller_rabin::solve (const long long &) :
				tests whether a certain integer is prime.
	*/

	struct miller_rabin {

		int BASE[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

		bool check (const long long &prime, const long long &base) {
			long long number = prime - 1;
			for (; ~number & 1; number >>= 1);
			long long result = llfpm (base, number, prime);
			for (; number != prime - 1 && result != 1 && result != prime - 1; number <<= 1)
				result = mul_mod (result, result, prime);
			return result == prime - 1 || (number & 1) == 1;
		}

		bool solve (const long long &number) {
			if (number < 2) return false;
			if (number < 4) return true;
			if (~number & 1) return false;
			for (int i = 0; i < 12 && BASE[i] < number; ++i)
				if (!check (number, BASE[i]))
					return false;
			return true;
		}

	};

}

#endif

