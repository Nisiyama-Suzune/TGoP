#include "miller-rabin.cpp"

#ifndef __NUMBER_POLLARD_RHO
#define __NUMBER_POLLARD_RHO

namespace number {

	/*	Pollard Rho :
			std::vector <long long> pollard_rho::solve (const long long &) :
				factorizes an integer.
	*/

	struct pollard_rho {

		miller_rabin is_prime;
		const long long threshold = 13E9;

		long long factorize (const long long &number, const long long &seed) {
			long long x = rand() % (number - 1) + 1, y = x;
			for (int head = 1, tail = 2; ; ) {
				x = mul_mod (x, x, number);
				x = (x + seed) % number;
				if (x == y)
					return number;
				long long answer = gcd (abs (x - y), number);
				if (answer > 1 && answer < number)
					return answer;
				if (++head == tail) {
					y = x;
					tail <<= 1;
				}
			}
		}

		void search (const long long &number, std::vector<long long> &divisor) {
			if (number > 1)  {
				if (is_prime.solve (number))
					divisor.push_back (number);
				else {
					long long factor = number;
					for (; factor >= number;
					        factor = factorize (number, rand () % (number - 1) + 1));
					search (number / factor, divisor);
					search (factor, divisor);
				}
			}
		}

		std::vector <long long> solve (const long long &number) {
			std::vector <long long> ans;
			if (number > threshold)
				search (number, ans);
			else {
				long long rem = number;
				for (long long i = 2; i * i <= rem; ++i)
					while (!(rem % i)) {
						ans.push_back (i);
						rem /= i;
					}
				if (rem > 1) ans.push_back (rem);
			}
			return ans;
		}

	};

}

#endif

