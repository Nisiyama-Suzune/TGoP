/*	Formula template :
		Various formulas for programming.
*/

#include <bits/stdc++.h>

namespace formula {

	/*	Zeller's congruence :
			converts between a specific date and its index. 
			(y must be greater than or equal to 0 in all cases.)
			(0 = Monday, 1 = Tuesday, ..., 6 = Sunday) (mod 7)
	*/

	int get_id (int y, int m, int d) {
		if (m < 3) { --y; m += 12; }
		return 365 * y + y / 4 - y / 100 + y / 400 + (153 * (m - 3) + 2) / 5 + d - 307;
	}

	std::tuple <int, int, int> date (int id) {
		int x = id + 1789995, n, i, j, y, m, d;
		n = 4 * x / 146097; x -= (146097 * n + 3) / 4;
		i = (4000 * (x + 1)) / 1461001; x -= 1461 * i / 4 - 31;
		j = 80 * x / 2447; d = x - 2447 * j / 80;
		x = j / 11;
		m = j + 2 - 12 * x; y = 100 * (n - 49) + i + x;
		return std::make_tuple (y, m, d);
	}

	/*	Lattice points below segment :
			solves for sigma [(a + b * i) / m] where 0 <= i < n.
	*/

	long long solve(long long n, long long a, long long b, long long m){
		if (b == 0) return n * (a / m);
		if (a >= m) return n * (a / m) + solve (n, a % m, b, m);
		if (b >= m) return (n - 1) * n / 2 * (b / m) + solve (n, a, b % m, m);
		return solve ((a + b * n) / m, (a + b * n) % m, m, b);
	}

}

using namespace formula;

int main () {
	return 0;
}

