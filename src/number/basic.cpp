#include <bits/stdc++.h>

#ifndef __NUMBER_BASIC
#define __NUMBER_BASIC

namespace number {

	/*	Basic constants & functions :
			long long inverse (const long long &x, const long long &mod) :
				returns the inverse of x modulo mod.
				i.e. x * inv (x) % mod = 1.
			int fpm (int x, int n, int mod) :
				returns x^n % mod. i.e. Fast Power with Modulo.
			void euclid (const long long &a, const long long &b,
			             long long &x, long long &y) :
				solves for ax + by = gcd (a, b).
			long long gcd (const long long &a, const long long &b) :
				solves for the greatest common divisor of a and b.
			long long mul_mod (const long long &a, const long long &b, const long long &mod) :
				returns a * b % mod.
			long long llfpm (const long long &x, const long long &n, const long long &mod) :
				returns x^n % mod.
	*/

	const double PI = acos (-1.);

	long long abs (const long long &x) { return x > 0 ? x : -x; }

	long long inverse (const long long &x, const long long &mod) {
		if (x == 1) return 1;
		return (mod - mod / x) * inverse (mod % x, mod) % mod;
	}

	int fpm (int x, int n, int mod) {
		register int ans = 1, mul = x;
		while (n) {
			if (n & 1) ans = int (1ll * ans * mul % mod);
			mul = int (1ll * mul * mul % mod);
			n >>= 1;
		}
		return ans;
	}

	void euclid (const long long &a, const long long &b,
	             long long &x, long long &y) {
		if (b == 0) x = 1, y = 0;
		else euclid (b, a % b, y, x), y -= a / b * x;
	}

	long long gcd (const long long &a, const long long &b) {
		if (!b) return a;
		long long x = a, y = b;
		while (x > y ? (x = x % y) : (y = y % x));
		return x + y;
	}

	long long mul_mod (const long long &a, const long long &b, const long long &mod) {
		long long d = (long long) floor (a * (double) b / mod + 0.5);
		long long ret = a * b - d * mod;
		if (ret < 0) ret += mod;
		return ret;
	}

	long long llfpm (const long long &x, const long long &n, const long long &mod) {
		long long ans = 1, mul = x, k = n;
		while (k) {
			if (k & 1) ans = mul_mod (ans, mul, mod);
			mul = mul_mod (mul, mul, mod);
			k >>= 1;
		}
		return ans;
	}

}

#endif

