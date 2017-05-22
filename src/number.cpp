/*	Number theory template :
	Deals with various aspects of integer, division, modulo, etc.
*/

#include <cmath>
#include <complex>

namespace number {

/*	Basic operation :
		long long inverse (const long long &x, const long long &mod) :
			returns the inverse of x modulo mod.
			i.e. x * inv (x) % mod = 1.
		int fpm (int x, int n, int mod) :
			returns x^n % mod. i.e. Fast Power with Modulo.
*/
		
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

	/* Discrete Fourier transform :
		int dft::prepare (int n) : readys the transformation with dimension n.
		void dft::main (complex *a, int n, int f) :
			transforms array a with dimension n to its frequency representation.
			transforms it back when f = 1.
	*/

	namespace dft {

		const int MAXN = 1E6;
		const double PI = acos (-1);

		typedef std::complex <double> complex;

		complex e[2][MAXN];

		int prepare (int n) {
			int len = 1;
			for (; len <= 2 * n; len <<= 1);
			for (int i = 0; i < len; i++) {
				e[0][i] = complex (cos (2 * PI * i / len), sin (2 * PI * i / len));
				e[1][i] = complex (cos (2 * PI * i / len), -sin (2 * PI * i / len));
			}
			return len;
		}

		void main (complex *a, int n, int f) {
			for (int i = 0, j = 0; i < n; i++) {
				if (i > j) std::swap (a[i], a[j]);
				for (int t = n >> 1; (j ^= t) < t; t >>= 1);
			}
			for (int i = 2; i <= n; i <<= 1)
				for (int j = 0; j < n; j += i)
					for (int k = 0; k < (i >> 1); k++) {
						complex A = a[j + k];
						complex B = e[f][n / i * k] * a[j + k + (i >> 1)];
						a[j + k] = A + B;
						a[j + k + (i >> 1)] = A - B;
					}
			if (f == 1) {
				for (int i = 0; i < n; i++)
					a[i] = complex (a[i].real () / n, a[i].imag ());
			}
		}
	}

	/* Number-theoretic transform :
		void ntt::main (int *a, int n, int f, int mod, int prt) :
			converts polynominal f (x) = a[0] * x^0 + a[1] * x^1 + ... + a[n - 1] * x^(n - 1)
				to a vector (f (prt^0), f (prt^1), f (prt^2), ..., f (prt^(n - 1))). (module mod)
			Converts back if f = 1.
			Requries specific mod and corresponding prt to work. (given in MOD and PRT)
		void ntt::crt (int *a, int mod) :
			makes up the results a from 3 primes to a certain module mod.
	*/

	namespace ntt {

		const int MAXN = 1E6;

		void main (int *a, int n, int f, int mod, int prt) {
			for (register int i = 0, j = 0; i < n; i++) {
				if (i > j) std::swap (a[i], a[j]);
				for (register int t = n >> 1; (j ^= t) < t; t >>= 1);
			}
			for (register int i = 2; i <= n; i <<= 1) {
				static int exp[MAXN];
				exp[0] = 1;
				exp[1] = fpm (prt, (mod - 1) / i, mod);
				if (f == 1) exp[1] = fpm (exp[1], mod - 2, mod);
				for (register int k = 2; k < (i >> 1); k++) {
					exp[k] = int (1ll * exp[k - 1] * exp[1] % mod);
				}
				for (register int j = 0; j < n; j += i) {
					for (register int k = 0; k < (i >> 1); k++) {
						register int &pA = a[j + k], &pB = a[j + k + (i >> 1)];
						register int A = pA, B = int (1ll * pB * exp[k] % mod);
						pA = (A + B) % mod;
						pB = (A - B + mod) % mod;
					}
				}
			}
			if (f == 1) {
				register int rev = fpm (n, mod - 2, mod);
				for (register int i = 0; i < n; i++) {
					a[i] = int (1ll * a[i] * rev % mod);
				}
			}
		}

		const int MOD[3] = {1045430273, 1051721729, 1053818881}, PRT[3] = {3, 6, 7};

		int crt (int *a, int mod) {
			static int inv[3][3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					inv[i][j] = (int) inverse (MOD[i], MOD[j]);
			static int x[3];
			for (int i = 0; i < 3; i++) {
				x[i] = a[i];
				for (int j = 0; j < i; j++) {
					int t = (x[i] - x[j] + MOD[i]) % MOD[i];
					if (t < 0) t += MOD[i];
					x[i] = int (1LL * t * inv[j][i] % MOD[i]);
				}
			}
			int sum = 1, ret = x[0] % mod;
			for (int i = 1; i < 3; i ++) {
				sum = int (1LL * sum * MOD[i - 1] % mod);
				ret += int (1LL * x[i] * sum % mod);
				if (ret >= mod) ret -= mod;
			}
			return ret;
		}

	}

}

#include <cstdio>

int main () {
	return 0;
}

