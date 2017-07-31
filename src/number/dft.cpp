#include "basic.cpp"

#ifndef __NUMBER_DFT
#define __NUMBER_DFT

namespace number {

	/*	Discrete Fourier transform :
			int dft::init (int n) : 
				initializes the transformation with dimension n.
				Returns the recommended size.
			void dft::solve (complex *a, int n, int f) :
				transforms array a with dimension n to its image representation.
				Transforms back when f = 1. (n should be 2^k)
	*/

	template <int MAXN = 1000000>
	struct dft {

		typedef std::complex <double> complex;

		complex e[2][MAXN];

		int init (int n) {
			int len = 1;
			for (; len <= 2 * n; len <<= 1);
			for (int i = 0; i < len; i++) {
				e[0][i] = complex (cos (2 * PI * i / len), sin (2 * PI * i / len));
				e[1][i] = complex (cos (2 * PI * i / len), -sin (2 * PI * i / len));
			}
			return len;
		}

		void solve (complex *a, int n, int f) {
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

	};

}

#endif

