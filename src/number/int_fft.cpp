#include "basic.cpp"

#ifndef __NUMBER_INT_FFT
#define __NUMBER_INT_FFT

namespace number {

	/*	Fast Fourier transform for integer :
		Usage : int_fft::solve (int a[N], int b[N]) : calculate a[] * b[].
			The result is in a[].
	*/

	template <const int N = 65536, int L = 15, int MOD = 1000003>
	struct int_fft {

		typedef std::complex <double> complex;
		int MASK; complex w[N];

		int_fft () {
			MASK = (1 << L) - 1;
			for (int i = 0; i < N; ++i)
				w[i] = complex (cos(2 * i * PI / N), sin(2 * i * PI / N));
		}

		void FFT (complex p[], int n) {
			for (int i = 1, j = 0; i < n - 1; ++i) {
				for (int s = n; j ^= s >>= 1, ~j & s;);
				if (i < j) swap(p[i], p[j]);
			}
			for (int d = 0; (1 << d) < n; ++d) {
				int m = 1 << d, m2 = m * 2, rm = n >> (d + 1);
				for (int i = 0; i < n; i += m2) {
					for (int j = 0; j < m; ++j) {
						complex &p1 = p[i + j + m], &p2 = p[i + j];
						complex t = w[rm * j] * p1;
						p1 = p2 - t, p2 = p2 + t;
					}
				}
			}
		}

		complex A[N], B[N], C[N], D[N];

		void solve (int a[N], int b[N]) {
			for (int i = 0; i < N; ++i) {
				A[i] = complex (a[i] >> L, a[i] & MASK);
				B[i] = complex (b[i] >> L, b[i] & MASK);
			}
			FFT(A, N), FFT(B, N);
			for (int i = 0; i < N; ++i) {
				int j = (N - i) % N;
				complex da = (A[i] - conj(A[j])) * complex(0, -0.5),
						db = (A[i] + conj(A[j])) * complex(0.5, 0),
						dc = (B[i] - conj(B[j])) * complex(0, -0.5),
						dd = (B[i] + conj(B[j])) * complex(0.5, 0);
				C[j] = da * dd + da * dc * complex(0, 1);
				D[j] = db * dd + db * dc * complex(0, 1);
			}
			FFT(C, N), FFT(D, N);
			for (int i = 0; i < N; ++i) {
				long long da = (long long)(C[i].imag() / N + 0.5) % MOD,
					 db = (long long)(C[i].real() / N + 0.5) % MOD,
					 dc = (long long)(D[i].imag() / N + 0.5) % MOD,
					 dd = (long long)(D[i].real() / N + 0.5) % MOD;
				a[i] = ((dd << (L * 2)) + ((db + dc) << L) + da) % MOD;
			}
		}

	};

}

#endif

