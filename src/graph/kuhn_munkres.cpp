#include "basic.cpp"

#ifndef __GRAPH_KHUN_MUNKRES
#define __GRAPH_KHUN_MUNKRES

namespace graph {

	/*	Kuhn–Munkres algorithm :
			weighted maximum matching algorithm for bipartition graphs (1-base). Complexity O (N^3).
			struct kuhn_munkres :
				Initialize : pass n as the size of both sets, w as the weight matrix.
				Usage : solve () for the maximum matching. The exact matching is in match[].
	*/

	template <int MAXN = 500>
	struct kuhn_munkres {

		int n, w[MAXN][MAXN];  

		int lx[MAXN], ly[MAXN], match[MAXN], way[MAXN], slack[MAXN];
		bool used[MAXN];

		void hungary(int x) {
			match[0] = x; 
			int j0 = 0;
			std::fill (slack, slack + n + 1, INF);
			std::fill (used, used + n + 1, false);
			do {
				used[j0] = true;
				int i0 = match[j0], delta = INF, j1 = 0;
				for (int j = 1; j <= n; ++j) 
					if (used[j] == false) {
						int cur = -w[i0][j] - lx[i0] - ly[j];
						if (cur < slack[j]) {
							slack[j] = cur;
							way[j] = j0;
						}
						if (slack[j] < delta) {
							delta = slack[j];
							j1 = j;
						}
					}
				for (int j = 0; j <= n; ++j) {
					if (used[j]) {
						lx[match[j]] += delta;
						ly[j] -= delta;
					}
					else slack[j] -= delta;
				}
				j0 = j1;
			} while (match[j0] != 0);
			do {
				int j1 = way[j0];
				match[j0] = match[j1];
				j0 = j1;
			} while (j0);
		}

		int solve() {
			for (int i = 1; i <= n; ++i)
				match[i] = lx[i] = ly[i] = way[i] = 0;
			for (int i = 1; i <= n; ++i) hungary (i);
			int sum = 0;
			for (int i = 1; i <= n; ++i) sum += w[match[i]][i];
			return sum;
		}

	};

}

#endif

