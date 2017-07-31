#include "basic.cpp"

#ifndef __GRAPH_HOPCOFT_KARP
#define __GRAPH_HOPCOFT_KARP

namespace graph {

	/*	Hopcoft-Carp algorithm :
			unweighted maximum matching for bipartition graphs with complexity O (m * n^0.5).
			struct hopcoft_carp :
				Usage : solve () for maximum matching. The matching is in matchx and matchy.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct hopcoft_carp {

		int n, m;

		int matchx[MAXN], matchy[MAXN], level[MAXN];

		bool dfs (edge_list <MAXN, MAXM> &e, int x) {
			for (int i = e.begin[x]; ~i; i = e.next[i]) {
				int y = e.dest[i];
				int w = matchy[y];
				if (w == -1 || (level[x] + 1 == level[w] && dfs (e, w))) {
					matchx[x] = y;
					matchy[y] = x;
					return true;
				}
			}
			level[x] = -1;
			return false;
		}

		int solve (edge_list <MAXN, MAXM> &e, int n, int m) {
			std::fill (matchx, matchx + n, -1);
			std::fill (matchy, matchy + m, -1);
			for (int answer = 0; ; ) {
				std::vector <int> queue;
				for (int i = 0; i < n; ++i) {
					if (matchx[i] == -1) {
						level[i] = 0;
						queue.push_back (i);
					} else {
						level[i] = -1;
					}
				}
				for (int head = 0; head < (int) queue.size(); ++head) {
					int x = queue[head];
					for (int i = e.begin[x]; ~i; i = e.next[i]) {
						int y = e.dest[i];
						int w = matchy[y];
						if (w != -1 && level[w] < 0) {
							level[w] = level[x] + 1;
							queue.push_back (w);
						}
					}
				}
				int delta = 0;
				for (int i = 0; i < n; ++i)
					if (matchx[i] == -1 && dfs (e, i)) delta++;
				if (delta == 0) return answer;
				else answer += delta;
			}
		}

	};

}

#endif

