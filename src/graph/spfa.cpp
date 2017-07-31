#include "basic.cpp"

#ifndef __GRAPH_SPFA
#define __GRAPH_SPFA

namespace graph {

	/*	SPFA :
			Shortest path fast algorithm. (with SLF and LLL)
			bool spfa::solve (const cost_edge_list &e, int n, int s) :
				dist[] gives the distance from s.
				last[] gives the previous vertex.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct spfa {

		int dist[MAXN], last[MAXN];

		int queue[MAXN], cnt[MAXN];
		bool inq[MAXN];

		bool solve (const cost_edge_list <MAXN, MAXM> &e, int n, int s) {
			std::fill (dist, dist + MAXN, INF);
			std::fill (last, last + MAXN, -1);
			std::fill (cnt, cnt + MAXN, 0);
			std::fill (inq, inq + MAXN, false);
			int p = 0, q = 1, size = 1;
			long long avg = 0;
			dist[s] = 0; queue[0] = s; inq[s] = true;
			while (p != q) {
				int n = queue[p]; p = (p + 1) % MAXN;
				if (1LL * dist[n] * size > avg) {
					queue[q] = n;
					q = (q + 1) % MAXN;
					continue;
				}
				inq[n] = false; avg -= dist[n]; --size;
				for (int i = e.begin[n]; ~i; i = e.next[i]) {
					int v = e.dest[i];
					if (dist[v] > dist[n] + e.cost[i]) {
						dist[v] = dist[n] + e.cost[i]; last[v] = n;
						if (!inq[v]) {
							if (++cnt[v] > n) return false;
							inq[v] = true; avg += dist[v]; --size;
							if (dist[v] < dist[queue[p]])
								queue[p = (p + MAXN - 1) % MAXN] = v;
							else {
								queue[q] = v;
								q = (q + 1) % MAXN;
							}
						}
					}
				}
			}
			return true;
		}

	};

}

#endif

