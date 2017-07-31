#include "basic.cpp"

#ifndef __GRAPH_DIJKSTRA
#define __GRAPH_DIJKSTRA

namespace graph {

	/*	Dijkstra :
			Shortest path algorithm.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct dijkstra {

		int dist[MAXN], last[MAXN];

		bool vis[MAXN];

		void solve (const cost_edge_list <MAXN, MAXM> &e, int s) {
			std::priority_queue <std::pair <int, int>, std::vector <std::pair <int, int> >,
			    std::greater <std::pair <int, int> > > queue;
			std::fill (dist, dist + MAXN, INF);
			std::fill (last, last + MAXN, -1);
			std::fill (vis, vis + MAXN, false);
			dist[s] = 0;
			queue.push (std::make_pair (0, s));
			while (!queue.empty ()) {
				int n = queue.top ().second; queue.pop (); vis[n] = true;
				for (int i = e.begin[n]; ~i; i = e.next[i]) {
					int v = e.dest[i];
					if (dist[v] > dist[n] + e.cost[i]) {
						dist[v] = dist[n] + e.cost[i]; last[v] = n;
						queue.push (std::make_pair (dist[v], v));
					}
				}
				while (!queue.empty () && vis[queue.top ().second]) queue.pop ();
			}
		}

	};

}

#endif

