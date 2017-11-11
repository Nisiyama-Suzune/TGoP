#include "basic.cpp"

#ifndef __GRAPH_MAXIMUM_FLOW
#define __GRAPH_MAXIMUM_FLOW

namespace graph {

	/*	Sparse graph maximum flow :
			int isap::solve (flow_edge_list &e, int n, int s, int t) :
				e : edge list.
				n : vertex size.
				s : source.
				t : sink.
	*/

	template <int MAXN = 1000, int MAXM = 100000>
	struct isap {

		int pre[MAXN], d[MAXN], gap[MAXN], cur[MAXN];

		int solve (flow_edge_list <MAXN, MAXM> &e, int n, int s, int t) {
			std::fill (pre, pre + n + 1, 0);
			std::fill (d, d + n + 1, 0);
			std::fill (gap, gap + n + 1, 0);
			for (int i = 0; i < n; i++) cur[i] = e.begin[i];
			gap[0] = n;
			int u = pre[s] = s, v, maxflow = 0;
			while (d[s] < n) {
				v = n;
				for (int i = cur[u]; ~i; i = e.next[i])
					if (e.flow[i] && d[u] == d[e.dest[i]] + 1) {
						v = e.dest[i];
						cur[u] = i;
						break;
					}
				if (v < n) {
					pre[v] = u;
					u = v;
					if (v == t) {
						int dflow = INF, p = t;
						u = s;
						while (p != s) {
							p = pre[p];
							dflow = std::min (dflow, e.flow[cur[p]]);
						}
						maxflow += dflow;
						p = t;
						while (p != s) {
							p = pre[p];
							e.flow[cur[p]] -= dflow;
							e.flow[e.inv[cur[p]]] += dflow;
						}
					}
				} else {
					int mindist = n + 1;
					for (int i = e.begin[u]; ~i; i = e.next[i])
						if (e.flow[i] && mindist > d[e.dest[i]]) {
							mindist = d[e.dest[i]];
							cur[u] = i;
						}
					if (!--gap[d[u]]) return maxflow;
					gap[d[u] = mindist + 1]++;
					u = pre[u];
				}
			}
			return maxflow;
		}

	};

	/*	Dense graph maximum flow :
			int dinic::solve (flow_edge_list &e, int n, int s, int t) :
				e : edge list.
				n : vertex size.
				s : source.
				t : sink.
	*/

	template <int MAXN = 1000, int MAXM = 100000>
	struct dinic {

		int n, s, t;

		int d[MAXN], w[MAXN], q[MAXN];

		int bfs (flow_edge_list <MAXN, MAXM> &e) {
			for (int i = 0; i < n; ++i) d[i] = -1;
			int l, r;
			q[l = r = 0] = s, d[s] = 0;
			for (; l <= r; l ++)
				for (int k = e.begin[q[l]]; k > -1; k = e.next[k])
					if (d[e.dest[k]] == -1 && e.flow[k] > 0) d[e.dest[k]] = d[q[l]] + 1, q[++r] = e.dest[k];
			return d[t] > -1 ? 1 : 0;
		}

		int dfs (flow_edge_list <MAXN, MAXM> &e, int u, int ext) {
			if (u == t) return ext;
			int k = w[u], ret = 0;
			for (; k > -1; k = e.next[k], w[u] = k) {
				if (ext == 0) break;
				if (d[e.dest[k]] == d[u] + 1 && e.flow[k] > 0) {
					int flow = dfs (e, e.dest[k], std::min (e.flow[k], ext));
					if (flow > 0) {
						e.flow[k] -= flow, e.flow[e.inv[k]] += flow;
						ret += flow, ext -= flow;
					}
				}
			}
			if (k == -1) d[u] = -1;
			return ret;
		}

		int solve (flow_edge_list <MAXN, MAXM> &e, int n, int s, int t) {
			int ans = 0;
			dinic::n = n; dinic::s = s; dinic::t = t;
			while (bfs (e)) {
				for (int i = 0; i < n; ++i) w[i] = e.begin[i];
				ans += dfs (e, s, INF);
			}
			return ans;
		}

	};

}

#endif

