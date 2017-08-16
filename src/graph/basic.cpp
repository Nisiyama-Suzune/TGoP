#include <bits/stdc++.h>

#ifndef __GRAPH_BASIC
#define __GRAPH_BASIC

namespace graph {

	const int INF = 1E9;

	/*	Edge list:
			Various kinds of edge list.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct edge_list {
		int size;
		int begin[MAXN], dest[MAXM], next[MAXM];
		void clear (int n) {
			size = 0;
			std::fill (begin, begin + n, -1);
		}
		edge_list (int n = MAXN) {
			clear (n);
		}
		void add_edge (int u, int v) {
			dest[size] = v; next[size] = begin[u]; begin[u] = size++;
		}
	};

	template <int MAXN = 100000, int MAXM = 100000>
	struct cost_edge_list {
		int size;
		int begin[MAXN], dest[MAXM], next[MAXM], cost[MAXM];
		void clear (int n) {
			size = 0;
			std::fill (begin, begin + n, -1);
		}
		cost_edge_list (int n = MAXN) {
			clear (n);
		}
		void add_edge (int u, int v, int c) {
			dest[size] = v; next[size] = begin[u]; cost[size] = c; begin[u] = size++;
		}
	};

	template <int MAXN = 100000, int MAXM = 100000>
	struct flow_edge_list {
		int size;
		int begin[MAXN], dest[MAXM], next[MAXM], flow[MAXM], inv[MAXM];
		void clear (int n) {
			size = 0;
			std::fill (begin, begin + n, -1);
		}
		flow_edge_list (int n = MAXN) {
			clear (n);
		}
		void add_edge (int u, int v, int f) {
			dest[size] = v; next[size] = begin[u]; flow[size] = f; inv[size] = size + 1; begin[u] = size++;
			dest[size] = u; next[size] = begin[v]; flow[size] = 0; inv[size] = size - 1; begin[v] = size++;
		}
	};

	template <int MAXN = 100000, int MAXM = 100000>
	struct cost_flow_edge_list {
		int size;
		int begin[MAXN], dest[MAXM], next[MAXM], cost[MAXM], flow[MAXM], inv[MAXM];
		void clear (int n) {
			size = 0;
			std::fill (begin, begin + n, -1);
		}
		cost_flow_edge_list (int n = MAXN) {
			clear (n);
		}
		void add_edge (int u, int v, int c, int f) {
			dest[size] = v; next[size] = begin[u]; cost[size] = c;
			flow[size] = f; inv[size] = size + 1; begin[u] = size++;
			dest[size] = u; next[size] = begin[v]; cost[size] = -c;
			flow[size] = 0; inv[size] = size - 1; begin[v] = size++;
		}
	};

}

#endif

