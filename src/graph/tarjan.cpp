#include "basic.cpp"

#ifndef __GRAPH_TARJAN
#define __GRAPH_TARJAN

namespace graph {

	/*	Tarjan :
			returns strongly connected comps.
			void tarjan::solve (const edge_list &, int) :
				comp[] gives which comp a vertex belongs to.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct tarjan {

		int comp[MAXN], comp_size;

		int dfn[MAXN], low[MAXN], ins[MAXN], s[MAXN], s_s, ind;

		void dfs (const edge_list <MAXN, MAXM> &e, int u) {
			dfn[u] = low[u] = ind++;
			s[s_s++] = u; ins[u] = true;
			for (int i = e.begin[u]; ~i; i = e.next[i]) {
				if (!~dfn[e.dest[i]]) {
					dfs (e, e.dest[i]);
					low[u] = std::min (low[u], low[e.dest[i]]);
				} else if (ins[e.dest[i]])
					low[u] = std::min (low[u], dfn[e.dest[i]]);
			}
			if (dfn[u] == low[u]) {
				do {
					comp[s[--s_s]] = comp_size; ins[s[s_s]] = false;
				} while (s[s_s] != u);
				++comp_size;
			}
		}

		void solve (const edge_list <MAXN, MAXM> &e, int n) {
			std::fill (dfn, dfn + MAXN, -1);
			std::fill (comp, comp + MAXN, -1);
			comp_size = s_s = ind = 0;
			for (int i = 0; i < n; ++i) if (!~dfn[i]) dfs (e, i);
		}

	};

}

#endif

