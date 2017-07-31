#include "basic.cpp"

#ifndef __GRAPH_TARJAN
#define __GRAPH_TARJAN

namespace graph {

	/*	Tarjan :
			returns strongly connected components.
			void tarjan::solve (const edge_list &, int) :
				component[] gives which component a vertex belongs to.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct tarjan {

		int component[MAXN], component_size;

		int dfn[MAXN], low[MAXN], ins[MAXN], s[MAXN], s_s, ind;

		void dfs (const edge_list <MAXN, MAXM> &e, int u) {
			dfn[u] = low[u] = ind++;
			s[s_s++] = u;
			for (int i = e.begin[u]; ~i; i = e.next[i]) {
				if (!~dfn[e.dest[i]]) {
					dfs (e, e.dest[i]);
					low[u] = std::min (low[u], low[e.dest[i]]);
				} else if (ins[e.dest[i]])
					low[u] = std::min (low[u], dfn[e.dest[i]]);
			}
			if (dfn[u] == low[u]) {
				do {
					component[s[--s_s]] = component_size; ins[s[s_s]] = false;
				} while (s[s_s] != u);
				++component_size;
			}
		}

		void solve (const edge_list <MAXN, MAXM> &e, int n) {
			std::fill (dfn, dfn + MAXN, -1);
			std::fill (component, component + MAXN, -1);
			component_size = s_s = ind = 0;
			for (int i = 0; i < n; ++i) if (!~dfn[i]) dfs (e, i);
		}

	};

}

#endif

