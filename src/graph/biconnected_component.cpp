#include "basic.cpp"

#ifndef __GRAPH_BICONNECTED_COMPONENT
#define __GRAPH_BICONNECTED_COMPONENT

namespace graph {

	/*	Vertex biconnected component :
			Divides the edges of an undirected graph into several vertex biconnected components.
			vertex_biconnected_component::solve (const edge_list <MAXN, MAXM> &, int n) :
				comp[] gives the index of the component each edge belongs to.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct vertex_biconnected_component {

		int dfn[MAXN], low[MAXN], comp[MAXN], comp_size;
		std::vector <int> stk;

		int tarjan (int u, int fu) {
			static int tot = 0;
			low[u] = dfn[u] = ++tot;
			for (int i = e.begin[u]; ~i; i = e.next[i]) {
				int v = e.dest[i];
				if (v == fu) continue;
				if (dfn[v] < dfn[u]) stk.push_back (i);
				if (!dfn[v]) {
					low[u] = min (low[u], tarjan (e, v, u));
					if (low[v] >= dfn[u]) {
						int t;
						do {
							t = stk.back (); stk.pop_back ();
							comp[stk.back ()] = comp_size;
						} while (t != i);
						++comp_size;
					}
				} else low[u] = min (low[u], dfn[v]);
			}
			return low[u];
		}

		void solve (const edge_list <MAXN, MAXM> &e, int n) {
			std::fill (dfn, dfn + n, 0);
			std::fill (low, low + n, 0);
			comp_size = 0;
			for (int i = 0; i < n; ++i) if (!dfn[i]) tarjan (e, i, -1);
		}

	};

	/*	Edge biconnected component :
			Divides the vertices of an undirected graph into several edge biconnected components.
			edge_biconnected_component::solve (const edge_list <MAXN, MAXM> &, int n) :
				comp[] gives the index of the component each vertex belongs to.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct edge_biconnected_component {

		int dfn[MAXN], low[MAXN], ins[MAXN], comp[MAXN], comp_size;
		std::vector <int> stk;

		int tarjan (const edge_list <MAXN, MAXM> &e, int u, int from) {
			static int tot = 0;
			low[u] = dfn[u] = ++tot;
			stk.push_back (u);
			for (int i = e.begin[u]; ~i; i = e.next[i]) {
				int v = e.dest[i];
				if (v == from) continue;
				low[u] = min (low[u], dfn[v] ? dfn[v] : tarjan (e, v, u));
			}
			if (low[u] == dfn[u]) {
				int t;
				do {
					t = stk.back(); stk.pop_back ();
					comp[t] = comp_size; ins[t] = false;
				} while (t != u);
				++comp_size;
			}
			return low[u];
		}

		void solve (const edge_list <MAXN, MAXM> &e, int n) {
			std::fill (dfn, dfn + n, 0);
			std::fill (low, low + n, 0);
			comp_size = 0;
			for (int i = 0; i < n; ++i) if (!dfn[i]) tarjan (e, i, -1);
		}

	};

}

#endif

