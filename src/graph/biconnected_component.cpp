#include "basic.cpp"

#ifndef __GRAPH_BICONNECTED_COMPONENT
#define __GRAPH_BICONNECTED_COMPONENT

namespace graph {

	/*	Vertex biconnected component :
			Divides the edges of an undirected graph into several vertex biconnected components.
			vertex_biconnected_component::solve (const edge_list <MAXN, MAXM> &, int n) :
				component[] gives the index of the component each edge belongs to.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct vertex_biconnected_component {

		int component[MAXM], component_size;

		int dfn[MAXN], low[MAXN], s[MAXN], s_s, ind;

		void dfs (const edge_list <MAXN, MAXM> &e, int u, int f) {
			dfn[u] = low[u] = ind++;
			for (int i = e.begin[u]; ~i; i = e.next[i])
				if (e.dest[i] != f && dfn[u] < dfn[e.dest[i]]) {
					s[s_s++] = i;
					if (!~dfn[u]) {
						dfs (e, e.dest[i], u);
						low[u] = std::min (low[u], low[e.dest[i]]);
						if (low[e.dest[i]] >= dfn[u]) {
							do {
								component[s[--s_s]] = component_size;
							} while (component[s[s_s]] != i);
							component_size++;
						}
					} else
						low[u] = std::min (low[u], dfn[e.dest[i]]);
				}
		}

		void solve (const edge_list <MAXN, MAXM> &e, int n) {
			std::fill (dfn, dfn + MAXN, -1);
			std::fill (component, component + MAXM, -1);
			component_size = s_s = ind = 0;
			for (int i = 0; i < n; ++i) if (!~dfn[i]) dfs (e, i, -1);
		}

	};

	/*	Edge biconnected component :
			Divides the vertices of an undirected graph into several edge biconnected components.
			edge_biconnected_component::solve (const edge_list <MAXN, MAXM> &, int n) :
				component[] gives the index of the component each vertex belongs to.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct edge_biconnected_component {

		int component[MAXN], component_size;

		int dfn[MAXN], low[MAXN], s[MAXN], s_s, ind;

		void dfs (const edge_list <MAXN, MAXM> &e, int u, int f) {
			dfn[u] = low[u] = ind++;
			s[s_s++] = u;
			for (int i = e.begin[u]; ~i; i = e.next[i])
				if (e.dest[i] != f) {
					if (!~dfn[e.dest[i]]) {
						dfs (e, e.dest[i], u);
						low[u] = std::min (low[u], low[e.dest[i]]);
						if (low[e.dest[i]] > dfn[u]) {
							do {
								component[s[--s_s]] = component_size;
							} while (s[s_s] != e.dest[i]);
							component_size++;
						}
					} else low[u] = std::min (low[u], dfn[e.dest[i]]);
				}
		}

		void solve (const edge_list <MAXN, MAXM> &e, int n) {
			std::fill (dfn, dfn + MAXN, -1);
			std::fill (component, component + MAXN, -1);
			component_size = s_s = ind = 0;
			for (int i = 0; i < n; ++i) if (!~dfn[i]) dfs (e, i, -1);
		}

	};

}

#endif

