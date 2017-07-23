/*	Graph template :
		Most algorithms on a graph.
*/

#include <bits/stdc++.h>

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
			dest[size] = u; next[size] = begin[v]; cost[size] = c;
			flow[size] = 0; inv[size] = size - 1; begin[v] = size++;
		}
	};

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

	/*	Weighted matching algorithm :
			maximum matching for general weighted graphs. Not stable.
			struct weighted_match :
				Usage : Set k to the size of vertices, w to the weight matrix.
				Note that k has to be even for the algorithm to work.
	*/

	template <int MAXN = 500>
	struct weighted_match {

		int k;
		long long w[MAXN][MAXN];
		int match[MAXN], path[MAXN], p[MAXN], len;
		long long d[MAXN];
		bool v[MAXN];

		bool dfs (int i) {
			path[len++] = i;
			if (v[i]) return true;
			v[i] = true;
			for (int j = 0; j < k; ++j) {
				if (i != j && match[i] != j && !v[j]) {
					int kok = match[j];
					if (d[kok] < d[i] + w[i][j] - w[j][kok]) {
						d[kok] = d[i] + w[i][j] - w[j][kok];
						if (dfs (kok)) return true;
					}
				}
			}
			--len;
			v[i] = false;
			return false;
		}

		long long solve () {
			if (k & 1) ++k;
			for (int i = 0; i < k; ++i) p[i] = i, match[i] = i ^ 1;
			int cnt = 0;
			for (;;) {
				len = 0;
				bool flag = false;
				std::fill (d, d + k, 0);
				std::fill (v, v + k, 0);
				for (int i = 0; i < k; ++i) {
					if (dfs (p[i])) {
						flag = true;
						int t = match[path[len - 1]], j = len - 2;
						while (path[j] != path[len - 1]) {
							match[t] = path[j];
							std::swap (t, match[path[j]]);
							--j;
						}
						match[t] = path[j];
						match[path[j]] = t;
						break;
					}
				}
				if (!flag) {
					if (++cnt >= 2) break;
					std::random_shuffle (p, p + k);
				}
			}
			long long ans = 0;
			for (int i = 0; i < k; ++i)
				ans += w[i][match[i]];
			return ans / 2;
		}

	};

	/*	Weighted blossom algorithm (vfleaking ver.) :
			maximum matching for general weighted graphs. Complexity O (n^3).
			Note that the base index is 1.
			struct weighted_blossom :
				Usage :
					Set n to the size of the vertices.
					Run init ().
					Set g[][].w to the weight of the edge.
					Run solve ().
					The first result is the answer, the second one is the number of matching pairs.
					Obtain the matching with match[].
	*/

	template <int MAXN = 500>
	struct weighted_blossom {

		struct edge {
			int u, v, w;
			edge (int u = 0, int v = 0, int w = 0): u (u), v (v), w (w) {}
		};
		int n, n_x;
		edge g[MAXN * 2 + 1][MAXN * 2 + 1];
		int lab[MAXN * 2 + 1];
		int match[MAXN * 2 + 1], slack[MAXN * 2 + 1], st[MAXN * 2 + 1], pa[MAXN * 2 + 1];
		int flower_from[MAXN * 2 + 1][MAXN + 1], S[MAXN * 2 + 1], vis[MAXN * 2 + 1];
		std::vector<int> flower[MAXN * 2 + 1];
		std::queue<int> q;
		int e_delta (const edge &e) {
			return lab[e.u] + lab[e.v] - g[e.u][e.v].w * 2;
		}
		void update_slack (int u, int x) {
			if (!slack[x] || e_delta (g[u][x]) < e_delta (g[slack[x]][x]))slack[x] = u;
		}
		void set_slack (int x) {
			slack[x] = 0;
			for (int u = 1; u <= n; ++u)
				if (g[u][x].w > 0 && st[u] != x && S[st[u]] == 0)update_slack (u, x);
		}
		void q_push (int x) {
			if (x <= n)q.push (x);
			else for (size_t i = 0; i < flower[x].size(); i++)q_push (flower[x][i]);
		}
		void set_st (int x, int b) {
			st[x] = b;
			if (x > n)for (size_t i = 0; i < flower[x].size(); ++i)
					set_st (flower[x][i], b);
		}
		int get_pr (int b, int xr) {
			int pr = find (flower[b].begin(), flower[b].end(), xr) - flower[b].begin();
			if (pr % 2 == 1) {
				reverse (flower[b].begin() + 1, flower[b].end());
				return (int)flower[b].size() - pr;
			} else return pr;
		}
		void set_match (int u, int v) {
			match[u] = g[u][v].v;
			if (u > n) {
				edge e = g[u][v];
				int xr = flower_from[u][e.u], pr = get_pr (u, xr);
				for (int i = 0; i < pr; ++i)set_match (flower[u][i], flower[u][i ^ 1]);
				set_match (xr, v);
				rotate (flower[u].begin(), flower[u].begin() + pr, flower[u].end());
			}
		}
		void augment (int u, int v) {
			for (;;) {
				int xnv = st[match[u]];
				set_match (u, v);
				if (!xnv)return;
				set_match (xnv, st[pa[xnv]]);
				u = st[pa[xnv]], v = xnv;
			}
		}
		int get_lca (int u, int v) {
			static int t = 0;
			for (++t; u || v; std::swap (u, v)) {
				if (u == 0)continue;
				if (vis[u] == t)return u;
				vis[u] = t;
				u = st[match[u]];
				if (u)u = st[pa[u]];
			}
			return 0;
		}
		void add_blossom (int u, int lca, int v) {
			int b = n + 1;
			while (b <= n_x && st[b])++b;
			if (b > n_x)++n_x;
			lab[b] = 0, S[b] = 0;
			match[b] = match[lca];
			flower[b].clear();
			flower[b].push_back (lca);
			for (int x = u, y; x != lca; x = st[pa[y]])
				flower[b].push_back (x), flower[b].push_back (y = st[match[x]]), q_push (y);
			reverse (flower[b].begin() + 1, flower[b].end());
			for (int x = v, y; x != lca; x = st[pa[y]])
				flower[b].push_back (x), flower[b].push_back (y = st[match[x]]), q_push (y);
			set_st (b, b);
			for (int x = 1; x <= n_x; ++x)g[b][x].w = g[x][b].w = 0;
			for (int x = 1; x <= n; ++x)flower_from[b][x] = 0;
			for (size_t i = 0; i < flower[b].size(); ++i) {
				int xs = flower[b][i];
				for (int x = 1; x <= n_x; ++x)
					if (g[b][x].w == 0 || e_delta (g[xs][x]) < e_delta (g[b][x]))
						g[b][x] = g[xs][x], g[x][b] = g[x][xs];
				for (int x = 1; x <= n; ++x)
					if (flower_from[xs][x])flower_from[b][x] = xs;
			}
			set_slack (b);
		}
		void expand_blossom (int b) {
			for (size_t i = 0; i < flower[b].size(); ++i)
				set_st (flower[b][i], flower[b][i]);
			int xr = flower_from[b][g[b][pa[b]].u], pr = get_pr (b, xr);
			for (int i = 0; i < pr; i += 2) {
				int xs = flower[b][i], xns = flower[b][i + 1];
				pa[xs] = g[xns][xs].u;
				S[xs] = 1, S[xns] = 0;
				slack[xs] = 0, set_slack (xns);
				q_push (xns);
			}
			S[xr] = 1, pa[xr] = pa[b];
			for (size_t i = pr + 1; i < flower[b].size(); ++i) {
				int xs = flower[b][i];
				S[xs] = -1, set_slack (xs);
			}
			st[b] = 0;
		}
		bool on_found_edge (const edge &e) {
			int u = st[e.u], v = st[e.v];
			if (S[v] == -1) {
				pa[v] = e.u, S[v] = 1;
				int nu = st[match[v]];
				slack[v] = slack[nu] = 0;
				S[nu] = 0, q_push (nu);
			} else if (S[v] == 0) {
				int lca = get_lca (u, v);
				if (!lca)return augment (u, v), augment (v, u), true;
				else add_blossom (u, lca, v);
			}
			return false;
		}
		bool matching() {
			std::fill (S + 1, S + 1 + n_x, -1);
			std::fill (slack + 1, slack + 1 + n_x, -1);
			q = std::queue<int>();
			for (int x = 1; x <= n_x; ++x)
				if (st[x] == x && !match[x])pa[x] = 0, S[x] = 0, q_push (x);
			if (q.empty())return false;
			for (;;) {
				while (q.size()) {
					int u = q.front();
					q.pop();
					if (S[st[u]] == 1) continue;
					for (int v = 1; v <= n; ++v)
						if (g[u][v].w > 0 && st[u] != st[v]) {
							if (e_delta (g[u][v]) == 0) {
								if (on_found_edge (g[u][v])) return true;
							} else update_slack (u, st[v]);
						}
				}
				int d = INF;
				for (int b = n + 1; b <= n_x; ++b)
					if (st[b] == b && S[b] == 1)d = std::min (d, lab[b] / 2);
				for (int x = 1; x <= n_x; ++x)
					if (st[x] == x && slack[x]) {
						if (S[x] == -1)d = std::min (d, e_delta (g[slack[x]][x]));
						else if (S[x] == 0)d = std::min (d, e_delta (g[slack[x]][x]) / 2);
					}
				for (int u = 1; u <= n; ++u) {
					if (S[st[u]] == 0) {
						if (lab[u] <= d)return 0;
						lab[u] -= d;
					} else if (S[st[u]] == 1)lab[u] += d;
				}
				for (int b = n + 1; b <= n_x; ++b)
					if (st[b] == b) {
						if (S[st[b]] == 0)lab[b] += d * 2;
						else if (S[st[b]] == 1)lab[b] -= d * 2;
					}
				q = std::queue<int>();
				for (int x = 1; x <= n_x; ++x)
					if (st[x] == x && slack[x] && st[slack[x]] != x && e_delta (g[slack[x]][x]) == 0)
						if (on_found_edge (g[slack[x]][x]))return true;
				for (int b = n + 1; b <= n_x; ++b)
					if (st[b] == b && S[b] == 1 && lab[b] == 0)expand_blossom (b);
			}
			return false;
		}
		std::pair <long long, int> solve () {
			std::fill (match + 1, match + n + 1, 0);
			n_x = n;
			int n_matches = 0;
			long long tot_weight = 0;
			for (int u = 0; u <= n; ++u)st[u] = u, flower[u].clear();
			int w_max = 0;
			for (int u = 1; u <= n; ++u)
				for (int v = 1; v <= n; ++v) {
					flower_from[u][v] = (u == v ? u : 0);
					w_max = std::max (w_max, g[u][v].w);
				}
			for (int u = 1; u <= n; ++u)lab[u] = w_max;
			while (matching())++n_matches;
			for (int u = 1; u <= n; ++u)
				if (match[u] && match[u] < u)
					tot_weight += g[u][match[u]].w;
			return std::make_pair (tot_weight, n_matches);
		}
		void init () {
			for (int u = 1; u <= n; ++u)
				for (int v = 1; v <= n; ++v)
					g[u][v] = edge (u, v, 0);
		}

	};

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
			for (int i = 0; i < n; i ++) d[i] = -1;
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

		void solve (flow_edge_list <MAXN, MAXM> &e, int n, int s, int t) {
			dinic::n = n; dinic::s = s; dinic::t = t;
			while (bfs (e)) {
				for (int i = 0; i < n; i ++) w[i] = e.begin[i];
				dfs (e, s, INF);
			}
		}

	};

	/*	Sparse graph minimum cost flow :
			std::pair <int, int> minimum_cost_flow::solve (cost_flow_edge_list &e,
			                                               int n, int s, int t) :
				e : edge list.
				n : vertex size.
				s : source.
				t : sink.
				returns the flow and the cost respectively.
	*/


	template <int MAXN = 1000, int MAXM = 100000>
	struct minimum_cost_flow {

		int n, source, target;
		int prev[MAXN];
		int dist[MAXN], occur[MAXN];

		bool augment (cost_flow_edge_list <MAXN, MAXM> &e) {
			std::vector <int> queue;
			std::fill (dist, dist + n, INF);
			std::fill (occur, occur + n, 0);
			dist[source] = 0;
			occur[source] = true;
			queue.push_back (source);
			for (int head = 0; head < (int)queue.size(); ++head) {
				int x = queue[head];
				for (int i = e.begin[x]; ~i; i = e.next[i]) {
					int y = e.dest[i];
					if (e.flow[i] && dist[y] > dist[x] + e.cost[i]) {
						dist[y] = dist[x] + e.cost[i];
						prev[y] = i;
						if (!occur[y]) {
							occur[y] = true;
							queue.push_back (y);
						}
					}
				}
				occur[x] = false;
			}
			return dist[target] < INF;
		}

		std::pair <int, int> solve (cost_flow_edge_list <MAXN, MAXM> &e, int n, int s, int t) {
			minimum_cost_flow::n = n;
			source = s; target = t;
			std::pair <int, int> answer = std::make_pair (0, 0);
			while (augment (e)) {
				int number = INF;
				for (int i = target; i != source; i = e.dest[e.inv[prev[i]]]) {
					number = std::min (number, e.flow[prev[i]]);
				}
				answer.first += number;
				for (int i = target; i != source; i = e.dest[e.inv[prev[i]]]) {
					e.flow[prev[i]] -= number;
					e.flow[e.inv[prev[i]]] += number;
					answer.second += number * e.cost[prev[i]];
				}
			}
			return answer;
		}

	};

	/*	Dense graph minimum cost flow :
			std::pair <int, int> zkw_flow::solve (cost_flow_edge_list &e,
			                                      int n, int s, int t) :
				e : edge list.
				n : vertex size.
				s : source.
				t : sink.
				returns the flow and the cost respectively.
	*/

	template <int MAXN = 1000, int MAXM = 100000>
	struct zkw_flow {

		int n, s, t, totFlow, totCost;
		int dis[MAXN], slack[MAXN], visit[MAXN];

		int modlable() {
			int delta = INF;
			for (int i = 0; i < n; i++) {
				if (!visit[i] && slack[i] < delta) delta = slack[i];
				slack[i] = INF;
			}
			if (delta == INF) return 1;
			for (int i = 0; i < n; i++) if (visit[i]) dis[i] += delta;
			return 0;
		}

		int dfs (cost_flow_edge_list <MAXN, MAXM> &e, int x, int flow) {
			if (x == t) {
				totFlow += flow;
				totCost += flow * (dis[s] - dis[t]);
				return flow;
			}
			visit[x] = 1;
			int left = flow;
			for (int i = e.begin[x]; ~i; i = e.next[i])
				if (e.flow[i] > 0 && !visit[e.dest[i]]) {
					int y = e.dest[i];
					if (dis[y] + e.cost[i] == dis[x]) {
						int delta = dfs (e, y, std::min (left, e.flow[i]));
						e.flow[i] -= delta;
						e.flow[e.inv[i]] += delta;
						left -= delta;
						if (!left) { visit[x] = false; return flow; }
					} else
						slack[y] = std::min (slack[y], dis[y] + e.cost[i] - dis[x]);
				}
			return flow - left;
		}

		std::pair <int, int> solve (cost_flow_edge_list <MAXN, MAXM> &e, int n, int s, int t) {
			zkw_flow::n = n; zkw_flow::s = s; zkw_flow::t = t;
			totFlow = 0; totCost = 0;
			std::fill (dis + 1, dis + t + 1, 0);
			do {
				do {
					std::fill (visit + 1, visit + t + 1, 0);
				} while (dfs (e, s, INF));
			} while (!modlable ());
			return std::make_pair (totFlow, totCost);
		}

	};

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

	/*	Dominator tree :
			void dominator_tree::solve (int s, int n, const edge_list <MAXN, MAXM> &succ) :
				solves for the immediate dominator (idom[]) of each node,
					idom[x] will be x if x does not have a dominator,
					and will be -1 if x is not reachable from s.
	*/

	template <int MAXN = 100000, int MAXM = 100000>
	struct dominator_tree {

		int dfn[MAXN], sdom[MAXN], idom[MAXN], id[MAXN], f[MAXN], fa[MAXN], smin[MAXN], stamp;

		void predfs (int x, const edge_list <MAXN, MAXM> &succ) {
			id[dfn[x] = stamp++] = x;
			for (int i = succ.begin[x]; ~i; i = succ.next[i]) {
				int y = succ.dest[i];
				if (dfn[y] < 0) {
					f[y] = x;
					predfs (y, succ);
				}
			}
		}


		int getfa (int x) {
			if (fa[x] == x) return x;
			int ret = getfa (fa[x]);
			if (dfn[sdom[smin[fa[x]]]] < dfn[sdom[smin[x]]]) {
				smin[x] = smin[fa[x]];
			}
			return fa[x] = ret;
		}

		void solve (int s, int n, const edge_list <MAXN, MAXM> &succ) {
			std::fill (dfn, dfn + n, -1);
			std::fill (idom, idom + n, -1);
			static edge_list <MAXN, MAXM> pred, tmp;
			pred.clear (n);
			for (int i = 0; i < n; ++i)
				for (int j = succ.begin[i]; ~j; j = succ.next[j])
					pred.add_edge (succ.dest[j], i);
			stamp = 0; tmp.clear (n); predfs (s, succ);
			for (int i = 0; i < stamp; ++i)
				fa[id[i]] = smin[id[i]] = id[i];
			for (int o = stamp - 1; o >= 0; --o) {
				int x = id[o];
				if (o) {
					sdom[x] = f[x];
					for (int i = pred.begin[x]; ~i; i = pred.next[i]) {
						int p = pred.dest[i];
						if (dfn[p] < 0) continue;
						if (dfn[p] > dfn[x]) {
							getfa (p);
							p = sdom[smin[p]];
						}
						if (dfn[sdom[x]] > dfn[p]) {
							sdom[x] = p;
						}
					}
					tmp.add_edge (sdom[x], x);
				}
				while (~tmp.begin[x]) {
					int y = tmp.dest[tmp.begin[x]];
					tmp.begin[x] = tmp.next[tmp.begin[x]];
					getfa (y);
					if (x != sdom[smin[y]]) {
						idom[y] = smin[y];
					} else {
						idom[y] = x;
					}
				}
				for (int i = succ.begin[x]; ~i; i = succ.next[i]) {
					if (f[succ.dest[i]] == x) {
						fa[succ.dest[i]] = x;
					}
				}
			}
			idom[s] = s;
			for (int i = 1; i < stamp; ++i) {
				int x = id[i];
				if (idom[x] != sdom[x]) {
					idom[x] = idom[idom[x]];
				}
			}
		}

	};

	/*	Stoer Wagner algorithm :
			Finds the minimum cut of an undirected graph.
			Usage :
				input : n, edge[][] : the size and weights of the graph.
				int stoer_wagner::solve () : returns the minimum cut.
	*/

	template <int MAXN = 500>
	struct stoer_wagner {
		int n, edge[MAXN][MAXN];

		int dist[MAXN];
		bool vis[MAXN], bin[MAXN];

		stoer_wagner () {
			memset (edge, 0, sizeof (edge));
			memset (bin, false, sizeof (bin));
		}

		int contract (int &s, int &t)  {
			memset (dist, 0, sizeof (dist));
			memset (vis, false, sizeof (vis));
			int i, j, k, mincut, maxc;
			for (i = 1; i <= n; i++) {
				k = -1; maxc = -1;
				for (j = 1; j <= n; j++)
					if (!bin[j] && !vis[j] && dist[j] > maxc) {
						k = j;  maxc = dist[j];
					}
				if (k == -1) return mincut;
				s = t; t = k;
				mincut = maxc;
				vis[k] = true;
				for (j = 1; j <= n; j++)
					if (!bin[j] && !vis[j])
						dist[j] += edge[k][j];
			}
			return mincut;
		}

		int solve () {
			int mincut, i, j, s, t, ans;
			for (mincut = INF, i = 1; i < n; i++) {
				ans = contract ( s, t );
				bin[t] = true;
				if (mincut > ans) mincut = ans;
				if (mincut == 0) return 0;
				for (j = 1; j <= n; j++)
					if (!bin[j])
						edge[s][j] = (edge[j][s] += edge[j][t]);
			}
			return mincut;
		}

	};

	/*	Blossom algorithm :
			Maximum match for general graph.
			Usage : int blossom::solve (int n, const edge_list &e).
			The matching is in match[].
	*/

	template <int MAXN = 500, int MAXM = 250000> 
	struct blossom {
		int match[MAXN], d[MAXN], fa[MAXN], c1[MAXN], c2[MAXN], v[MAXN], q[MAXN];
		int *qhead, *qtail;
		struct {
			int fa[MAXN];
			void init (int n) {
				for(int i = 1; i <= n; i++)
					fa[i] = i;
			}
			int find (int x) {
				if(fa[x] != x) fa[x] = find(fa[x]);
				return fa[x];
			}
			void merge (int x, int y) {
				x = find(x);
				y = find(y);
				fa[x] = y;
			}
		} ufs;

		void solve(int x, int y) {
			if(x == y) return;
			if(d[y] == 0) {
				solve(x, fa[fa[y]]);
				match[fa[y]] = fa[fa[y]];
				match[fa[fa[y]]] = fa[y];
			}
			else if(d[y] == 1) {
				solve(match[y], c1[y]);
				solve(x, c2[y]);
				match[c1[y]] = c2[y];
				match[c2[y]] = c1[y];
			}
		}

		int lca (int x, int y, int root) {
			x = ufs.find(x); y = ufs.find(y);
			while (x != y && v[x] != 1 && v[y] != 0) {
				v[x] = 0; v[y] = 1;
				if (x != root) x = ufs.find (fa[x]);
				if (y != root) y = ufs.find (fa[y]);
			}
			if (v[y] == 0) std::swap(x, y);
			for (int i = x; i != y; i = ufs.find (fa[i])) v[i] = -1;
			v[y] = -1;
			return x;
		}

		void contract(int x, int y, int b) {
			for(int i = ufs.find(x); i != b; i = ufs.find(fa[i])) {
				ufs.merge (i, b);
				if(d[i] == 1) {
					c1[i] = x; c2[i] = y;
					*qtail++ = i;
				}
			}
		}

		bool bfs (int root, int n, const edge_list <MAXN, MAXM> &e) {
			ufs.init (n);
			std::fill (d, d + MAXN, -1);
			std::fill (v, v + MAXN, -1);
			qhead = qtail = q;
			d[root] = 0;
			*qtail++ = root;
			while(qhead < qtail) {
				for (int loc = *qhead++, i = e.begin[loc]; ~i; i = e.next[i]) {
					int dest = e.dest[i];
					if(match[dest] == -2 || ufs.find(loc) == ufs.find(dest)) continue;
					if(d[dest] == -1)
						if(match[dest] == -1)
						{
							solve(root, loc);
							match[loc] = dest;
							match[dest] = loc;
							return 1;
						} else {
							fa[dest] = loc; fa[match[dest]] = dest;
							d[dest] = 1; d[match[dest]] = 0;
							*qtail++ = match[dest];
						}
					else if (d[ufs.find(dest)] == 0) {
						int b = lca(loc, dest, root);
						contract(loc, dest, b);
						contract(dest, loc, b);
					}
				}
			}
			return 0;
		}

		int solve (int n, const edge_list <MAXN, MAXM> &e)
		{
			std::fill (fa, fa + n, 0);
			std::fill (c1, c1 + n, 0);
			std::fill (c2, c2 + n, 0);
			std::fill (match, match + n, -1);
			int re = 0;
			for(int i = 0; i < n; i++) 
				if(match[i] == -1)
					if (bfs (i, n, e)) ++re;
					else match[i] = -2;
			return re;
		}

	};

}

using namespace graph;

int main () {
	return 0;
}

