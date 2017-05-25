/*	Graph template :
		Most algorithms on a graph.
*/

#include <algorithm>
#include <array>
#include <queue>
#include <vector>

namespace graph {

	const int INF = 1E9;
	
	/*	Edge list:
			Various kinds of edge list.
	*/

	template <int MAXN = 1E5, int MAXM = 1E5>
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

	template <int MAXN = 1E5, int MAXM = 1E5>
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

	template <int MAXN = 1E5, int MAXM = 1E5>
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

	/*	Hopcoft-Carp algorithm :
			maximum matching with complexity O (m * n^0.5).
			struct hopcoft_carp :
				Initialize : pass n, m as the size of both vertex sets, e as the reference of the edge_list.
				Usage : solve () for maximum matching. The matching is in matchx and matchy.
	*/

	template <int MAXN = 1E5, int MAXM = 1E5>
	struct hopcoft_carp {

		int n, m;
		edge_list <MAXN, MAXM> &e;

		int matchx[MAXN], matchy[MAXN], level[MAXN];

		bool dfs (int x) {
			for (int i = e.begin[x]; ~i; i = e.next[i]) {
				int y = e.dest[i];
				int w = matchy[y];
				if (w == -1 || (level[x] + 1 == level[w] && dfs (w))) {
					matchx[x] = y;
					matchy[y] = x;
					return true;
				}
			}
			level[x] = -1;
			return false;
		}

		int solve () {
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
					if (matchx[i] == -1 && dfs (i)) delta++;
				if (delta == 0) return answer;
				else answer += delta;
			}
		}

	};

	/*	Kuhn–Munkres algorithm :
			weighted maximum matching algorithm. Complexity O (N^3).
			struct kuhn_munkres :
				Initialize : pass nx, ny as the size of both sets, w as the weight matrix.
				Usage : solve () for the minimum matching. The matching is in link.
	*/

	template <int MAXN = 500> 
	struct kuhn_munkres {

		int nx, ny;
		int w[MAXN][MAXN];

		int lx[MAXN], ly[MAXN], visx[MAXN], visy[MAXN], slack[MAXN], link[MAXN];

		int dfs (int x) {
			visx[x] = 1;
			for (int y = 0; y < ny; y ++) {
				if (visy[y]) continue;
				int t = lx[x] + ly[y] - w[x][y];
				if (t == 0) {
					visy[y] = 1;
					if (link[y] == -1 || dfs (link[y])) {
						link[y] = x;
						return 1;
					}
				} else slack[y] = std::max (slack[y], t);
			}
			return 0;
		}

		int solve () {
			int i, j;
			std::fill (link, link + ny, -1);
			std::fill (ly, ly + ny, 0);
			for (i = 0; i < nx; i++)
				for (j = 0, lx[i] = INF; j < ny; j++)
					lx[i] = std::min (lx[i], w[i][j]);
			for (int x = 0; x < nx; x++) {
				for (i = 0; i < ny; i++) slack[i] = -INF;
				while (true) {
					std::fill (visx, visx + nx, 0);
					std::fill (visy, visy + ny, 0);
					if (dfs (x)) break;
					int d = INF;
					for (i = 0; i < ny; i++)
						if (!visy[i] && d < slack[i]) d = slack[i];
					for (i = 0; i < nx; i++)
						if (visx[i]) lx[i] -= d;
					for (i = 0; i < ny; i++)
						if (visy[i]) ly[i] += d;
						else slack[i] -= d;
				}
			}
			int res = 0;
			for (i = 0; i < ny; i ++)
				if (link[i] > -1) res += w[link[i]][i];
			return res;
		}

	};

	/*	Weighted matching algorithm :
			maximum match for graphs. Not stable.
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
					if (++cnt >= 1) break;
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
			maximum match for graphs. Complexity O (n^3).
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
					if (S[st[u]] == 1)continue;
					for (int v = 1; v <= n; ++v)
						if (g[u][v].w > 0 && st[u] != st[v]) {
							if (e_delta (g[u][v]) == 0) {
								if (on_found_edge (g[u][v]))return true;
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

	template <int MAXN = 1E3, int MAXM = 1E5>
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

	template <int MAXN = 1E3, int MAXM = 1E5>
	struct dinic {

		flow_edge_list <MAXN, MAXM> &e;
		int n, s, t;

		int d[MAXN], w[MAXN], q[MAXN];

		int bfs() {
			for (int i = 0; i < n; i ++) d[i] = -1;
			int l, r;
			q[l = r = 0] = s, d[s] = 0;
			for (; l <= r; l ++)
				for (int k = e.begin[q[l]]; k > -1; k = e.next[k])
					if (d[e.dest[k]] == -1 && e.flow[k] > 0) d[e.dest[k]] = d[q[l]] + 1, q[++r] = e.dest[k];
			return d[t] > -1 ? 1 : 0;
		}

		int dfs (int u, int ext) {
			if (u == t) return ext;
			int k = w[u], ret = 0;
			for (; k > -1; k = e.next[k], w[u] = k) {
				if (ext == 0) break;
				if (d[e.dest[k]] == d[u] + 1 && e.flow[k] > 0) {
					int flow = dfs (e.dest[k], std::min (e.flow[k], ext));
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
			dinic::e = e; dinic::n = n; dinic::s = s; dinic::t = t;
			while (bfs ()) {
				for (int i = 0; i < n; i ++) w[i] = e.begin[i];
				dfs (s, INF);
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


	template <int MAXN = 1E3, int MAXM = 1E5>
	struct minimum_cost_flow {

		cost_flow_edge_list <MAXN, MAXM> &e;
		int n, source, target;
		int prev[MAXN];
		int dist[MAXN], occur[MAXN];

		bool augment() {
			std::vector <int> queue;
			std::fill (dist, dist + n, INT_MAX);
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
			return dist[target] < INT_MAX;
		}

		std::pair <int, int> solve (cost_flow_edge_list <MAXN, MAXM> &e, int n, int s, int t) {
			minimum_cost_flow::e = e; minimum_cost_flow::n = n;
			source = s; target = t;
			std::pair <int, int> answer = std::make_pair (0, 0);
			while (augment()) {
				int number = INT_MAX;
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

	template <int MAXN = 1E3, int MAXM = 1E5>
	struct zkw_flow {

		cost_flow_edge_list <MAXN, MAXM> &e;
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

		int dfs (int x, int flow) {
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
						int delta = dfs (y, std::min (left, e.flow[i]));
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
			zkw_flow::e = e; zkw_flow::n = n; zkw_flow::s = s; zkw_flow::t = t;
			totFlow = 0; totCost = 0;
			std::fill (dis + 1, dis + t + 1, 0);
			do {
				do {
					std::fill (visit + 1, visit + t + 1, 0);
				} while (dfs (s, INF));
			} while (!modlable ());
			return std::make_pair (totFlow, totCost);
		}

	};

}

#include <cstdio>

using namespace graph;

int main () {
	return 0;
}

