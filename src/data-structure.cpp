#include <bits/stdc++.h>

namespace data_structure {

	/*	KD-tree :
			queries the k-th closest point in O (k * n ^ (1 - 1 / k)).
			Stores the data in p[].
			Call init (n, k).
			Call min_kth (d, k) / max_kth (d, k).
	*/

	template <int MAXN = 200000, int MAXK = 2>
	struct kd_tree {

		int k, size;

		struct point {
			int data[MAXK], id;
		} p[MAXN];

		struct kd_node {
			int l, r;
			point p, dmin, dmax;
			kd_node() {}
			kd_node (const point &rhs) : l (0), r (0), p (rhs), dmin (rhs), dmax (rhs) {}
			inline void merge (const kd_node &rhs, int k) {
				for (register int i = 0; i < k; i++) {
					dmin.data[i] = std::min (dmin.data[i], rhs.dmin.data[i]);
					dmax.data[i] = std::max (dmax.data[i], rhs.dmax.data[i]);
				}
			}
			inline long long min_dist (const point &rhs, int k) const {
				register long long ret = 0;
				for (register int i = 0; i < k; i++) {
					if (dmin.data[i] <= rhs.data[i] && rhs.data[i] <= dmax.data[i]) continue;
					ret += std::min (1ll * (dmin.data[i] - rhs.data[i]) * (dmin.data[i] - rhs.data[i]),
					                 1ll * (dmax.data[i] - rhs.data[i]) * (dmax.data[i] - rhs.data[i]));
				}
				return ret;
			}
			long long max_dist (const point &rhs, int k) {
				long long ret = 0;
				for (register int i = 0; i < k; i++) {
					int tmp = std::max (std::abs (dmin.data[i] - rhs.data[i]),
					                    std::abs (dmax.data[i] - rhs.data[i]));
					ret += 1ll * tmp * tmp;
				}
				return ret;
			}
		} tree[MAXN * 4];

		struct result {
			long long dist;
			point d;
			result() {}
			result (const long long &dist, const point &d) : dist (dist), d (d) {}
			bool operator > (const result &rhs) const {
				return dist > rhs.dist || (dist == rhs.dist && d.id < rhs.d.id);
			}
			bool operator < (const result &rhs) const {
				return dist < rhs.dist || (dist == rhs.dist && d.id > rhs.d.id);
			}
		};

		inline long long sqrdist (const point &a, const point &b) {
			register long long ret = 0;
			for (register int i = 0; i < k; i++) {
				ret += 1ll * (a.data[i] - b.data[i]) * (a.data[i] - b.data[i]);
			}
			return ret;
		}

		inline int alloc() {
			tree[size].l = tree[size].r = 0;
			return size++;
		}

		void build (const int &depth, int &rt, const int &l, const int &r) {
			if (l > r) return;
			register int middle = (l + r) >> 1;
			std::nth_element (p + l, p + middle, p + r + 1,
			[ = ] (const point & a, const point & b) {
				return a.data[depth] < b.data[depth];
			});
			tree[rt = alloc()] = kd_node (p[middle]);
			if (l == r) return;
			build ((depth + 1) % k, tree[rt].l, l, middle - 1);
			build ((depth + 1) % k, tree[rt].r, middle + 1, r);
			if (tree[rt].l) tree[rt].merge (tree[tree[rt].l], k);
			if (tree[rt].r) tree[rt].merge (tree[tree[rt].r], k);
		}

		std::priority_queue<result, std::vector<result>, std::greater<result> > heap_l;
		std::priority_queue<result, std::vector<result>, std::less<result> > heap_r;

		void _min_kth (const int &depth, const int &rt, const int &m, const point &d) {
			result tmp = result (sqrdist (tree[rt].p, d), tree[rt].p);
			if ((int)heap_l.size() < m) {
				heap_l.push (tmp);
			} else if (tmp < heap_l.top()) {
				heap_l.pop();
				heap_l.push (tmp);
			}
			int x = tree[rt].l, y = tree[rt].r;
			if (x != 0 && y != 0 && sqrdist (d, tree[x].p) > sqrdist (d, tree[y].p)) std::swap (x, y);
			if (x != 0 && ((int)heap_l.size() < m || tree[x].min_dist (d, k) < heap_l.top().dist)) {
				_min_kth ((depth + 1) % k, x, m, d);
			}
			if (y != 0 && ((int)heap_l.size() < m || tree[y].min_dist (d, k) < heap_l.top().dist)) {
				_min_kth ((depth + 1) % k, y, m, d);
			}
		}

		void _max_kth (const int &depth, const int &rt, const int &m, const point &d) {
			result tmp = result (sqrdist (tree[rt].p, d), tree[rt].p);
			if ((int)heap_r.size() < m) {
				heap_r.push (tmp);
			} else if (tmp > heap_r.top()) {
				heap_r.pop();
				heap_r.push (tmp);
			}
			int x = tree[rt].l, y = tree[rt].r;
			if (x != 0 && y != 0 && sqrdist (d, tree[x].p) < sqrdist (d, tree[y].p)) std::swap (x, y);
			if (x != 0 && ((int)heap_r.size() < m || tree[x].max_dist (d, k) >= heap_r.top().dist)) {
				_max_kth ((depth + 1) % k, x, m, d);
			}
			if (y != 0 && ((int)heap_r.size() < m || tree[y].max_dist (d, k) >= heap_r.top().dist)) {
				_max_kth ((depth + 1) % k, y, m, d);
			}
		}

		void init (int n, int k) {
			this -> k = k; size = 0;
			int rt = 0;
			build (0, rt, 0, n - 1);
		}

		point min_kth (const point &d, const int &m) {
			heap_l = decltype (heap_l) ();
			_min_kth (0, 0, m, d);
			return heap_l.top ().d;
		}

		point max_kth (const point &d, const int &m) {
			heap_r = decltype (heap_r) ();
			_max_kth (0, 0, m, d);
			return heap_r.top ().d;
		}

	};

	/*	Link Cut Tree :
		Needs formatting.

	struct MsgNode {
		int leftColor, rightColor, answer;
		MsgNode() {
			leftColor = -1;
			rightColor = -1;
			answer = 0;
		}
		MsgNode (int c) {
			leftColor = rightColor = c;
			answer = 1;
		}
		MsgNode operator + (const MsgNode &p)const {
			if (answer == 0) return p;
			if (p.answer == 0) return *this;
			MsgNode ret;
			ret.leftColor = leftColor;
			ret.rightColor = p.rightColor;
			ret.answer = answer + p.answer - (rightColor == p.leftColor);
			return ret;
		}
	} d[MAXN], g[MAXN];
	int n, m, c[MAXN][2], f[MAXN], p[MAXN], s[MAXN], flag[MAXN];
	bool r[MAXN];
	void init (int x, int value) {
		d[x] = g[x] = MsgNode (value);
		c[x][0] = c[x][1] = 0;
		f[x] = p[x] = flag[x] = -1;
		s[x] = 1;
	}
	void update (int x) {
		s[x] = s[c[x][0]] + s[c[x][1]] + 1;
		g[x] = MsgNode();
		if (c[x][0 ^ r[x]]) g[x] = g[x] + g[c[x][0 ^ r[x]]];
		g[x] = g[x] + d[x];
		if (c[x][1 ^ r[x]]) g[x] = g[x] + g[c[x][1 ^ r[x]]];
	}
	void makesame (int x, int c) {
		flag[x] = c;
		d[x] = MsgNode (c);
		g[x] = MsgNode (c);
	}
	void pushdown (int x) {
		if (r[x]) {
			std::swap (c[x][0], c[x][1]);
			r[c[x][0]] ^= 1;
			r[c[x][1]] ^= 1;
			std::swap (g[c[x][0]].leftColor, g[c[x][0]].rightColor);
			std::swap (g[c[x][1]].leftColor, g[c[x][1]].rightColor);
			r[x] = false;
		}
		if (flag[x] != -1) {
			if (c[x][0]) makesame (c[x][0], flag[x]);
			if (c[x][1]) makesame (c[x][1], flag[x]);
			flag[x] = -1;
		}
	}
	void rotate (int x, int k) {
		pushdown (x); pushdown (c[x][k]);
		int y = c[x][k]; c[x][k] = c[y][k ^ 1]; c[y][k ^ 1] = x;
		if (f[x] != -1) c[f[x]][c[f[x]][1] == x] = y;
		f[y] = f[x]; f[x] = y; f[c[x][k]] = x; std::swap (p[x], p[y]);
		update (x); update (y);
	}
	void splay (int x, int s = -1) {
		pushdown (x);
		while (f[x] != s) {
			if (f[f[x]] != s) rotate (f[f[x]], (c[f[f[x]]][1] == f[x]) ^ r[f[f[x]]]);
			rotate (f[x], (c[f[x]][1] == x) ^ r[f[x]]);
		}
		update (x);
	}
	void access (int x) {
		int y = 0;
		while (x != -1) {
			splay (x); pushdown (x);
			f[c[x][1]] = -1; p[c[x][1]] = x;
			c[x][1] = y; f[y] = x; p[y] = -1;
			update (x); x = p[y = x];
		}
	}
	void setroot (int x) {
		access (x);
		splay (x);
		r[x] ^= 1;
		std::swap (g[x].leftColor, g[x].rightColor);
	}
	void link (int x, int y) {
		setroot (x);
		p[x] = y;
	}
	void cut (int x, int y) {
		access (x); splay (y, -1);
		if (p[y] == x) p[y] = -1;
		else {
			access (y);
			splay (x, -1);
			p[x] = -1;
		}
	}

	*/

}

using namespace data_structure;

int main () {
	return 0;
}

