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

	/*	Link-cut Tree :
			Dynamic tree that supports path operations.
			Usage :
				Maintain query values in msg.
				Maintain modifications in tag.
				Change merge () and push () accordingly.
	*/

	template <int MAXN = 100000>
	struct lct {

		struct msg {
			int size;
			explicit msg (int size = 0) : size (size) {}
		};

		struct tag {
			int r;
			explicit tag (int r = 0) : r (r) {}
		};

		struct node {
			int c[2];
			int f, p;
			msg m;
			tag t;
			node () {
				c[0] = c[1] = f = p = -1;
				m = msg ();
				t = tag ();
			}
		} n[MAXN];

		msg merge (const msg &a, const msg &b) {
			return msg (a.size + b.size);
		}

		msg merge (const msg &a, int b) {
			return msg (a.size + 1);
		}
			
		tag merge (const tag &a, const tag &b) {
			return tag (a.r ^ b.r);
		}

		void push (int x, const tag &t) {
			if (t.r) std::swap (n[x].c[0], n[x].c[1]);
			n[x].t = merge (n[x].t, t);
		}

		void update (int x) {
			n[x].m = merge (msg (), x);
			if (~n[x].c[0]) n[x].m = merge (n[x].m, n[n[x].c[0]].m);
			if (~n[x].c[1]) n[x].m = merge (n[x].m, n[n[x].c[1]].m);
		}

		void push_down (int x) {
			if (~n[x].c[0]) push (n[x].c[0], n[x].t);
			if (~n[x].c[1]) push (n[x].c[1], n[x].t);
			n[x].t = tag ();
		}

		void rotate (int x, int k) {
			push_down (x); push_down (n[x].c[k]);
			int y = n[x].c[k]; n[x].c[k] = n[y].c[k ^ 1]; n[y].c[k ^ 1] = x;
			if (n[x].f != -1) n[n[x].f].c[n[n[x].f].c[1] == x] = y;
			n[y].f = n[x].f; n[x].f = y; n[n[x].c[k]].f = x; std::swap (n[x].p, n[y].p);
			update (x); update (y);
		}

		void splay (int x, int s = -1) {
			push_down (x);
			while (n[x].f != s) {
				if (n[n[x].f].f != s) rotate (n[n[x].f].f, n[n[n[x].f].f].c[1] == n[x].f);
				rotate (n[x].f, n[n[x].f].c[1] == x);
			}
			update (x);
		}

		void access (int x) {
			int u = x, v = -1;
			while (u != -1) {
				splay (u); push_down (u);
				if (~n[u].c[1]) n[n[u].c[1]].f = -1, n[n[u].c[1]].p = u;
				n[u].c[1] = v;
				if (~v) n[v].f = u, n[v].p = -1;
				update (u); u = n[v = u].p;
			}
			splay (x);
		}

		void setroot (int x) {
			access (x);
			push (x, tag (1));
		}

		void link (int x, int y) {
			setroot (x);
			n[x].p = y;
		}

		void cut (int x, int y) {
			access (x); splay (y, -1);
			if (n[y].p == x) n[y].p = -1;
			else {
				access (y); splay (x, -1);
				n[x].p = -1;
			}
		}

		void directed_link (int x, int y) {
			access (x);
			n[x].p = y;
		}

		void directed_cut (int x) {
			access (x); 
			if (~n[x].c[0]) n[n[x].c[0]].f = -1;
			n[x].c[0] = -1;
			update (x);
		}

	};

}

using namespace data_structure;

int main () {
	return 0;
}

