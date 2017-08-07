#include <bits/stdc++.h>

#ifndef __DATA_STRUCTURE_SPLAY_TREE
#define __DATA_STRUCTURE_SPLAY_TREE

namespace data_structure {

	/*	Splay Tree :
			Solver for sequence problems.
			Maintain msg and tag, and update accordingly.
			This sample is collected from BZOJ 1500.
	*/

	const int INF = 1E9;

	template <int MAXN = 510000>
		struct splay_tree {

			//	TODO : Maintain messages here.

			struct msg {
				int size;
				int l_max;
				int r_max;
				int sum_max;
				int sum;
				explicit msg (int size = 0, int l_max = 0, int r_max = 0, int sum_max = 0, int sum = 0) :
					size (size), l_max (l_max), r_max (r_max), sum_max (sum_max), sum (sum) {}
			};

			//	TODO : Maintain tags here.

			struct tag {
				bool r;
				int mod;
				explicit tag (bool r = false, int mod = INF) : r (r), mod (mod) {}
			};

			struct node {
				int c[2], f;
				msg m;
				tag t;
				node () {
					c[0] = c[1] = f = -1;
					m = msg ();
					t = tag ();
				}

				//	TODO : Maintain node values here.

				int val;

			} n[MAXN];

			int root;

			//	TODO : msg & tag operations. 

			msg merge (const msg &a, const msg &b) {
				return msg (a.size + b.size,
						std::max (a.l_max, a.sum + b.l_max),
						std::max (b.r_max, a.r_max + b.sum),
						std::max (std::max (a.sum_max, b.sum_max), a.r_max + b.l_max),
						a.sum + b.sum);
			}

			msg gen (int a) {
				return msg (1, n[a].val, n[a].val, n[a].val, n[a].val);
			}

			tag merge (const tag &a, const tag &b) {
				if (b.mod != INF) return tag (a.r ^ b.r, b.mod);
				return tag (a.r ^ b.r, a.mod);
			}

			void push (int x, const tag &t) {
				if (t.mod != INF) {
					n[x].val = t.mod;
					n[x].m.l_max = (t.mod >= 0 ? t.mod * n[x].m.size : t.mod);
					n[x].m.r_max = (t.mod >= 0 ? t.mod * n[x].m.size : t.mod);
					n[x].m.sum_max = (t.mod >= 0 ? t.mod * n[x].m.size : t.mod);
					n[x].m.sum = t.mod * n[x].m.size;
				}
				if (t.r) {
					std::swap (n[x].c[0], n[x].c[1]);
					std::swap (n[x].m.l_max, n[x].m.r_max);
				}
				n[x].t = merge (n[x].t, t); 
			}

			//	Splay tree operations.

			void push_down (int x) {
				if (~n[x].c[0]) push (n[x].c[0], n[x].t);
				if (~n[x].c[1]) push (n[x].c[1], n[x].t);
				n[x].t = tag ();
			}

			void update (int x) {
				n[x].m = gen (x);
				if (~n[x].c[0]) n[x].m = merge (n[n[x].c[0]].m, n[x].m);
				if (~n[x].c[1]) n[x].m = merge (n[x].m, n[n[x].c[1]].m);
			}

			void rotate (int x, int k) {
				push_down (x); push_down (n[x].c[k]);
				int y = n[x].c[k]; n[x].c[k] = n[y].c[k ^ 1]; n[y].c[k ^ 1] = x;
				if (n[x].f != -1) n[n[x].f].c[n[n[x].f].c[1] == x] = y;
				n[y].f = n[x].f; n[x].f = y; if (~n[x].c[k]) n[n[x].c[k]].f = x; 
				update (x); update (y);
			}

			void splay (int x, int s = -1) {
				push_down (x);
				while (n[x].f != s) {
					if (n[n[x].f].f != s) rotate (n[n[x].f].f, n[n[n[x].f].f].c[1] == n[x].f);
					rotate (n[x].f, n[n[x].f].c[1] == x);
				}
				update (x);
				if (s == -1) root = x;
			}

			//	Alloc & free.

			int buf[MAXN], size;

			splay_tree () {
				root = -1; size = MAXN;
				for (int i = 0; i < MAXN; ++i) buf[i] = i;
			}

			int alloc () {
				while (size == 0);
				n[buf[--size]] = node ();
				return buf[size];
			}

			void del (int x) {
				buf[size++] = x;
			}

			void del_tree (int x) {
				if (!~x) return;
				del (x);
				del_tree (n[x].c[0]);
				del_tree (n[x].c[1]);
			}

			//	TODO : Put your own operations here. REMEMBER TO PUSH_DOWN!!

		};

}

#endif

