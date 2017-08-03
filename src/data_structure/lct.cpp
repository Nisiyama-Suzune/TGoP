#include <bits/stdc++.h>

#ifndef __DATA_STRUCTURE_LCT
#define __DATA_STRUCTURE_LCT

namespace data_structure {

	/*	Link-cut Tree :
			Dynamic tree that supports path operations.
			Usage :
				Maintain query values in msg.
				Maintain modifications in tag.
				Edit merge (), gen () and push () accordingly.
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

		msg merge (const msg &a, const msg &b) {	// Merge two messages.
			return msg (a.size + b.size);
		}

		msg gen (int a) {
			return msg (1);
		}
			
		tag merge (const tag &a, const tag &b) {	// Merge two tags.
			return tag (a.r ^ b.r);
		}

		void push (int x, const tag &t) {
			if (t.r) std::swap (n[x].c[0], n[x].c[1]);
			n[x].t = merge (n[x].t, t);	// Remember to update messages manually.
		}

		void update (int x) {
			n[x].m = gen (x);
			if (~n[x].c[0]) n[x].m = merge (n[n[x].c[0]].m, n[x].m);
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
			n[y].f = n[x].f; n[x].f = y; if (~n[x].c[k]) n[n[x].c[k]].f = x;
			std::swap (n[x].p, n[y].p);
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

#endif

