
namespace graph {

	const int MAXN = 1E5, MAXM = 1E5;

	struct edge_list {
		int size;
		int begin[MAXN], dest[MAXM], next[MAXM];
		void clear (int n) {
			size = 0;
			for (int i = 0; i < n; i++)
				begin[i] = -1;
		}
		edge_list (int n = MAXN) {
			clear (n);
		}
		void add_edge (int u, int v) {
			dest[size] = v; next[size] = begin[u]; begin[u] = size++;
		}
	};
}

#include <cstdio>

int main () {
	return 0;
}

