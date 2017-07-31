#include "basic.cpp"

#ifndef __GRAPH_MINIMUM_COST_FLOW
#define __GRAPH_MINIMUM_COST_FLOW

namespace graph {

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

}

#endif

