/*	String algorithm :
		Algorithms regarding string.
*/

#include <bits/stdc++.h>

namespace string {

	/*	KMP algorithm :
			void kmp::build (const std::string & str) :
				initializes and builds the failure array. Complexity O (n).
			int kmp::find (const std::string &str) :
				finds the first occurence of match in str. Complexity O (n).

			Note : match is cylic when L % (L - 1 - fail[L - 1]) == 0 &&
				L / (L - 1 - fail[L - 1]) > 1, where L = match.size ().
	*/

	template <int MAXN = 1E6>
	struct kmp {

		std::string match;
		int fail[MAXN];

		void build (const std::string &str) {
			match = str; fail[0] = -1;
			for (int i = 1; i < (int) str.size (); ++i) {
				int j = fail[i - 1];
				while (~j && str[i] != str[j + 1]) j = fail[j];
				fail[i] = str[i] == str[j + 1] ? j + 1 : -1;
			}
		}

		int find (const std::string &str) {
			for (int i = 0, j = 0; i < (int) str.size (); ++i, ++j) {
				if (j == match.size ()) return i - match.size ();
				while (~j && str[i] != match[j]) j = fail[j];
			}
			return str.size ();
		}

	};

	/*	Suffix automaton :
			void suffix_automaton::init () :
				initializes the automaton with an empty string.
			void suffix_automaton::extend (int token) :
				extends the string with token. Complexity O (1).

			head : the first state.
			tail : the last state.
				Terminating states can be reached via visiting the ancestors of tail.
			state::len : the longest length of the string in the state.
			state::parent : the parent link.
			state::dest : the automaton link.
	*/

	template <int MAXN = 1E6, int MAXC = 26>
	struct suffix_automaton {

		state *head, *tail;

		struct state {
			int len;
			state *parent, *dest[MAXC];
			state (int len = 0) : len (len), parent (NULL) {
				memset (dest, 0, sizeof (dest));
			}
			state *extend (state *, int token);
		} node_pool[MAXN * 2], *tot_node, *null = new state();

		state *state::extend (state *start, int token) {
			state *p = this;
			state *np = this -> dest[token] ? null : new (tot_node++) state (this -> len + 1);
			while (p && !p -> dest[token])
				p -> dest[token] = np, p = p -> parent;
			if (!p) np -> parent = start;
			else {
				state *q = p -> dest[token];
				if (p -> len + 1 == q -> len) {
					np -> parent = q;
				} else {
					state *nq = new (tot_node++) state (*q);
					nq -> len = p -> len + 1;
					np -> parent = q -> parent = nq;
					while (p && p -> dest[token] == q) {
						p -> dest[token] = nq, p = p -> parent;
					}
				}
			}
			return np == null ? np -> parent : np;
		}

		void init () {
			tot_node = node_pool;
			head = tail = new (tot_node++) state();
		}

		suffix_automaton () {
			init ();
		}

		void extend (int token) {
			tail = tail -> extend (head, token);
		}

	};


	/*	Palindromic tree :
			void palindromic_tree::init () : initializes the tree.
			bool palindromic_tree::extend (int) : extends the string with token.
				returns whether the tree has generated a new node.
				Complexity O (log MAXC).

			odd, even : the root of two trees.
			last : the node representing the last char.
			node::len : the palindromic string length of the node.
	*/

	template <int MAXN = 1E6, int MAXC = 26>
	struct palindromic_tree {

		struct node {
			node *child[MAXC], *fail;
			int len;
			node (int len) : fail (NULL), len (len) {
				memset (child, NULL, sizeof (child));
			}
		} node_pool[MAXN * 2], *tot_node;

		int size, text[MAXN];

		node *odd, *even, *last;

		node *match (node *now) {
			for (; text[size - now -> len - 1] != text[size]; now = now -> fail);
			return now;
		}

		bool extend (int token) {
			text[++size] = token;
			node *now = match (last);
			if (now -> child[token])
				return last = now -> child[token], false;
			last = now -> child[token] = new (tot_node++) node (now -> len + 2);
			if (now == odd) last -> fail = even;
			else {
				now = match (now -> fail);
				last -> fail = now -> child[token];
			}
			return true;
		}

		void init() {
			text[size = 0] = -1;
			tot_node = node_pool;
			last = even = new (tot_node++) node (0); odd = new (tot_node++) node (-1);
			even -> fail = odd;
		}

		palindromic_tree () {
			init ();
		}

	};

}

using namespace string;

int main () {
	return 0;
}

