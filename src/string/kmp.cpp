#include <bits/stdc++.h>

#ifndef __STRING_KMP
#define __STRING_KMP

namespace string {

	/*	KMP algorithm :
			void kmp::build (const std::string & str) :
				initializes and builds the failure array. Complexity O (n).
			int kmp::find (const std::string &str) :
				finds the first occurence of match in str. Complexity O (n).

			Note : match is cylic when L % (L - 1 - fail[L - 1]) == 0 &&
				L / (L - 1 - fail[L - 1]) > 1, where L = match.size ().
	*/

	template <int MAXN = 1000000>
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
			for (int i = 0, j = -1; i < (int) str.size (); j += str[i] == match[j + 1], ++i) {
				if (j == match.size () - 1) return i - match.size ();
				while (~j && str[i] != match[j + 1]) j = fail[j];
			}
			if (j == match.size () - 1) return str.size () - match.size ();
			return str.size ();
		}

	};

}

#endif


