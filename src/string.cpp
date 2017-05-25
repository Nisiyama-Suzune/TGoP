/*	String algorithm :
		Algorithms regarding string.
*/

#include <map>
#include <string>
#include <vector>

namespace string {

	/*	KMP algorithm :
			void kmp::build (const std::string & str) :
				initializes and builds the failure array.
			int kmp::find (const std::string &str) :
				finds the first occurence of match in str.

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

}

#include <cstdio>

using namespace string;

int main () {
	return 0;
}

