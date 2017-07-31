#include <bits/stdc++.h>

int main () {
	/*	C++11 supports regular expressions, see below for an example.	*/

	char str[] = "The thirty-three thieves thought that they thrilled the throne throughout Thursday.";
	std::regex pattern ("(th|Th)[\\w]*", std::regex_constants::optimize | std::regex_constants::ECMAScript);
	std::match_results <char *> match;
	std::regex_constants::match_flag_type flag = std::regex_constants::match_default;
	int begin = 0, end = strlen (str);
	while (std::regex_search (str + begin, str + end, match, pattern, flag)) {
		std::cout << match[0] << " " << match[1] << std::endl;
		begin += match.position (0) + 1;
		flag |= std::regex_constants::match_prev_avail;
	}
}

