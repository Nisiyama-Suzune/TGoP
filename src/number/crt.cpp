#include "basic.cpp"

#ifndef __NUMBER_CRT
#define __NUMBER_CRT

namespace number {

	/*	Chinese remainder theroem :
			bool crt::solve (const std::vector <std::pair<long long, long long> > &input,
			                std::pair<long long, long long> &output) :
				solves for an integer set x = output.first + k * output.second
					that satisfies x % input[i].second = input[i].first.
				Returns whether a solution exists.
	*/

	struct crt {

		long long fix (const long long &a, const long long &b) {
			return (a % b + b) % b;
		}

		bool solve (const std::vector <std::pair <long long, long long> > &input,
		            std::pair <long long, long long> &output) {
			output = std::make_pair (1, 1);
			for (int i = 0; i < (int) input.size (); ++i) {
				long long number, useless;
				euclid (output.second, input[i].second, number, useless);
				long long divisor = gcd (output.second, input[i].second);
				if ((input[i].first - output.first) % divisor) {
					return false;
				}
				number *= (input[i].first - output.first) / divisor;
				number = fix (number, input[i].second);
				output.first += output.second * number;
				output.second *= input[i].second / divisor;
				output.first = fix (output.first, output.second);
			}
			return true;
		}

	};

}

#endif

