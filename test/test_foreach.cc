#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;
using namespace boost;

int main()
{
	vector<string> a;
	a.push_back("123");
	a.push_back("abc");
	a.push_back("czl");
	BOOST_FOREACH(string s, a) {
		cout << s << endl;
	}
	BOOST_FOREACH(string s, a) {
		s = "111";
	}
	BOOST_FOREACH(string s, a) {
		cout << s << endl;
	}
	BOOST_FOREACH(string & s, a) {
		s = "111";
	}
	BOOST_FOREACH(string s, a) {
		cout << s << endl;
	}
	return 0;
}
