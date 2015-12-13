#include "../czl_common.h"
#include "gtest/gtest.h"

using namespace std;
using namespace czl_bio;

TEST(StringUtility, erase) {
	string e0 = " \t";
	string e1 = "atcg";
	string e2 = "abcdefghijklmnopqrstuvwxyz";
	string str0 = "  I'm a \thero. ";
	string str1 = "aATCGactggctagctaGGAC";
	string str2 = "ABCDEFGHIJKLMNO{QSLHSADFasfsfa;AAAdsljfbzcABFDG";
	string str3 = "    \t  \t";
	string str4 = "aaaaaaaaaaaaaabbbbbbbbbddddddddddddd";
	StringUtility::erase_all(str0, e0);
	ASSERT_STREQ(str0.c_str(), "I'mahero.");
	StringUtility::erase_all(str1, e1);
	ASSERT_STREQ(str1.c_str(), "ATCGGGAC");
	StringUtility::erase_all(str2, e2);
	ASSERT_STREQ(str2.c_str(), "ABCDEFGHIJKLMNO{QSLHSADF;AAAABFDG");
	StringUtility::erase_all(str3, e0);
	ASSERT_STREQ(str3.c_str(), "");
	StringUtility::erase_all(str4, e2);
	ASSERT_STREQ(str4.c_str(), "");
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
