#include "../czl_bio_align.h"
#include "../czl_io.h"
#include "gtest/gtest.h"

using namespace std;
using namespace czl_bio;

typedef IOBLAST<int, int> IOBLAST1;

TEST(ClassAlignTest, First) {
    cout << "Test czl_bio_align\n";
    string in_dir("input");
    string blast_m6 = in_dir + "/in.blast.m6";
    ifstream fin;
    Align<int, int> align;
    map<string, int> name_to_id;
    vector<string> name_v;
    if (is_gz(blast_m6)) {
    } else {
        fin.open(blast_m6.c_str());
        int r=0;
        r = IOBLAST1::get_one(fin, align, name_to_id, name_v);
        // - first line is:
        // - 4740991 9 98.40  188  3  0  15  202  971612  971799  4e-89  331
        // -   70TA19AG68TC28
        // - query align length = 188
        // - target align length = 188
        Range<int> trange = align.get_range_ref(0);
        Range<int> qrange = align.get_range_ref(1);
        AlignStat tstat = align.get_stat_ref(0);
        AlignStat qstat = align.get_stat_ref(1);
        ASSERT_STREQ(name_v[0].c_str(), "4740991");
        ASSERT_STREQ(name_v[1].c_str(), "9");
        ASSERT_EQ(align.get_seq_id(0), 1);
        ASSERT_EQ(align.get_seq_id(1), 0);
        ASSERT_EQ(qrange.get_begin(), 14);
        ASSERT_EQ(qrange.get_end(), 202);
        ASSERT_EQ(qrange.get_strand(), '+');
        ASSERT_EQ(trange.get_begin(), 971611);
        ASSERT_EQ(trange.get_end(), 971799);
        ASSERT_EQ(trange.get_strand(), '+');
        ASSERT_EQ(tstat.get_score(), 4e-89);
        ASSERT_EQ(qstat.get_score(), 331.0);
        ASSERT_EQ(qstat.get_len(), 188);
        ASSERT_EQ(qstat.get_match(), 185);
        ASSERT_EQ(qstat.get_mismatch(), 3);
        ASSERT_EQ(qstat.get_gap_n(), 0);
    }
}

TEST(ClassAlignTest, LAST) {
    cout << "Test czl_bio_align\n";
    string in_dir("input");
    string blast_m6 = in_dir + "/in.blast.m6";
    ifstream fin;
    Align<int, int> align;
    map<string, int> name_to_id;
    vector<string> name_v;
    if (is_gz(blast_m6)) {
    } else {
        int n=0;
		int r=0;
        fin.open(blast_m6.c_str());
		Align<int, int> align;
        while ( !(r=IOBLAST1::get_one(fin, align, name_to_id, name_v)) ) {
            n++;
        }
        ASSERT_EQ(n, 1000);
        // - last line:
        // - 6456324  7  92.98  57  4  0  1  57  25119176  25119232  1e-14  84.2
        // -   7CT9GC28TC4GC5
        Range<int> & trange = align.get_range_ref(0);
        Range<int> & qrange = align.get_range_ref(1);
        AlignStat & tstat = align.get_stat_ref(0);
        AlignStat & qstat = align.get_stat_ref(1);
		ASSERT_STREQ(name_v[align.get_seq_id(1)].c_str(), "6456324");
		ASSERT_STREQ(name_v[align.get_seq_id(0)].c_str(), "7");
		ASSERT_EQ(qrange.get_begin(), 1-1);
		ASSERT_EQ(qrange.get_end(), 57);
		ASSERT_EQ(qrange.get_strand(), '+');
		ASSERT_EQ(trange.get_begin(), 25119176-1);
		ASSERT_EQ(trange.get_end(), 25119232);
		ASSERT_EQ(trange.get_strand(), '+');
		ASSERT_EQ(tstat.get_score(), 1e-14);
		ASSERT_EQ(qstat.get_score(), 84.2);
		ASSERT_EQ(qstat.get_len(), 57);
		ASSERT_EQ(qstat.get_match(), 57-4);
		ASSERT_EQ(qstat.get_mismatch(), 4);
		ASSERT_EQ(qstat.get_gap_n(), 0);
	}
}

TEST(Class_IOBLAT_Test, First) {
    cout << "Test IOBLAT\n";
    string in_dir("input");
    string blat_psl = in_dir + "/in.blat.psl";
    ifstream fin;
    Align<int, int> align;
    map<string, int> name_to_id;
    vector<string> name_v;
    vector<int> len_v;
    if (is_gz(blat_psl)) {
    } else {
        fin.open(blat_psl.c_str());
        int r=0;
        r = IOBLAT<int,int>::get_one(fin, align, name_to_id, name_v, len_v);
        // - first line is:
        // - 118  1  0  0  0  0  0  0  -  c8_g1_i1  119  0  119  
        //     gi|688619196|ref|XM_003201735.3|  3973  2161  2280  
        //     1  119,  0,  2161,
        Range<int> trange = align.get_range_ref(0);
        Range<int> qrange = align.get_range_ref(1);
        AlignStat tstat = align.get_stat_ref(0);
        AlignStat qstat = align.get_stat_ref(1);
        ASSERT_STREQ(name_v[align.get_seq_id(0)].c_str(), "gi|688619196|ref|XM_003201735.3|");
        ASSERT_STREQ(name_v[align.get_seq_id(1)].c_str(), "c8_g1_i1");
        ASSERT_EQ(qrange.get_begin(), 0);
        ASSERT_EQ(qrange.get_end(), 119);
        ASSERT_EQ(qrange.get_strand(), '-');
        ASSERT_EQ(trange.get_begin(), 2161);
        ASSERT_EQ(trange.get_end(), 2280);
        ASSERT_EQ(trange.get_strand(), '+');
        ASSERT_EQ(tstat.get_score(), 0);
        ASSERT_EQ(qstat.get_score(), 0);
        ASSERT_EQ(qstat.get_len(), 119);
        ASSERT_EQ(qstat.get_match(), 118);
        ASSERT_EQ(qstat.get_mismatch(), 1);
        ASSERT_EQ(qstat.get_gap_n(), 0);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


