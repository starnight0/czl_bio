#ifndef CZL_BAM_H
#define CZL_BAM_H

#include "bam/bam.h"
#include "czl_common.hpp"

using namespace std;

namespace czl_bio {

class Bam {
    public:

    static int split_by_chr_pos(string & bam_file, int is_sorted_by_coord, string & out_dir, int len, int step, string * out_bed_file=NULL);
    static int split_by_chr(string & bam_file, string & out_dir, vector<string> * out_chroms=NULL);
    static int sort_by_name(string & bam_file, string & out_bam_file, string & tmp_dir, int buf_n);
    static int sort_by_name(const char *bam_file, const char *out_bam_file, const char *tmp_dir, int buf_n);
    static int sort_by_name(vector<bam1_t*> & out_bams, int first_pos=0);
    static int sort_by_name(vector<bam1_t*> & bams, vector<bam1_t*> & out_bams, int first_pos=0);
    static int sort_by_name(vector<bam1_t*>::iterator out_bams_begin, vector<bam1_t*>::iterator out_bams_end, int first_pos=0);
    static int sort_by_name(bam1_t** bams, int n__bam, bam1_t** out_bams, int first_pos=0);
    static int sort_by_name(bam1_t** bams, int n__bam, int first_pos=0);
};

int bam_to_fastq(const string & bam_file, const string fastq_files[3]);
int bam_to_fasta_qual(const string & bam_file, const string fasta_files[3], const string qual_files[3]);

};
#endif
