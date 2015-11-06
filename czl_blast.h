#ifndef CZL_BLAST_H
#define CZL_BLAST_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "common.hh"

using namespace std;

namespace czl_bio {

/// @brief class BlastAlign
// {{{
class BlastAlign {
public:
    int id_t;
    int qid_t;
    int tid_t;
//  string qname_t; 
//  string tname_t;
    float iden_perc_t;
    int align_len_t, mismatch_t, gap_open_t;
    int qbegin_t, qend_t, tbegin_t, tend_t;
    double evalue_t, bit_t;

	int8_t strand_t;
	int qlen_t;
	int tlen_t;
    int match_t;
    int qgap_len_t;
    int tgap_len_t;
    unsigned int flag_t;
	void* data_t; // for custom use

    BlastAlign();
    BlastAlign(const BlastAlign & a);
    BlastAlign & operator =(const BlastAlign& a);

	int cal_after_read();
	float get_qcov_perc();
	float get_tcov_perc();
	float get_qgap_perc();
	float get_tgap_perc();
	int read_int_name(istream & in);
	int read_int_name(string & line);
	int read(istream & in, map<string, int> & name_to_id);
	int read(istream & in, map<string, int> & name_to_id, string & qname, string & tname);
	int read(string & line, map<string, int> & name_to_id);
	int read(string & line, map<string, int> & name_to_id, string & qname, string & tname);
	int write_int_name(string & line);
	int write_int_name(ostream & out);
	int write(string & line, vector<string> & id_to_name);
	int write(ostream & out, vector<string> & id_to_name);

	// read and write in binary mode
	int read_bin(istream & in);
	int write_bin(ostream & out);

    virtual ~BlastAlign();
};
// }}}

};
#endif
