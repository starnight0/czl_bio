#include "czl_bio_seq.h"
#include <boost/algorithm/string.hpp>

using namespace boost;

namespace czl_bio {

void seq_init_codon_table()
{
	Int i;

	for (i=0; i<N_CN0; i++) {
		MAP_CODON_TO_ID[CODON[i]] = i;
		MAP_CODON_TO_AA1_GENERAL[CODON[i]] = CODON_AA1_GENERAL[i];
	}
	for (i=0; i<N_NT; i++) {
		MAP_NT1_TO_ID[NT1[i]] = i;
	}
	for (i=0; i<N_AA; i++) {
		MAP_AA1_TO_ID[AA1[i]] = i;
	}
}

char seq_to_complement1(char nt1)
{
	char nt2;
	switch(nt1) {
		case 'T': nt2='A'; break;
		case 't': nt2='a'; break;
		case 'C': nt2='G'; break;
		case 'c': nt2='g'; break;
		case 'A': nt2='T'; break;
		case 'a': nt2='t'; break;
		case 'G': nt2='C'; break;
		case 'g': nt2='c'; break;
		case 'U': nt2='A'; break;
		case 'u': nt2='a'; break;
		default: 
			stringstream ss;
		   	ss << nt1 << " is invalid character for nucleic acid";
			Msg::warn(ss);
			nt2 = '?';
	}
	return nt2;
}

Int seq_to_complement(string & seq)
{
	for (Int i=0; i<seq.size(); i++) {
		seq[i] = seq_to_complement1(seq[i]);
	}
	return 0;
}

Int seq_to_complement_copy(string & seq, string & tag)
{
	for (Int i=0; i<seq.size(); i++) {
		tag.push_back( seq_to_complement1(seq[i]) );
	}
	return 0;
}

Int seq_to_rev_complement(string & seq)
{
	for (Int i=0; i<seq.size(); i++) {
		seq[i] = seq_to_complement1(seq[i]);
	}
	reverse(seq.begin(), seq.end());
	return 0;
}

Int seq_to_rev_complement_copy(string & seq, string & tag)
{
	for (Int i=seq.size()-1; i>=0; i--) {
		tag.push_back(seq_to_complement1(seq[i]));
	}
	return 0;
}

char seq_codon_to_aa1(string & codon)
{
	string codon_up = boost::to_upper_copy(codon);
	return MAP_CODON_TO_AA1_GENERAL[codon_up];
}

char seq_codon_to_id(string & codon)
{
	string codon_up = boost::to_upper_copy(codon);
	return MAP_CODON_TO_ID[codon];
}

Int seq_translate(string & seq, string & tag)
{
	for (Int i=0; i<seq.size(); i+=3) {
		string s = seq.substr(i, 3);
		boost::to_upper(s);
		tag.push_back(seq_codon_to_aa1(s));
	}
}

Int seq_aa1_to_id(char aa1)
{
	aa1 = toupper(aa1);
	if (MAP_AA1_TO_ID.find(aa1)==MAP_AA1_TO_ID.end()) return -1;
	return MAP_AA1_TO_ID[aa1];
}

Int seq_aa1_to_id(string & aa1)
{
	return seq_aa1_to_id(aa1[0]);
}

Int seq_nt1_to_id(char nt1)
{
	nt1 = toupper(nt1);
	if (MAP_NT1_TO_ID.find(nt1)==MAP_NT1_TO_ID.end()) return -1;
	else return MAP_NT1_TO_ID[nt1];
}

Int seq_nt1_to_id(string & nt1)
{
	return seq_nt1_to_id(nt1[0]);
}

/// read a protein alignment score matrix from file
void seq_load_aa_score_mat(string & file, vector< vector<float> >& out_mat)
// {{{
{
	ifstream fin(file.c_str());
	vector<Int> aa1s;
	Int r=0;
	Int n=0;
	while (!fin.eof()) {
		string line;
		getline(fin, line);
		if (line.empty()) continue;
		if (line[0]=='#') continue;

		vector<string> tab;
		boost::split(tab, line, boost::is_any_of("\t"));
		if (r==0) {
			for (Int i=1; i<tab.size(); i++) {
				if (!tab[i].empty()) {
					Int id = seq_aa1_to_id(tab[i]);
					if (id > n) n = id;
				}
			}
			r=1;
		} else {
			Int id = seq_aa1_to_id(tab[0]);
			if (id > n) n = id;
		}
	}
	n++;
	out_mat.resize(n);
	for (Int i=0; i<n; i++) out_mat[i].resize(n);
	fin.close();

	r=0;
	fin.open(file.c_str());
	while (!fin.eof()) {
		string line;
		getline(fin, line);
		if (line.empty()) continue;
		if (line[0]=='#') continue;

		vector<string> tab;
		boost::split(tab, line, boost::is_any_of("\t"));
		if (r==0) {
			for (Int i=1; i<tab.size(); i++) {
				if (tab[i].empty()) aa1s.push_back(-1);
				else aa1s.push_back(seq_aa1_to_id(tab[i]));
			}
			r=1;
		} else {
			Int id1 = seq_aa1_to_id(tab[0]);
			for (Int i=1; i<tab.size(); i++) {
				if (aa1s[i-1]<0) continue;
				Int id2 = aa1s[i-1];
				out_mat[id1][id2] = atof(tab[i].c_str());
			}
		}
	}
	fin.close();
}
// }}}

void seq_load_aa_score_mat(string & file, Int *out_n, float ***out_mat)
// {{{
{
	ifstream fin(file.c_str());
	vector<Int> aa1s;
	Int r=0;
	Int n=0;
	while (!fin.eof()) {
		string line;
		getline(fin, line);
		if (line.empty()) continue;
		if (line[0]=='#') continue;

		vector<string> tab;
		boost::split(tab, line, boost::is_any_of("\t"));
		if (r==0) {
			for (Int i=1; i<tab.size(); i++) {
				if (!tab[i].empty()) {
					Int id = seq_aa1_to_id(tab[i]);
					if (id > n) n = id;
				}
			}
			r=1;
		} else {
			Int id = seq_aa1_to_id(tab[0]);
			if (id > n) n = id;
		}
	}
	n++;
	(*out_n) = n;

	(*out_mat) = new float*[n];
	for (Int i=0; i<n; i++) (*out_mat)[i] = new float[n];
	fin.close();

	r=0;
	fin.open(file.c_str());
	while (!fin.eof()) {
		string line;
		getline(fin, line);
		if (line.empty()) continue;
		if (line[0]=='#') continue;

		vector<string> tab;
		boost::split(tab, line, boost::is_any_of("\t"));
		if (r==0) {
			for (Int i=1; i<tab.size(); i++) {
				if (tab[i].empty()) aa1s.push_back(-1);
				else aa1s.push_back(seq_aa1_to_id(tab[i]));
			}
			r=1;
		} else {
			Int id1 = seq_aa1_to_id(tab[0]);
			for (Int i=1; i<tab.size(); i++) {
				if (aa1s[i-1]<0) continue;
				Int id2 = aa1s[i-1];
				(*out_mat)[id1][id2] = atof(tab[i].c_str());
			}
		}
	}
	fin.close();
}
// }}}
///

/**
  * @brief    class BioSeq
  */
string BioSeq::get_seq(Int begin, Int len)
{
	string out;
	get_seq(out, begin, len);
	return out;
}

Int BioSeq::get_seq(string & out, Int begin, Int len)
{
	if (begin<0 || begin>seq_m.size()) return 1;
	if (len<0) {
		out = seq_m.substr(begin, seq_m.size()-begin);
	} else {
		if (begin+len > seq_m.size()) len = seq_m.size()-begin;
		out = seq_m.substr(begin, len);
	}
	return 0;
}

//Int BioSeq::get_seq(string & out, Int begin=0, Int len=0)
//{
//	Int l; 
//	char *cs=NULL;
//	out.clear();
//	if (begin + len >= total_len_m ) { len = total_len_m-1 - begin ; }
//	if (begin==0 && len==0) { ///< get the whole sequence
//		if ( seq_m.size() <= total_len_m ) {
//			out = seq_m;
//		} else {
//			fs_m.clear();
//			fs_m.seekg(fs_pos_m);
//			l=0;
//			cs = new char[len_m];
//			while ( l + len_m < total_len_m ) {
//				fs_m.read(len_m, len_m);
//				out += cs;
//				l += len_m;
//			}
//			if ( l < total_len_m ) {
//				fs_m.read(cs, total_len_m-l);
//				out += cs;
//				l = total_len_m;
//			}
//			delete[] cs;
//		}
//	} else {
//		Int end = begin+len;
//		if ( begin >= begin_m ) {
//			if ( end < begin_m + len_m ) {
//				out = seq_m.substr( begin-begin_m, len );
//			} else {
//				Int b = begin-begin_m;
//				l = len_m - b;
//				out = seq_m.substr( b, l );
//				cs = new char[len_m];
//				l = begin_m + len_m;
//				fs_m.clear();
//				fs_m.seekg(fs_pos_m + l );
//				while ( l + len_m < end ) {
//					fs_m.read(l, len_m);
//					out += cs;
//					l += len_m;
//				}
//				if ( l < end ) {
//					fs_m.read(cs, end-l);
//					out += cs;
//					l = end;
//				}
//				delete[] cs;
//			}
//		} else {
//			cs = new char[len_m];
//			/// read into the seq buffer
//			begin_m = begin;
//			fs_m.clear();
//			fs_m.seekg(fs_pos_m + begin_m );
//			fs_m.read(cs, len_m);
//			seq_m = cs;
//			///<
//			if ( len <= len_m ) {
//				out = seq_m.substr(0, len);
//			} else {
//				out = seq_m;
//				l = begin_m + len_m;
//			}
//			out = seq_m.substr( b, l );
//			l = begin_m + len_m;
//			while ( l + len_m < end ) {
//				fs_m.read(l, len_m);
//				out += cs;
//				l += len_m;
//			}
//			if ( l < end ) {
//				fs_m.read(cs, end-l);
//				out += cs;
//				l = end;
//			}
//			delete[] cs;
//		}
//	}
//}

void BioSeq::set(BioSeq & seq)
{
	BioObj::set(seq);
	seq_m = seq.seq_m;
}

void BioSeq::set(const BioSeq & seq)
{
	BioObj::set(seq);
	seq_m = seq.seq_m;
}

BioSeq & BioSeq::operator = (BioSeq & seq)
{
	BioObj::set(seq);
	seq_m = seq.seq_m;
}

BioSeq & BioSeq::operator = (const BioSeq & seq)
{
	BioObj::set(seq);
	seq_m = seq.seq_m;
}

/**
  * @brief    class SeqRegion
  */
// {{{
SeqRegion::SeqRegion(SeqRegion & seq_region)
{
	set(seq_region);
}

SeqRegion::SeqRegion(const SeqRegion & seq_region)
{
	set(seq_region);
}

string SeqRegion::get_seq()
{
	string out_seq = "";
	string tmp;
	BioSeq *seq_pt = dynamic_cast<BioSeq*>(obj_pt_m);
	for (Int i=0; i<ranges_m.size(); i++) {
		Int b = get_begin(i);
		Int e = get_end(i);
		if ( !seq_pt->get_seq(tmp, b, e) ) {
			out_seq += tmp;
		} else {
			return "";
		}
	}
	return out_seq;
}

Int SeqRegion::get_seq(string & out_seq)
{
	out_seq = "";
	string tmp;
	BioSeq *seq_pt = dynamic_cast<BioSeq*>(obj_pt_m);
	for (Int i=0; i<ranges_m.size(); i++) {
		Int b = get_begin(i);
		Int e = get_end(i);
		if ( !seq_pt->get_seq(tmp, b, e) ) {
			out_seq += tmp;
		} else {
			return -1;
		}
	}
	return 0;
}

Int SeqRegion::get_seq_len(Int i)
{
	if (i<0 || i>=ranges_m.size()) { return -1; }
	else {
		return ranges_m[i].second - ranges_m[i].first;
	}
}

Int SeqRegion::get_seq_len()
{
	Int l=0;
	for (Int i=0; i<ranges_m.size(); i++) {
		l += ranges_m[i].second - ranges_m[i].first;
	}
	return l;
}

SeqRegion & SeqRegion::set(SeqRegion & seq_region)
{
	ObjRegion::set(seq_region);
	strand_m = seq_region.strand_m;
	return *this;
}

SeqRegion & SeqRegion::set(const SeqRegion & seq_region)
{
	ObjRegion::set(seq_region);
	strand_m = seq_region.strand_m;
	return *this;
}

SeqRegion & SeqRegion::operator =(SeqRegion & seq_region)
{
	return set(seq_region);
}

SeqRegion & SeqRegion::operator =(const SeqRegion & seq_region)
{
	return set(seq_region);
}
// }}}

};
