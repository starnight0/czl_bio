#include "czl_blast.h"

using namespace czl_bio;
/// @brief class BlastAlign
// {{{

// class BlastAlign {
// public:
//     int id_t;
//     int qid_t;
//     int tid_t;
//     float iden_perc_t;
//     int align_len_t, mismatch_t, gap_open_t;
//     int qbegin_t, qend_t, tbegin_t, tend_t;
//     double evalue_t, bit_t;
//
//     uint8_t strand_t;
//     float qcov_perc_t;
//     int match_t;
//     int qgap_len_t;
//     int tgap_len_t;
//     uint32_t flag_t;
//     void* data_t; // for custom use
//};

BlastAlign::BlastAlign() {
    id_t = -1;
    qid_t = -1;
    tid_t = -1;
    iden_perc_t = 0;
    align_len_t = 0;
    mismatch_t = 0;
    gap_open_t = 0;
    qbegin_t = 0;
    qend_t = 0;
    tbegin_t = 0;
    tend_t = 0;
    evalue_t = 0;
    bit_t = 0;

	qlen_t = 0;
	tlen_t = 0;
    strand_t = 0;
    match_t = 0;
    qgap_len_t = 0;
    tgap_len_t = 0;
    flag_t = 0;
	data_t = NULL;
}

BlastAlign::BlastAlign(const BlastAlign & a) 
{
    id_t = a.id_t;
    qid_t = a.qid_t;
    tid_t = a.tid_t;
    iden_perc_t = a.iden_perc_t;
    align_len_t = a.align_len_t;
    mismatch_t = a.mismatch_t;
    gap_open_t = a.gap_open_t;
    qbegin_t = a.qbegin_t;
    qend_t = a.qend_t;
    tbegin_t = a.tbegin_t;
    tend_t = a.tend_t;
    evalue_t = a.evalue_t;
    bit_t = a.bit_t;

    qlen_t = a.qlen_t;
    tlen_t = a.tlen_t;
    strand_t = a.strand_t;
    match_t = a.match_t;
    qgap_len_t = a.qgap_len_t;
    tgap_len_t = a.tgap_len_t;
    flag_t = a.flag_t;
    data_t = a.data_t;
}

BlastAlign::~BlastAlign() 
{
}

BlastAlign & BlastAlign::operator =(const BlastAlign& a)
{
    id_t = a.id_t;
    qid_t = a.qid_t;
    tid_t = a.tid_t;
    iden_perc_t = a.iden_perc_t;
    align_len_t = a.align_len_t;
    mismatch_t = a.mismatch_t;
    gap_open_t = a.gap_open_t;
    qbegin_t = a.qbegin_t;
    qend_t = a.qend_t;
    tbegin_t = a.tbegin_t;
    tend_t = a.tend_t;
    evalue_t = a.evalue_t;
    bit_t = a.bit_t;

    qlen_t = a.qlen_t;
    tlen_t = a.tlen_t;
    strand_t = a.strand_t;
    match_t = a.match_t;
    qgap_len_t = a.qgap_len_t;
    tgap_len_t = a.tgap_len_t;
    flag_t = a.flag_t;
    data_t = a.data_t;
    return (*this);
}

int BlastAlign::cal_after_read()
{
	int qlen = qend_t - qbegin_t;
	int tlen = tend_t - tbegin_t;
	match_t = ( iden_perc_t * align_len_t + 0.5) / 100;
	qgap_len_t = align_len_t - qlen;
	tgap_len_t = align_len_t - tlen;
	return 0;
}

float BlastAlign::get_qcov_perc()
{
	return float(qend_t - qbegin_t)*100/qlen_t;
}

float BlastAlign::get_tcov_perc()
{
	return float(tend_t - tbegin_t)*100/tlen_t;
}

float BlastAlign::get_qgap_perc()
{
	return float(qgap_len_t)*100/qlen_t;
}

float BlastAlign::get_tgap_perc()
{
	return float(tgap_len_t)*100/tlen_t;
}

/// @breif read from istream, query or target name should be integer
int BlastAlign::read_int_name(istream & in)
// {{{
{
    if (in.good()) {
		in >> qid_t >> tid_t;
		in >> iden_perc_t >> align_len_t >> mismatch_t >> gap_open_t;
		in >> qbegin_t >> qend_t >> tbegin_t >> qend_t;
		in >> evalue_t >> bit_t;
		int qs, ts;
		if (qbegin_t<qend_t) qs = 1;
		else {
			int t = qbegin_t;
			qbegin_t = qend_t;
			qend_t = t;
			qs = -1;
		}
		if (tbegin_t<tend_t) ts = 1;
		else {
			int t = tbegin_t;
			tbegin_t = tend_t;
			tend_t = t;
			ts = -1;
		}
		strand_t = qs*ts;
		qbegin_t--;
		tbegin_t--;
		cal_after_read();
	} else {
		return -1;
	}
}
// }}}

/// @breif read from string, query or target name should be integer
int BlastAlign::read_int_name(string & line)
// {{{
{
	stringstream ss(line);
	ss >> qid_t >> tid_t;
	ss >> iden_perc_t >> align_len_t >> mismatch_t >> gap_open_t;
	ss >> qbegin_t >> qend_t >> tbegin_t >> qend_t;
	ss >> evalue_t >> bit_t;
	int qs, ts;
    if (qbegin_t<qend_t) qs = 1;
    else {
        int t = qbegin_t;
        qbegin_t = qend_t;
        qend_t = t;
        qs = -1;
    }
    if (tbegin_t<tend_t) ts = 1;
    else {
        int t = tbegin_t;
        tbegin_t = tend_t;
        tend_t = t;
        ts = -1;
    }
    strand_t = qs*ts;
    qbegin_t--;
    tbegin_t--;
	cal_after_read();
	return 0;
}
// }}}

int BlastAlign::read(istream & in, map<string, int> & name_to_id)
// {{{
{
    string line;
    getline(in, line);
    if (in.good()) {
        return read(line, name_to_id);
    } else {
        return -1;
    }
}
// }}}

int BlastAlign::read(istream & in, map<string, int> & name_to_id, string & qname, string & tname)
// {{{
{
    string line;
    getline(in, line);
    if (in.good()) {
        return read(line, name_to_id, qname, tname);
    } else {
        return -1;
    }
}
// }}}

// read a blast from a string 
int BlastAlign::read(string & line, map<string, int> & name_to_id) 
// {{{
{
	string qname, tname;
	return read(line, name_to_id, qname, tname);
}
// }}}

int BlastAlign::read(string & line, map<string, int> & name_to_id, string & qname, string & tname) 
// {{{
{
    int i=0;
    int n=0;
    int j;
    // qname
    n++;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    qname = line.substr(i, j-i);
	if (name_to_id.find(qname)==name_to_id.end()) {
		qid_t = -1;
	} else {
		qid_t = name_to_id[qname];
	}
    // tname
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    tname = line.substr(i, j-i);
	if (name_to_id.find(tname)==name_to_id.end()) {
		tid_t = -1;
	} else {
		tid_t = name_to_id[tname];
	}

    string s;
    // iden perc
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    iden_perc_t = atof(s.c_str());
    // align length
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    align_len_t = atol(s.c_str());
    // mismatch
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    mismatch_t = atol(s.c_str());
    // gap_open 
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    gap_open_t = atol(s.c_str());
    // query begin
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    qbegin_t = atol(s.c_str());
    // query end
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    qend_t = atol(s.c_str());
    int qs = 0;
    if (qbegin_t<qend_t) qs = 1;
    else {
        int t = qbegin_t;
        qbegin_t = qend_t;
        qend_t = t;
        qs = -1;
    }
    // target begin
    int ts = 0;
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    tbegin_t = atol(s.c_str());
    // target end
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    tend_t = atol(s.c_str());
    if (tbegin_t<tend_t) ts = 1;
    else {
        int t = tbegin_t;
        tbegin_t = tend_t;
        tend_t = t;
        ts = -1;
    }
    strand_t = qs*ts;
    qbegin_t--;
    tbegin_t--;

    // E value
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
        return n;
    }
    s = line.substr(i, j-i);
    evalue_t = atof(s.c_str());
    // bit
    n++;
    i = ++j;
    j = line.find("\t", i);
    if (j==string::npos) {
		j = line.size();
    }
    s = line.substr(i, j-i);
    bit_t = atol(s.c_str());
	cal_after_read();

    return 0;
}
// }}}

int BlastAlign::write_int_name(string & line)
// {{{
{
	stringstream out;
	out << qid_t << "\t" << tid_t << "\t";
	out << iden_perc_t << "\t" << align_len_t << "\t" << mismatch_t << "\t" << gap_open_t << "\t";
	out << qbegin_t+1 << "\t" << qend_t << "\t";
	if (strand_t>=0) {
		out << tbegin_t+1 << "\t" << tend_t << "\t";
	} else {
		out << tend_t << "\t" << tbegin_t+1 << "\t";
	}
	out << evalue_t << "\t" << bit_t;
	line = out.str();
	return 0;
}
// }}}

int BlastAlign::write_int_name(ostream & out)
// {{{
{
	if ( !out.good() ) return -1;
	out << qid_t << "\t" << tid_t << "\t";
	out << iden_perc_t << "\t" << align_len_t << "\t" << mismatch_t << "\t" << gap_open_t << "\t";
	out << qbegin_t+1 << "\t" << qend_t << "\t";
	if (strand_t>=0) {
		out << tbegin_t+1 << "\t" << tend_t << "\t";
	} else {
		out << tend_t << "\t" << tbegin_t+1 << "\t";
	}
	out << evalue_t << "\t" << bit_t;
	if ( !out.good() ) return -2;
	return 0;
}
// }}}

int BlastAlign::write(string & line, vector<string> & id_to_name) 
// {{{
{
	if ( qid_t<0 || qid_t >= id_to_name.size() ) return 1;
	if ( tid_t<0 || tid_t >= id_to_name.size() ) return 2;
	stringstream out;
	out << id_to_name[qid_t] << "\t" << id_to_name[tid_t] << "\t";
	out << iden_perc_t << "\t" << align_len_t << "\t" << mismatch_t << "\t" << gap_open_t << "\t";
	out << qbegin_t+1 << "\t" << qend_t << "\t";
	if (strand_t>=0) {
		out << tbegin_t+1 << "\t" << tend_t << "\t";
	} else {
		out << tend_t << "\t" << tbegin_t+1 << "\t";
	}
	out << evalue_t << "\t" << bit_t;
	line = out.str();
	return 0;
}
// }}}

int BlastAlign::write(ostream & out, vector<string> & id_to_name) 
// {{{
{
	if ( qid_t<0 || qid_t >= id_to_name.size() ) return 1;
	if ( tid_t<0 || tid_t >= id_to_name.size() ) return 2;
	if ( !out.good() ) return -1;
	out << id_to_name[qid_t] << "\t" << id_to_name[tid_t] << "\t";
	out << iden_perc_t << "\t" << align_len_t << "\t" << mismatch_t << "\t" << gap_open_t << "\t";
	out << qbegin_t+1 << "\t" << qend_t << "\t";
	if (strand_t>=0) {
		out << tbegin_t+1 << "\t" << tend_t << "\t";
	} else {
		out << tend_t << "\t" << tbegin_t+1 << "\t";
	}
	out << evalue_t << "\t" << bit_t;
	if ( !out.good() ) return -2;
	return 0;
}
// }}}

int BlastAlign::read_bin(istream & in)
// {{{
{
	if (in.good()) {
		in.read((char*)&id_t, sizeof(int));
		if (in.gcount()<sizeof(int)) return 1;
		in.read((char*)&qid_t, sizeof(int));
		in.read((char*)&qlen_t, sizeof(int));
		in.read((char*)&tid_t, sizeof(int));
		in.read((char*)&tlen_t, sizeof(int));
		in.read((char*)&iden_perc_t, sizeof(float));
		in.read((char*)&align_len_t, sizeof(int));
		in.read((char*)&mismatch_t, sizeof(int));
		in.read((char*)&gap_open_t, sizeof(int));
		in.read((char*)&qbegin_t, sizeof(int));
		in.read((char*)&qend_t, sizeof(int));
		in.read((char*)&tbegin_t, sizeof(int));
		in.read((char*)&tend_t, sizeof(int));
		in.read((char*)&strand_t, sizeof(int8_t));
		in.read((char*)&evalue_t, sizeof(float));
		in.read((char*)&bit_t, sizeof(float));
		in.read((char*)&flag_t, sizeof(unsigned int));
		if (in.gcount()<sizeof(unsigned int)) return 1;
		cal_after_read();
	} else {
		return -1;
	}
}
// }}}

int BlastAlign::write_bin(ostream & out)
// {{{
{
	if (out.good()) {
		out.write((char*)&id_t, sizeof(int));
		out.write((char*)&qid_t, sizeof(int));
		out.write((char*)&qlen_t, sizeof(int));
		out.write((char*)&tid_t, sizeof(int));
		out.write((char*)&tlen_t, sizeof(int));
		out.write((char*)&iden_perc_t, sizeof(float));
		out.write((char*)&align_len_t, sizeof(int));
		out.write((char*)&mismatch_t, sizeof(int));
		out.write((char*)&gap_open_t, sizeof(int));
		out.write((char*)&qbegin_t, sizeof(int));
		out.write((char*)&qend_t, sizeof(int));
		out.write((char*)&tbegin_t, sizeof(int));
		out.write((char*)&tend_t, sizeof(int));
		out.write((char*)&strand_t, sizeof(int8_t));
		out.write((char*)&evalue_t, sizeof(float));
		out.write((char*)&bit_t, sizeof(float));
		out.write((char*)&flag_t, sizeof(unsigned int));
	} else {
		return -1;
	}
}
// }}}

// }}}

