#include "czl_io.hpp"

using namespace std;

namespace czl_bio {

//class BioSeq;
//class BioFileSeq;

/// class Fasta
// {{{
/**
  * @brief   get all sequence from a fasta file
  * @details
  *   get all sequence from istream, directly into memory
  *
  * @param[in]   fin      The istream of the fasta input file
  * @param[out]  seqs     The vector of BioSeq contain all output sequences
  * @return    If failed, no sequences in 'seqs'
  * @warnning  All lines start with '#' will be ignored in fasta file
  */ 
    /*
void Fasta::get_all_seq(ifstream & fin, vector<BioSeq> & seqs)
{
    string line;
    Int n=0;
    BioSeq seq;

    while (!fin.eof()) {
        getline(fin, line);
        boost::trim(line);
        if (line.size()==0) continue;
        if (line[0]=='#') continue;

        if (line[0]=='>') {
            string id = line.substr(1);
            boost::trim(id);
            seq.set_id(id);
            seq.set_seq("");

            seqs.push_back(seq);
            n++;
        } else {
            boost::erase_all(line, " ");
            boost::erase_all(line, "\t");
            seqs[n-1].app_seq(line);
        }
    }
}
*/

/**
  * @brief   get all sequence from a fasta file
  * @details
  *   get all sequence from istream, and put them into temporary files, 
  *   each sequence in one file
  *
  * @param[in]   fin      The istream of the fasta input file
  * @param[out]  seqs     The vector of BioFileSeq contain all output sequences
  * @param[in]   out_dir  The directory including all temporary sequence files in BioFileSeq
  * @return    If failed, no sequences in 'seqs'
  * @warnning  All lines start with '#' will be ignored in fasta file
  */ 
/*
void Fasta::get_all_seq(ifstream & fin, vector<BioFileSeq> & seqs, string & out_dir)
{
    string line;
    Int n=0;
    BioFileSeq seq;

    while (!fin.eof()) {
        getline(fin, line);
        boost::trim(line);
        if (line.size()==0) continue;
        if (line[0]=='#') continue;

        if (line[0]=='>') {
            string id = line.substr(1);
            boost::trim(id);
            seq.set_id(id);
            seq.file_m = out_dir+"/"+id;
            seq.total_len_m = 0;

            seqs.push_back(seq);
            if ( n>0 && seqs[n-1].fs_m.is_open() ) { 
                seqs[n-1].fs_m.close();
            }
            seqs[n].fs_m.open(seq.file_m.c_str(), ios::out|ios::binary);
            n++;
        } else {
            boost::erase_all(line, " ");
            boost::erase_all(line, "\t");
            seqs[n-1].fs_m.write( line.c_str(), line.size() );
            seqs[n-1].total_len_m += line.size();
        }
    }
    if ( n>0 && seqs[n-1].fs_m.is_open() ) { 
        seqs[n-1].fs_m.close();
    }
//      Int l=0;
//      while (!fin.eof()) {
//          getline(fin, line);
//          boost::trim(line);
//          if (line.size()==0) continue;
//          if (line[0]=='#') continue;

//          if (line[0]=='>') {
//              if (n>0) { seqs.push_back(seq); }
//              string id = line.substr(1);
//              boost::trim(id);
//              seq.set_id(id);
//              seq.seq_m   = "";
//              seq.begin_m = 0;
//              seq.len_m   = mem_len;
//              fs_m = fs;
//              seq.fs_pos_m = fs_m.tellg();
//              seqs.push_back(seq);
//              l=0;
//              n++;
//          } else {
//              boost::erase_all(line, " ");
//              boost::erase_all(line, "\t");
//              if ( l + line.size() < mem_len ) {
//                  seqs[n-1].seq_m += line;
//                  l+=line.size();
//              } else if ( l < mem_len ) {
//                  seqs[n-1].seq_m += line.substr(0, mem_len-l);
//                  l = mem_len;
//              }
//              seqs[n-1].fs_m.write( line.c_str(), line.size() );
//          }
//      }
}
*/

/*
ifstream & Fasta::get_a_seq(ifstream & fin, BioSeq & seq)
{
    char c;
    short si=0;
    string line;
    seq.set_seq("");
    while (!fin.eof()) {
        fin >> c;
        if (c=='>') {
            if (si==0) {
                getline(fin, line);
                boost::trim(line);
                seq.set_id(line);
                si=1;
            } else {
                fin.unget();
                break;
            }
        } else if (!isspace(c)) {
            fin.unget();
            getline(fin, line);
            boost::trim(line);
            boost::erase_all(line, " ");
            boost::erase_all(line, "\t");
            seq.app_seq(line);
        }
    }
    return fin;
}
*/

/*
ifstream & Fasta::get_a_seq(ifstream & fin, BioFileSeq & seq, string & file)
{
    char c;
    short si=0;
    string line;
    seq.file_m = file;
    seq.fs_m.open(file.c_str(), ios::binary);
    if ( !seq.fs_m.is_open() ) { return fin; }

    boost::regex e("\\s");
    Int l;
    while (!fin.eof()) {
        fin >> c;
        if (c=='>') {
            if (si==0) {
                getline(fin, line);
                boost::trim(line);
                seq.set_id(line);
                l=0;
                si=1;
            } else {
                fin.unget();
                break;
            }
        } else if (!isspace(c)) {
            fin.unget();
            getline(fin, line);
            boost::trim(line);
            boost::erase_all_regex(line, e);
            seq.fs_m.write( line.c_str(), line.size() );
            seq.total_len_m += line.size();
        }
    }
    return fin;
}
*/

int Fasta::get_a_seq(istream & fin, string & id, string & seq)
{
    char c;
    short si=0;
	int r = 0;
    string line;
    seq.clear();
    while (!fin.eof()) {
		c='\0';
        fin >> c;
        if (c=='>') {
            if (si==0) {
                getline(fin, line);
                StringUtility::trim(line, " \t\n\r");
				if ( line.empty() ) {
					id.clear();
					r = 1; // seq name is empty
					break;
				}
                id = line;
                si=1;
            } else {
                fin.unget();
                break;
            }
        } else if ( !(c=='\0' || isspace(c)) ) {
            line.clear();
            getline(fin, line);
            StringUtility::trim(line, " \t\r\n");
            StringUtility::erase_all(line, " \t\r\n");
            seq+=c;
            seq+=line;
        }
    }
	if ( seq.empty() ) {
		r = 2;
	}
    return r;
}

/*
ofstream & Fasta::put_all_seq(ofstream & fout, vector<BioSeq> & seqs)
{
    for (Int i=0; i<seqs.size(); i++) {
        put_a_seq(fout, seqs[i]);
    }
    return fout;
}
*/

int Fasta::put_a_seq(ostream & fout, const string & id, const string & seq)
{
    fout << ">" << id << "\n" << seq << "\n";
	if ( fout.fail() ) return 1;
	else return 0;
}

/*
ofstream & Fasta::put_a_seq(ofstream & fout, BioSeq & seq)
{
    fout << ">" << seq.get_id() << "\n";
    string s;
    seq.get_seq(s);
    fout << s << "\n";
    return fout;
}
*/

int Fasta::split_by_id(string const & fasta_file, string const & out_prefix, vector<string> & out_ids)
{
	int r = 0;
    ifstream fin(fasta_file.c_str());
    if ( !fin.is_open() ) {
        return -1;
    }
    string line, out_file;
    ofstream fout;
    string id, seq;
    while (!fin.eof()) {
		if ( get_a_seq(fin, id, seq) ) break;
		out_file = out_prefix + id + ".fa";
		fout.open(out_file.c_str());
		if ( fout.fail() ) {
			msg.error( string("Fail to open File ") + out_file + " " + CZL_DBG_INFO);
			r = -2;
		}
		put_a_seq(fout, id, seq);
		out_ids.push_back(id);
		fout.close();
    }
    fin.close();
	if ( r!=0 ) {
		for ( size_t i=0; i<=out_ids.size(); i++ ) {
			out_file = out_prefix + out_ids[i] + ".fa";
			if ( File::exists_file(out_file) ) {
				remove(out_file.c_str());
			}
		}
	}
	return 0;
}

/*
 * @brief split fasta file into fasta files, each with only one sequence of length 'length'
 * @param  fasta_file  input fasta file
 * @param  out_prefix  output prefix
 * @param  length      each sequence length
 * @param  step        begin position step
 * @param  out_bed_file  output bed files for splited position
 * @return 0  if success, other if failed
 */
int Fasta::split_by_id_pos(string const & fasta_file, string const & out_prefix, int length, int step, string * out_bed_file)
{
    ifstream fin(fasta_file.c_str());
    if ( !fin.is_open() ) {
        cerr << "Can't open FILE " << fasta_file << endl;
        return 1;
    }
    string id, seq;
    string file;
    ofstream fout_bed;
    if (out_bed_file!=NULL) {
        fout_bed.open(out_bed_file->c_str());
        if ( !fout_bed.is_open() ) {
            cerr << "Can't opne FILE " << (*out_bed_file) << endl;
            return 2;
        }
    }
    while (!fin.eof()) {
        get_a_seq(fin, id, seq);
        int begin = 0;
        int end;
        while (begin < seq.size()) {
            end = begin+length;
            if (end > seq.size()) end = seq.size();
            stringstream ss, name;
            name << id << "_" << begin << "_" << begin+length << ".fa";
            ss << out_prefix << name.str();
            ofstream fout(ss.str().c_str());
            if ( !fout.is_open() ) {
                cerr << "Can't opne FILE " << ss.str() << endl;
                return 3;
            }
            put_a_seq(fout, name.str(), seq.substr(begin, end-begin));
            fout.close();
            fout_bed << id << "\t" << begin << "\t" << begin+length << "\n";
            begin+=step;
        }
        id.clear();
        seq.clear();
    }
    fin.close();
    if (out_bed_file!=NULL) {
        fout_bed.close();
    }

    return 0;
}

int Fasta::get_seq_num(string const & in_file, size_t & n)
{
	int r = 0;
    ifstream fin(in_file.c_str());
    if ( !fin.is_open() ) {
        msg.error(string("Can't open INFILE ") + in_file + CZL_DBG_INFO);
        return 1;
    }
	n = 0;
	string id, seq;
    while (!fin.eof()) {
        if ( get_a_seq(fin, id, seq) ) { r = 2; break; }
		n++;
	}
	fin.close();
	return r;
}

int Fasta::get_seq_bp(string const & in_file, size_t & n)
{
	int r = 0;
    ifstream fin(in_file.c_str());
    if ( !fin.is_open() ) {
        msg.error(string("Can't open INFILE ") + in_file + CZL_DBG_INFO);
        return 1;
    }
	n = 0;
	string id, seq;
    while (!fin.eof()) {
        if ( get_a_seq(fin, id, seq) ) { r = 2; break; }
		n+=seq.size();
	}
	fin.close();
	return r;
}

int Fasta::split_by_seq_num(string const & in_file, string const & out_prefix, size_t n)
{
	int r=0;
    ifstream fin(in_file.c_str());
    if ( !fin.is_open() ) {
        msg.error(string("Can't open INFILE ") + in_file + CZL_DBG_INFO);
        return 1;
    }
    string id, seq;
    string out_file;
	size_t m=0, k=0;
    ofstream fout;
    while (!fin.eof()) {
        if ( get_a_seq(fin, id, seq) ) { r = 3; break; }
		if ( ! fout.is_open() ) {
			out_file = out_prefix + itos(m) + ".fa";
			fout.open(out_file.c_str());
			if ( fout.fail() ) {
				msg.error(string("Can't open OUTFILE ") + out_file + CZL_DBG_INFO);
				return 2;
			}
		}
		if ( put_a_seq(fout, id, seq) ) { r=4; break; }
		k++;
		if ( k == n ) {
			fout.close();
			m++;
			k=0;
		}
	}
	fin.close();
	if ( fout.is_open() ) fout.close();
	if ( r!=0 ) {
		for ( size_t i=0; i<=m; i++ ) {
			out_file = out_prefix + itos(m) + ".fa";
			if ( File::exists_file(out_file) ) {
				remove(out_file.c_str());
			}
		}
	}
	return r;
}

int Fasta::split_by_file_num_eq_seq_num(string const & in_file, string const & out_prefix, size_t n)
{
	int r=0;
	size_t seq_n = 0;
	if ( r = get_seq_num(in_file, seq_n) ) return r;

	size_t per_n;
	per_n = seq_n/n;

    ifstream fin(in_file.c_str());
    if ( !fin.is_open() ) {
        msg.error(string("Can't open INFILE ") + in_file + CZL_DBG_INFO);
        return 1;
    }
    string id, seq;
    string out_file;
	size_t m=0, k=0;
    ofstream fout;
    while (!fin.eof()) {
        if ( get_a_seq(fin, id, seq) ) { r = 3; break; }
		if ( ! fout.is_open() ) {
			out_file = out_prefix + itos(m) + ".fa";
			fout.open(out_file.c_str());
			if ( fout.fail() ) {
				msg.error(string("Can't open OUTFILE ") + out_file + CZL_DBG_INFO);
				return 2;
			}
		}
		if ( put_a_seq(fout, id, seq) ) { r=4; break; }
		k++;
		if ( m<n-1 && k == per_n ) {
			fout.close();
			m++;
			k=0;
		}
	}
	fin.close();
	if ( fout.is_open() ) fout.close();
	if ( r!=0 ) {
		for ( size_t i=0; i<=m; i++ ) {
			out_file = out_prefix + itos(m) + ".fa";
			if ( File::exists_file(out_file) ) {
				remove(out_file.c_str());
			}
		}
	}

	return r;
}

int Fasta::split_by_file_num_eq_bp(string const & in_file, string const & out_prefix, size_t n)
{
	int r=0;
	size_t bp = 0;
	if ( r = get_seq_bp(in_file, bp) ) return r;

	size_t per_bp;
    if ( bp%n==0 ) {
		per_bp = bp/n;
	} else {
		per_bp = bp/(n-1);
	}

	if ( per_bp==0 ) per_bp = 100;

	r = split_by_bp(in_file, out_prefix, per_bp);

	return r;
}


int Fasta::split_by_bp(string const & in_file, string const & out_prefix, size_t n)
{
	int r=0;
    ifstream fin(in_file.c_str());
    if ( fin.fail() ) {
        msg.error(string("Can't open INFILE ") + in_file + CZL_DBG_INFO);
        return 1;
    }
    string id, seq;
    string out_file;
	size_t m=0, k=0;
    ofstream fout;
    while (!fin.eof()) {
        if ( get_a_seq(fin, id, seq) ) { r = 3; break; }
		if ( !fout.is_open() ) {
			out_file = out_prefix + itos(m) + ".fa";
			fout.open(out_file.c_str());
			if ( fout.fail() ) {
				msg.error(string("Can't open OUTFILE ") + out_file + CZL_DBG_INFO);
				r=2;
				break;
			}
		}
		if ( put_a_seq(fout, id, seq) ) { r=4; break; }
		k += seq.size();
		if ( k >= n ) {
			fout.close();
			m++;
			k=0;
		}
	}
	fin.close();
	if ( fout.is_open() ) fout.close();
	if ( r!=0 ) {
		for ( size_t i=0; i<=m; i++ ) {
			out_file = out_prefix + itos(m) + ".fa";
			if ( File::exists_file(out_file) ) {
				remove(out_file.c_str());
			}
		}
	}
	return r;
}
// }}}

// class Fastq
// {{{
int Fastq::get_a_seq(istream & s, string & id, string & seq,
        vector<int8_t> & qual)
{
    string line;
    int r = -1;
    while ( !s.eof() ) {
        getline(s, line);
        if ( line.empty() ) continue;
        if (line[0]=='@') {
            id = line.substr(1);
            line.clear();
            getline(s, seq);
            StringUtility::trim(seq, " \t\n\r");
            if (seq.empty()) { r = 1; break; }
            getline(s, line);
            if (line[0]!='+') { r = 2; break; }
            line.clear();
            getline(s, line);
            if (line.size()!=seq.size()) { 
                cerr<< "E: length of QUAL is not the same as SEQ. "
                    << CZL_DBG_INFO << endl;
                r = 3;
                break; 
            }
            if ( Fastq::str_to_qual(line, qual) ) {
                r=4;
            }
            line.clear();
            r = 0;
            break;
        }
        line.clear();
    }
    return r;
}

int Fastq::get_a_seq(istream & s, string & id, string & seq, string & qual)
{
    string line;
    int r = -1;
    while ( !s.eof() ) {
        getline(s, line);
        if ( line.empty() ) continue;
        if (line[0]=='@') {
            id = line.substr(1);
            line.clear();
            getline(s, seq);
            StringUtility::trim(seq, " \t\n\r");
            if (seq.empty()) { r = 1; break; }
            getline(s, line);
            if (line[0]!='+') { r = 2; break; }
            qual.clear();
            getline(s, qual);
            if (qual.size()!=seq.size()) { 
                cerr<< "E: length of QUAL is not the same as SEQ. "
                    << CZL_DBG_INFO << endl;
                r = 3;
                break; 
            }
            r = 0;
            break;
        }
        line.clear();
    }
    return r;
}

int Fastq::put_a_seq(ostream & s, string & id, string & seq,
        vector<int8_t> & qual)
{
    s << '@' << id << '\n' << seq << "\n+\n";
    string qual_s;
    for (int i=0; i<qual.size(); i++) {
        s << static_cast<char>(qual[i]+33);
    }
    if (s.fail()) return 1;
    else return 0;
}

int Fastq::put_a_seq(ostream & s, string & id, string & seq, string & qual)
{
    s << '@' << id << '\n' << seq << "\n+\n" << qual;
    if (s.fail()) return 1;
    else return 0;
}

int Fastq::str_to_qual(string const & s, vector<int8_t> & qual)
{
    qual.resize(s.size());
    for (int i=0; i<s.size(); i++) {
        qual[i] = s[i]-33;
        // if quality < 0
        if (qual[i]<0) {
            cerr<< "E: QUAL is below Zero at position " << i << " "
                << CZL_DBG_INFO << endl;
            return 1;
        }
    }
    return 0;
}

int Fastq::qual_to_str(vector<int8_t> const & qual, string & s)
{
    s.clear();
    for (int i=0; i<qual.size(); i++) {
        s.push_back(static_cast<char>(qual[i]+33));
    }
    return 0;
}
// }}}

/// class VCF
//// {{{
//string VCF::get_version(ifstream & fin)
//{
//    string line;
//    vector<string> sv;
//    Int i;
//
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.size()==0) continue;
//        if (line[0]=='#' && line[1]=='#') {
//            i = line.find_first_of("=", 2);
//            string name;
//            string value;
//            if ( i != std::string::npos ) {
//                name = line.substr(2, i-2);
//                boost::trim(name);
//                boost::to_lower(name);
//                value = line.substr(i+1);
//                boost::trim(value);
//                if ( name == "fileformat") {
//                    return value;
//                }
//            }
//        } else {
//            break;
//        }
//    }
//
//    return "";
//}
//
//string VCF::get_version(const string & file)
//{
//    ifstream fin(file.c_str());
//    return get_version(fin);
//    fin.close();
//}
//
//void VCF::get_all_variant_one_sample(ifstream & fin, GTP *gtp, vector< Variant* > & variants)
//{
//    string ver = get_version(fin);
//}
//
//void VCF::get_all_variant_one_sample(ifstream & fin, string & version, GTP *gtp, vector< Variant* > & variants)
//// {{{
//{
//    string line;
//    vector<string> sv;
//
//    if (version == "VCFv4.1") {
//        while (!fin.eof()) {
//            getline(fin, line);
//            boost::trim(line);
//            if (line.size()==0) continue;
//            if (line[0]=='#') continue;
//            /// chrom_id  pos  variant_id  ref  alt  qual  filter  info  format  sample1  sample2  ...
//            boost::split(sv, line, boost::is_any_of("\t"));
//            if ( sv[4] == "." ) continue;
//
//            Variant *pv=new Variant;
//            Chrom *chr = gtp->get_chrom_by_symbol(sv[0]);
//            Int b = atoi(sv[1].c_str()) - 1;
//            boost::trim(sv[3]);
//            Int l = sv[3].size();
//            Int e = b+l;
//            pv->get_seq_region_pt()->set_seq_pt(chr);
//            pv->get_seq_region_pt()->ins_range(b, e);
//            pv->ins_mut_seq(sv[3]);
//            vector<string> mut;
//            boost::split(mut, sv[4], boost::is_any_of(","));
//            BOOST_FOREACH(string & s, mut) {
//                boost::trim(s);
//                if (s != "X" && !s.empty() ) { pv->ins_mut_seq(s); }
//            //  gtp->get_variant().push_back(pv);
//            //  gtp->set_variant_by_id_range_mut(chr->get_id(), b, e, s, pv);
//            }
//            float qual   = atof(sv[5].c_str());  ///< ignore
//            pv->set_qual(qual);
//            string filter = sv[6];  ///< ignore
//            /// read the INFO part
//            // {{{
//            vector<string> infos;
//            boost::split(infos, sv[7], boost::is_any_of(";"));
//            BOOST_FOREACH(string & s, infos) {
//                vector<string> s1;
//                boost::split(s1, s, boost::is_any_of("="));
//                string name = s1[0];
//                boost::to_upper(name);
//                if ( name == "DP" ) {
//                //  pv->set_total_read( atol(value.c_str()) );
//                //  vi.set_map_qual( atol(value.c_str()) );
//                } else if ( name == "INDEL" ) {
//                    pv->set_type(VT_INDEL);
//                } else if ( name == "IS" ) { ///< for indel
//                } else if ( name == "DP4" ) {
//                    string value = s1[1];
//                    vector<string> s2;
//                    boost::split(s2, value, boost::is_any_of(","));
//                    for (Int i2=0; i2<s2.size(); ) {
//                        Int f = atol(s2[i2++].c_str());
//                        Int r = atol(s2[i2++].c_str());
//                        pv->ins_read_num(f, r);
//                    }
//                } else if ( name == "MQ" ) {
//                    string value = s1[1];
//                    pv->set_map_qual( atol(value.c_str()) );
//                } else {
//                }
//            }
//            // }}}
//            ///
//            /// read infomation for each sample according to the format
//            vector<string> formats;
//            boost::split(formats, sv[8], boost::is_any_of(":"));
//            for (Int i=0, j=9; j<sv.size(); i++, j++) {
//                vector<string> values;
//                boost::split(values, sv[j], boost::is_any_of(":"));
//                for (Int k=0; k<formats.size(); k++) {
//                    string fm = formats[k];
//                    if ( fm == "GT" ) {
//                        vector<string> s2;
//                        boost::split(s2, values[k], boost::is_any_of("/|"));
//                        pv->set_genotype( atoi(s2[0].c_str()), atoi(s2[1].c_str()) );
//                    } else if ( fm == "DP" ) {
//                    } else if ( fm == "PL" ) {
//                    } else if ( fm == "GQ" ) {
//                    } else {
//                    }
//                }
//                break;
//            }
//            ///
//            variants.push_back(pv);
//        }
//    }
//}
//// }}}
//
//void VCF::get_all_variant_only_seq_region(ifstream & fin, string & version, GTP *gtp, vector<Variant*> & variants)
//// {{{
//{
//    string line;
//    vector<string> sv;
//
//    if (version == "VCFv4.1") {
//        while (!fin.eof()) {
//            getline(fin, line);
//            boost::trim(line);
//            if (line.size()==0) continue;
//            if (line[0]=='#') continue;
//            /// chrom_id  pos  variant_id  ref  alt  qual  filter  info  format  sample1  sample2  ...
//            boost::split(sv, line, boost::is_any_of("\t"));
//            if ( sv[4] == "." ) continue;
//
//            Variant *pv=new Variant;
//            Chrom *chr = gtp->get_chrom_by_symbol(sv[0]);
//            Int b = atoi(sv[1].c_str()) - 1;
//            boost::trim(sv[3]);
//            Int l = sv[3].size();
//            Int e = b+l;
//            pv->get_seq_region_pt()->set_seq_pt(chr);
//            pv->get_seq_region_pt()->ins_range(b, e);
//            pv->ins_mut_seq(sv[3]);
//            vector<string> mut;
//            boost::split(mut, sv[4], boost::is_any_of(","));
//            BOOST_FOREACH(string & s, mut) {
//                boost::trim(s);
//                if (s != "X") { pv->ins_mut_seq(s); }
//            //  gtp->get_variant().push_back(pv);
//            //  gtp->set_variant_by_id_range_mut(chr->get_id(), b, e, s, pv);
//            }
//            float qual   = atof(sv[5].c_str());
//            pv->set_qual(qual);
//        //  string filter = sv[6];  ///< ignore
//            /// read the INFO part
//            // {{{
//            vector<string> infos;
//            boost::split(infos, sv[7], boost::is_any_of(";"));
//            BOOST_FOREACH(string & s, infos) {
//                vector<string> s1;
//                boost::split(s1, s, boost::is_any_of("="));
//                string name = s1[0];
//                boost::to_upper(name);
//                if ( name.compare("DP")==0 ) {
//                //  pv->set_total_read( atol(value.c_str()) );
//                //  vi.set_map_qual( atol(value.c_str()) );
//                } else if ( name.compare("INDEL")==0 ) {
//                    pv->set_type(VT_INDEL);
//            //  } else if ( name.compare("IS")==0 ) { ///< for indel
//            //  } else if ( name.compare("DP4")==0 ) {
//                } else if ( name.compare("MQ")==0 ) {
//                    string value = s1[1];
//                    pv->set_map_qual( atol(value.c_str()) );
//                } else {
//                }
//            }
//            // }}}
//            ///
//            /// read infomation for each sample according to the format
//            vector<string> formats;
//            boost::split(formats, sv[8], boost::is_any_of(":"));
//            for (Int i=0, j=9; j<sv.size(); i++, j++) {
//                VariantInfo *pvi = new VariantInfo;
//                pv->ins_info(pvi);
//                vector<string> values;
//                boost::split(values, sv[j], boost::is_any_of(":"));
//                for (Int k=0; k<formats.size(); k++) {
//                    string fm = formats[k];
//                    if ( fm.compare("GT")==0 ) {
//                        vector<string> s2;
//                        boost::split(s2, values[k], boost::is_any_of("/|"));
//                        pvi->set_genotype( atoi(s2[0].c_str()), atoi(s2[1].c_str()) );
//                    } else if ( fm.compare("DP")==0 ) {
//                        pvi->ins_read_num( atol(values[k].c_str()), 0 );
//                        for (Int i1=1; i1<pv->get_site_num(); i1++) {
//                            pvi->ins_read_num(0,0);
//                        }
//                //  } else if ( fm == "PL" ) {
//                //  } else if ( fm == "GQ" ) {
//                    } else {
//                    }
//                }
//            }
//            ///
//            variants.push_back(pv);
//        }
//    }
//}
//// }}}
//
//void VCF::get_all_variant_only_seq_region(ifstream & fin, GTP *gtp, vector<Variant*> & variants)
//{
//    string ver = get_version(fin);
//    get_all_variant_only_seq_region(fin, ver, gtp, variants);
//}
//
//void VCF::get_all_variant(istream & fin, string & version, GTP *gtp, vector<Variant*> & variants)
//{
//    string line;
//    vector<string> sv;
//
//    if (version == "VCFv4.1") {
//        while (!fin.eof()) {
//            getline(fin, line);
//            boost::trim(line);
//            if (line.size()==0) continue;
//            if (line[0]=='#') continue;
//            /// chrom_id  pos  variant_id  ref  alt  qual  filter  info  format  sample1  sample2  ...
//            boost::split(sv, line, boost::is_any_of("\t"));
//            if ( sv[4] == "." ) continue;
//
//            Variant *pv=new Variant;
//            Chrom *chr = gtp->get_chrom_by_symbol(sv[0]);
//            Int b = atoi(sv[1].c_str()) - 1;
//            boost::trim(sv[3]);
//            Int l = sv[3].size();
//            Int e = b+l;
//            pv->get_seq_region_pt()->set_seq_pt(chr);
//            pv->get_seq_region_pt()->ins_range(b, e);
//            pv->ins_mut_seq(sv[3]);
//            vector<string> mut;
//            boost::split(mut, sv[4], boost::is_any_of(","));
//            BOOST_FOREACH(string & s, mut) {
//                boost::trim(s);
//                if (s != "X") { pv->ins_mut_seq(s); }
//            //  gtp->get_variant().push_back(pv);
//            //  gtp->set_variant_by_id_range_mut(chr->get_id(), b, e, s, pv);
//            }
//            float qual   = atof(sv[5].c_str());
//            pv->set_qual(qual);
//        //  string filter = sv[6];  ///< ignore
//            /// read the INFO part
//            // {{{
//            vector<string> infos;
//            boost::split(infos, sv[7], boost::is_any_of(";"));
//            BOOST_FOREACH(string & s, infos) {
//                vector<string> s1;
//                boost::split(s1, s, boost::is_any_of("="));
//                string name = s1[0];
//                boost::to_upper(name);
//                if ( name.compare("DP")==0 ) {
//                //  pv->set_total_read( atol(value.c_str()) );
//                //  vi.set_map_qual( atol(value.c_str()) );
//                } else if ( name.compare("INDEL")==0 ) {
//                    pv->set_type(VT_INDEL);
//            //  } else if ( name.compare("IS")==0 ) { ///< for indel
//            //  } else if ( name.compare("DP4")==0 ) {
//                } else if ( name.compare("MQ")==0 ) {
//                    string value = s1[1];
//                    pv->set_map_qual( atol(value.c_str()) );
//                } else {
//                }
//            }
//            // }}}
//            ///
//            /// read infomation for each sample according to the format
//            vector<string> formats;
//            boost::split(formats, sv[8], boost::is_any_of(":"));
//            for (Int i=0, j=9; j<sv.size(); i++, j++) {
//                VariantInfo *pvi = new VariantInfo;
//                for (Int i1=0; i1<pv->get_site_num(); i1++) {
//                    pvi->ins_read_num(0,0);
//                }
//                pv->ins_info(pvi);
//                vector<string> values;
//                boost::split(values, sv[j], boost::is_any_of(":"));
//                for (Int k=0; k<formats.size(); k++) {
//                    if (k>=values.size()) break;
//                    string fm = formats[k];
//                    if ( fm.compare("GT")==0 ) {
//                        vector<string> s2;
//                        boost::split(s2, values[k], boost::is_any_of("/|"));
//                        Int t1, t2;
//                        if (s2[0]==".") {
//                            t1 = -1;
//                        } else {
//                            t1 = atoi(s2[0].c_str());
//                        }
//                        if (s2[1]==".") {
//                            t2 = -1;
//                        } else {
//                            t2 = atoi(s2[1].c_str());
//                        }
//                        pvi->set_genotype( t1, t2 );
//                    } else if ( fm == "AD" ) {
//                        vector<string> s2;
//                        boost::split(s2, values[k], boost::is_any_of(",;"));
//                        for (Int i1=0; i1<s2.size(); i1++) {
//                            pvi->set_forw_read_num( i1, atol(s2[i1].c_str()) );
//                        }
//                    } else if ( fm.compare("DP")==0 ) {
//                //  } else if ( fm == "PL" ) {
//                //  } else if ( fm == "GQ" ) {
//                    } else {
//                    }
//                }
//            }
//            ///
//            variants.push_back(pv);
//        }
//    }
//}
//
//int VCF::split_by_chr(const string & vcf_file, const string & out_dir, vector<string> & out_chroms)
//// {{{
//{
//    ifstream fin(vcf_file.c_str());
//    if ( !fin.is_open() ) {
//        stringstream ss;
//        ss << "In VCF::split_by_chrom(vcf_file=" << vcf_file << ", out_dir=" << out_dir << "): ";
//        ss << "Can't opne FILE " << vcf_file << endl;
//        Msg::error(ss);
//    }
//    string line, file;
//    file = out_dir + "/header";
//    ofstream fout(file.c_str());
//    short r = 0;
//    vector<string> tab;
//    string chrom;
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            if (r==0) {
//                fout << line << "\n";
//            }
//        } else {
//            if (r==0) {
//                fout.close();
//                r=1;
//            }
//            boost::split(tab, line, boost::is_any_of("\t"));
//            boost::trim(tab[0]);
//            if (tab[0]!=chrom) {
//                if (fout.is_open()) fout.close();
//                chrom = tab[0];
//                file = out_dir +"/" + chrom;
//                fout.open(file.c_str());
//                out_chroms.push_back(chrom);
//            }
//            fout << line << "\n";
//        }
//    }
//    if (fout.is_open()) fout.close();
//    return 0;
//}
//// }}}
//
//int VCF::split_by_pos(const string & vcf_file, const string & out_dir, Int len, string * out_pos_file)
//{
//    ifstream fin(vcf_file.c_str());
//    if ( !fin.is_open() ) {
//        stringstream ss;
//        ss << "In VCF::split_by_pos(vcf_file=" << vcf_file << ", out_dir=" << out_dir << ", len=" << len << ", out_pos_file=" << out_pos_file << "): ";
//        ss << "Can't opne FILE " << vcf_file << endl;
//        Msg::error(ss);
//    }
//    ofstream fout_pos;
//    if (out_pos_file !=NULL) {
//        fout_pos.open(out_pos_file->c_str());
//        if ( !fout_pos.is_open() ) {
//            stringstream ss;
//            ss << "In VCF::split_by_pos(vcf_file=" << vcf_file << ", out_dir=" << out_dir << ", len=" << len << ", out_pos_file=" << out_pos_file << "): ";
//            ss << "Can't opne FILE " << out_pos_file << endl;
//            Msg::error(ss);
//        }
//    }
//    string line, file;
//    file = out_dir + "/header";
//    ofstream fout(file.c_str());
//    short r = 0;
//    vector<string> tab;
//    string chrom;
//    Int pos1 = len;
//    Int pos;
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            if (r==0) {
//                fout << line << "\n";
//            }
//        } else {
//            if (r==0) {
//                fout.close();
//                r=1;
//            }
//            boost::split(tab, line, boost::is_any_of("\t"));
//            boost::trim(tab[0]);
//            boost::trim(tab[1]);
//            pos = atol(tab[1].c_str());
//            if (tab[0]!=chrom || pos > pos1) {
//                if (fout.is_open()) fout.close();
//                chrom = tab[0];
//                Int pos0 = pos/len * len;
//                pos1 = pos0 + len;
//                stringstream ss1;
//                ss1 << out_dir << "/" << chrom << "_" << pos0 << "_" << pos1;
//                file = ss1.str();
//                fout.open(file.c_str());
//                if (out_pos_file!=NULL) {
//                    fout_pos << chrom << "\t" << pos0 << "\t" << pos1 << "\n";
//                }
//            }
//            fout << line << "\n";
//        }
//    }
//    if (fout.is_open()) fout.close();
//    return 0;
//}
//
//int VCF::split_by_pos(const string & vcf_file, const string & out_dir, Int len, string & out_pos_file)
//{
//    return split_by_pos(vcf_file, out_dir, len, &out_pos_file);
//}
//
//struct TmpA {
//    string chr;    
//    int64_t pos;
//    string data;
//};
///**
// * @brief  split VCF file by id and position
// * @param  vcf_file    VCF file to split
// * @param  out_prefix  out prefix for split files, each out file is named $out_prefix+$chr+'_'+$begin+'_'+$end
// * @param  length      length of each region
// * @param  step        step between nearby regions
// * @return  0 if success, other if fail
// */
//int VCF::split_by_chr_pos(const string & vcf_file, int is_sorted, const string & out_prefix, int64_t length, int64_t step, string * out_bed_file)
//{
//    ifstream fin(vcf_file.c_str());
//    if ( !fin.is_open() ) {
//        cerr << "Can't open FILE " << vcf_file << endl;
//        return 1;
//    }
//    ofstream fout_bed;
//    if (out_bed_file !=NULL) {
//        fout_bed.open(out_bed_file->c_str());
//        if ( !fout_bed.is_open() ) {
//            cerr << "Can't opne FILE " << *out_bed_file << endl;
//            return 2;
//        }
//    }
//    string line, file;
//    file = out_prefix + "header";
//    ofstream fout(file.c_str());
//    short r = 0;
//    vector<string> tab;
//    list<TmpA*> vs;
//    string chr0;
//    int64_t begin, end, len;
//    int64_t pos;
//    int i, j;
//    if (is_sorted) {
//        begin = 0;
//        end = begin+length;
//        while (!fin.eof()) {
//            getline(fin, line);
//            boost::trim(line);
//            if (line.empty()) {
//                continue;
//            } else if (line[0]=='#') {
//                if (r==0) {
//                    fout << line << "\n";
//                }
//            } else {
//                if (r==0) {
//                    fout.close();
//                    r=1;
//                }
//                i = line.find('\t');
//                string chr = line.substr(0, i);
//                j=i+1;
//                i = line.find('\t', j);
//                pos = atol(line.substr(j, i-j).c_str());
//                if (chr!=chr0 || pos >= end) {
//                    if (fout.is_open()) fout.close();
//                    if (chr==chr0) {
//                        begin += step;
//                    } else {
//                        chr0 = chr;
//                        begin = 0;
//                        BOOST_FOREACH(TmpA *a, vs) {
//                            delete a;
//                        }
//                        vs.clear();
//                    }
//                    end = begin+length;
//                    stringstream name, ss1;
//                    name << chr << "_" << begin << "_" << end;
//                    ss1 << out_prefix << name.str();
//                    fout.open(ss1.str().c_str());
//                    if (out_bed_file!=NULL) {
//                        fout_bed << chr << "\t" << begin << "\t" << end << "\n";
//                    }
//
//                    BOOST_FOREACH(TmpA *a, vs) {
//                        if (a->pos>=begin && a->pos<end) {
//                            fout << a->data << "\n";
//                        } else {
//                            break;
//                        }
//                    }
//                    while (!vs.empty() && vs.front()->pos < begin+step) {
//                        delete(vs.front());
//                        vs.pop_front();
//                    }
//                }
//                fout << line << "\n";
//                if ( pos >= begin+step ) {
//                    TmpA *a = new TmpA;
//                    a->chr  = chr;
//                    a->pos  = pos;
//                    a->data = line;
//                    vs.push_back(a);
//                }
//            }
//        }
//        if (fout.is_open()) {
//            BOOST_FOREACH(TmpA *a, vs) {
//                if (a->pos>=begin && a->pos<end) {
//                    fout << a->data << "\n";
//                } else {
//                    break;
//                }
//            }
//            while (!vs.empty()) {
//                delete(vs.front());
//            }
//            vs.clear();
//            fout.close();
//        }
//    }
//    if (out_bed_file!=NULL) {
//        fout_bed.close();
//    }
//
//    return 0;
//}
//// }}}

/// class Cufflink
//// {{{
//void Cufflink::get_all_one_sample_tracking(istream & fin, GTP *gtp, vector<Express*> * exps)
//{
//    string line;
//    /// genes.fpkm_tracking filet header
//    /// racking_id    class_code    nearest_ref_id    gene_id    gene_short_name    
//    /// tss_id    locus    length    coverage    FPKM    
//    /// FPKM_conf_lo    FPKM_conf_hi    FPKM_status
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            continue;
//        } else {
//            vector<string> tab;
//            boost::split(tab, line, boost::is_any_of("\t"));
//            boost::trim(tab[6]);
//            vector<string> tab1;
//            boost::split(tab1, tab[6], boost::is_any_of(":-"));
//            for (Int i=0; i<tab1.size(); i++) { boost::trim(tab1[i]); }
//            Chrom *chr = gtp->get_chrom_by_symbol(tab1[0]);
//            Int b = atol(tab1[1].c_str())-1;
//            Int e = atol(tab1[2].c_str());
//            Express *exp_pt = new Express(chr);
//            exp_pt->ins_range(b, e);
//
//            double fpkm = atof(tab[9].c_str());
//            double fpkm_low = atof(tab[10].c_str());
//            double fpkm_high = atof(tab[11].c_str());
//            exp_pt->set_fpkm(fpkm, fpkm_low, fpkm_high);
//            exps->push_back(exp_pt);
//        }
//    }
//}
//
//void Cufflink::get_all_one_sample_tracking(istream & fin, GTP *gtp, vector<Express*> & exps)
//{
//    get_all_one_sample_tracking(fin, gtp, &exps);
//}
//
//void Cufflink::get_all_one_sample_tracking(istream & fin, vector<FPKMInfo*> & fpkm_infos)
//{
//    string line;
//    /// genes.fpkm_tracking filet header
//    /// racking_id    class_code    nearest_ref_id    gene_id    gene_short_name    
//    /// tss_id    locus    length    coverage    FPKM    
//    /// FPKM_conf_lo    FPKM_conf_hi    FPKM_status
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            continue;
//        } else {
//            vector<string> tab;
//            boost::split(tab, line, boost::is_any_of("\t"));
//            boost::trim(tab[6]);
//            vector<string> tab1;
//            boost::split(tab1, tab[6], boost::is_any_of(":-"));
//            for (Int i=0; i<tab1.size(); i++) { boost::trim(tab1[i]); }
//       //   Chrom *chr = gtp->get_chrom_by_symbol(tab1[0]);
//            Int b = atol(tab1[1].c_str())-1;
//            Int e = atol(tab1[2].c_str());
//       //   Express *exp_pt = new Express(chr);
//       //   exp_pt->ins_range(b, e);
//            FPKMInfo *info = new FPKMInfo;
//            fpkm_infos.push_back(info);
//
//            info->chrom = tab1[0];
//            info->begin = b;
//            info->end   = e;
//            info->strand= 0;
//
//            boost::trim(tab[0]);
//            info->gene_id = tab[0];
//            for (Int i=9; i<tab.size(); ) {
//                info->fpkm = atof(tab[i++].c_str());
//                info->lo_fpkm = atof(tab[i++].c_str());
//                info->hi_fpkm = atof(tab[i++].c_str());
//                string status = tab[i++];
//                if (status == "OK") {
//                    info->status = 'Y';
//                } else {
//                    info->status = 'N';
//                }
//                break;
//            }
//        }
//    }
//}
//
//void Cufflink::get_all_tracking(istream & fin, GTP *gtp, vector< vector<Express*> >* exps)
//{
//    string line;
//    /// genes.fpkm_tracking filet header
//    /// racking_id    class_code    nearest_ref_id    gene_id    gene_short_name    
//    /// tss_id    locus    length    coverage    FPKM    
//    /// FPKM_conf_lo    FPKM_conf_hi    FPKM_status
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            continue;
//        } else {
//            vector<string> tab;
//            boost::split(tab, line, boost::is_any_of("\t"));
//            boost::trim(tab[6]);
//            vector<string> tab1;
//            boost::split(tab1, tab[6], boost::is_any_of(":-"));
//            for (Int i=0; i<tab1.size(); i++) { boost::trim(tab1[i]); }
//            Chrom *chr = gtp->get_chrom_by_symbol(tab1[0]);
//            Int b = atol(tab1[1].c_str())-1;
//            Int e = atol(tab1[2].c_str());
//            Express *exp_pt = new Express(chr);
//            exp_pt->ins_range(b, e);
//
//            for (Int i=9, j=0; i<tab1.size(); j++) {
//                if ( j >= exps->size() ) {
//                    vector<Express*> v;
//                    exps->push_back(v);
//                }
//                double fpkm = atof(tab[i++].c_str());
//                double fpkm_low = atof(tab[i++].c_str());
//                double fpkm_high = atof(tab[i++].c_str());
//                i++;
//                exp_pt->set_fpkm(fpkm, fpkm_low, fpkm_high);
//                (*exps)[j].push_back(exp_pt);
//            }
//        }
//    }
//}
//
//void Cufflink::get_all_tracking(istream & fin, GTP *gtp, vector< vector<Express*> >& exps)
//{
//    get_all_tracking(fin, gtp, &exps);
//}
//
//void Cufflink::get_all_fpkm_from_one_sample_gtf(istream & fin, vector< FPKMInfo*> & fpkm_infos)
//{
//    string line;
//    /// GTF::type is 'transcript'
//    /// GTF::info: transcript_id, FPKM, conf_lo, conf_hi
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            continue;
//        } else {
//            vector<string> tab;
//            boost::split(tab, line, boost::is_any_of("\t"));
//            if (tab[2] != "transcript") continue;
//
//            boost::trim(tab[8]);
//            vector<string> tab1;
//            boost::split(tab1, tab[8], boost::is_any_of(";"));
//            for (Int i=0; i<tab1.size(); i++) { boost::trim(tab1[i]); }
//
//            FPKMInfo * info = new FPKMInfo;
//            fpkm_infos.push_back(info);
//            boost::trim(tab[0]);
//            info->chrom = tab[0];
//            info->begin = atoi(tab[3].c_str())-1;
//            info->end   = atoi(tab[4].c_str());
//            if (tab[6]=="+") {
//                info->strand = 1;
//            } else if (tab[6]=="-") {
//                info->strand = -1;
//            } else {
//                info->strand = 0;
//            }
//
//            for (Int j=0; j<tab1.size(); j++) {
//                vector<string> tab2;
//                boost::split(tab2, tab1[j], boost::is_any_of(" ="));
//                for (Int i=0; i<tab2.size(); i++) { 
//                    boost::trim_if(tab2[i], boost::is_any_of(" \"")); 
//                }
//                if (tab2[0]=="transcript_id") {
//                    info->tscript_id = tab2[1];
//                } else if (tab2[0]=="gene_id") {
//                    info->gene_id = tab2[1];
//                } else if (tab2[0]=="FPKM") {
//                    info->fpkm = atof(tab2[1].c_str());
//                } else if (tab2[0]=="conf_lo") {
//                    info->lo_fpkm = atof(tab2[1].c_str());
//                } else if (tab2[0]=="conf_hi") {
//                    info->hi_fpkm = atof(tab2[1].c_str());
//                } else if (tab2[0]=="frac") {
//                    info->frac = atof(tab2[1].c_str());
//                } else if (tab2[0]=="cov") {
//                    info->cov = atof(tab2[1].c_str());
//                }
//            }
//        }
//    }
//}
//
//void Cufflink::split_tracking_by_chrom(string & tracking_file, string & out_dir, vector<string> & out_chroms)
//{
//    split_tracking_by_chrom(tracking_file, out_dir, &out_chroms);
//}
//
//void Cufflink::split_tracking_by_chrom(string & tracking_file, string & out_dir, vector<string> * out_chroms)
//{
//    ifstream fin(tracking_file.c_str());
//    if ( !fin.is_open() ) {
//        stringstream ss;
//        ss << "In Cufflink::split_tracking_by_chrom(" << endl;
//        ss << "  tracking_file=" << tracking_file << ", out_dir=" << out_dir << ", out_chrom=" << out_chroms << "): ";
//        ss << "Can't opne FILE " << tracking_file << endl;
//        Msg::error(ss);
//    }
//    string line, file;
//    file = out_dir + "/header";
//    ofstream fout(file.c_str());
//    short r = 0;
//    vector<string> tab;
//    string chrom;
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            if (r==0) {
//                fout << line << endl;
//            }
//        } else {
//            if (r==0) {
//                fout.close();
//                r=1;
//            }
//            boost::split(tab, line, boost::is_any_of("\t"));
//            boost::trim(tab[6]);
//            vector<string> tab1;
//            boost::split(tab1, tab[6], boost::is_any_of(":-"));
//            boost::trim(tab1[6]);
//            if (tab1[0]!=chrom) {
//                if (fout.is_open()) fout.close();
//                chrom = tab1[0];
//                file = out_dir +"/" + chrom;
//                fout.open(file.c_str());
//                if (out_chroms!=NULL) {
//                    out_chroms->push_back(chrom);
//                }
//            }
//            fout << line << endl;
//        }
//    }
//    if (fout.is_open()) fout.close();
//}
//
///*
//void Cufflink::split_tracking_by_pos(string & tracking_file, string & out_dir, Int len, vector<string> * out_poss)
//{
//    ifstream fin(tracking_file.c_str());
//    if ( !fin.is_open() ) {
//        stringstream ss;
//        ss << "In Cufflink::split_tracking_by_chrom(" << endl;
//        ss << "  tracking_file=" << tracking_file << ", out_dir=" << out_dir << ", out_chrom=" << out_chroms << "): ";
//        ss << "Can't opne FILE " << tracking_file << endl;
//        Msg::error(ss);
//    }
//    string line, file;
//    file = out_dir + "/header";
//    ofstream fout(file.c_str());
//    short r = 0;
//    vector<string> tab;
//    string chrom;
//    Int pos1 = len;
//    Int pos;
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            if (r==0) {
//                fout << line << endl;
//            }
//        } else {
//            if (r==0) {
//                fout.close();
//                r=1;
//            }
//            boost::split(tab, line, boost::is_any_of("\t"));
//            boost::trim(tab[6]);
//            vector<string> tab1;
//            boost::split(tab1, tab[6], boost::is_any_of(":-"));
//            boost::trim(tab1[6]);
//            if (tab1[0]!=chrom) {
//                if (fout.is_open()) fout.close();
//                chrom = tab1[0];
//                file = out_dir +"/" + chrom;
//                fout.open(file.c_str());
//                if (out_chroms!=NULL) {
//                    out_chroms->push_back(chrom);
//                }
//            }
//            fout << line << endl;
//        }
//    }
//    if (fout.is_open()) fout.close();
//}
//*/
//
//void Cufflink::split_gtf_by_chrom(string & gtf_file, string & out_dir, vector<string> & out_chroms)
//{
//    split_gtf_by_chrom(gtf_file, out_dir, &out_chroms);
//}
//
//void Cufflink::split_gtf_by_chrom(string & gtf_file, string & out_dir, vector<string> * out_chroms)
//// {{{
//{
//    ifstream fin(gtf_file.c_str());
//    if ( !fin.is_open() ) {
//        stringstream ss;
//        ss << "In Cufflink::split_gtf_by_chrom(" << endl;
//        ss << "  gtf_file=" << gtf_file << ", out_dir=" << out_dir << ", out_chrom=" << out_chroms << "): ";
//        ss << "Can't opne FILE " << gtf_file << endl;
//        Msg::error(ss);
//    }
//    string line, file;
//    file = out_dir + "/header";
//    ofstream fout(file.c_str());
//    short r = 0;
//    vector<string> tab;
//    string chrom;
//    while (!fin.eof()) {
//        getline(fin, line);
//        boost::trim(line);
//        if (line.empty()) {
//            continue;
//        } else if (line[0]=='#') {
//            if (r==0) {
//                fout << line << endl;
//            }
//        } else {
//            if (r==0) {
//                fout.close();
//                r=1;
//            }
//            boost::split(tab, line, boost::is_any_of("\t"));
//            boost::trim(tab[0]);
//            if (tab[0]!=chrom) {
//                if (fout.is_open()) fout.close();
//                chrom = tab[0];
//                file = out_dir +"/" + chrom;
//                fout.open(file.c_str());
//                if (out_chroms!=NULL) {
//                    out_chroms->push_back(chrom);
//                }
//            }
//            fout << line << endl;
//        }
//    }
//    if (fout.is_open()) fout.close();
//}
//// }}}
//
//bool less_fpkminfo_gene_id(Cufflink::FPKMInfo* a, Cufflink::FPKMInfo* b)
//{
//    if ( a->gene_id < b->gene_id ) return true;
//    else if ( a->gene_id > b->gene_id ) return false;
//    
//    return false;
//}
//
//bool less_fpkminfo_pos(Cufflink::FPKMInfo* a, Cufflink::FPKMInfo* b)
//{
//    if ( a->chrom < b->chrom ) return true;
//    else if ( a->chrom > b->chrom ) return false;
//
//    if ( a->begin < b->begin ) return true;
//    else if ( a->begin > b->begin ) return false;
//
//    if ( a->end < b->end ) return true;
//    else if ( a->end > b->end ) return false;
//
//    if ( a->strand < b->strand ) return true;
//    else if ( a->strand > b->strand ) return false;
//
//    return false;
//}
//// }}}

/// class CZLDB
/*
void CZLDB::split_by_pos(string &db_dir, string &out_dir, Int split_len, string * out_pos_file)
{
    string line;
    vector<string> sv;
    Int i;
    while (!fin.eof()) {
        getline(fin, line);
        boost::trim(line);
        if (line.size()==0) continue;
        if (line[0]=='#' && line[1]=='#') {
        }
    }
}
*/

bool is_gz(const char *fname)
{
    int i=0;
    ifstream fin(fname, ios::binary);
    uint8_t byte1=0, byte2=0;
    fin.read((char*)&byte1, sizeof(byte1));
    if (fin.fail()) return false; 
    fin.read((char*)&byte2, sizeof(byte2));
    if (fin.fail()) return false; 
    fin.close();
    return byte1==0x1f && byte2==0x8b;
}

bool is_gz(const string & fname)
{
    return is_gz(fname.c_str());
}

};
