/** 
 * @file     io.hh
 * @details  This file include all IO function and classes
 * @author   Zelin
 */
#ifndef CZL_IO_H
#define CZL_IO_H

#include "czl_common.h"
#include "czl_bio_base.h"
//#include "czl_bio_annot.h"
#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string_regex.hpp"
#include "boost/foreach.hpp"
#include "boost/filesystem.hpp"
#include "boost/lexical_cast.hpp"
#include <gzstream.h>

using namespace std;

namespace czl_bio {

//class BioSeq;
//class BioFileSeq;

class Fasta {
public:
//    static void get_all_seq(ifstream & fin, vector<BioSeq> & seqs);
//    static void get_all_seq(ifstream & fin, vector<BioFileSeq> & seqs, string & file);
//    static ifstream & get_a_seq(ifstream & fin, BioSeq & seqs);
//    static ifstream & get_a_seq(ifstream & fin, BioFileSeq & seqs, string & file);
    static istream & get_a_seq(istream & fin, string & id, string & seq);

//    static ofstream & put_all_seq(ofstream & fout, vector<BioSeq> & seqs);
    static ostream & put_a_seq(ostream & fout, const string & id, const string & seq);
//    static ofstream & put_a_seq(ofstream & fout, BioSeq & seq);

    static int split_by_id(string & file, string & out_dir, vector<string> & out_ids);
    static int split_by_id_pos(string & fasta_file, string & out_prefix, int length, int step, string * out_bed_file=NULL);
};

/*
class VCF {
public:
    static string get_version(ifstream & fin);
    static string get_version(const string & file);

    static void get_all_variant_one_sample(ifstream & fin, GTP *gtp, vector<Variant*> & variants);
    static void get_all_variant_one_sample(ifstream & fin, string & version, GTP *gtp, vector< Variant* > & variants);
    static void get_all_variant_only_seq_region(ifstream & fin, GTP *gtp, vector<Variant*> & variants);
    static void get_all_variant_only_seq_region(ifstream & fin, string & version, GTP *gtp, vector<Variant*> & variants);
    static void get_all_variant(istream & fin, string & version, GTP *gtp, vector<Variant*> & variants);
//    static void get_a_variant(istream & fin, Variant & v);

    static int split_by_chr(const string & file, const string & out_dir, vector<string> & out_chroms);
    static int split_by_pos(const string & vcf_file, const string & out_dir, Int len, string * out_pos_file=NULL);
    static int split_by_pos(const string & vcf_file, const string & out_dir, Int len, string & out_pos_file);
    static int split_by_chr_pos(const string & vcf_file, int is_sorted, const string & out_prefix, int64_t length, int64_t step, string * out_bed_file=NULL);
};

class Cufflink {
public:

    class FPKMInfo {
    public:
        string tscript_id;
        string gene_id;
        float fpkm;
        float lo_fpkm;
        float hi_fpkm;
        float frac;
        float cov;
        char status;
        string chrom;
        Int begin;
        Int end;
        short strand;
    };


    static void get_all_one_sample_tracking(istream & fin, GTP *gtp, vector<Express*> * exps);
    static void get_all_one_sample_tracking(istream & fin, GTP *gtp, vector<Express*> & exps);
    static void get_all_one_sample_tracking(istream & fin, vector<FPKMInfo*> & fpkm_infos);

    static void get_all_tracking(istream & fin, GTP *gtp, vector< vector<Express*> >* exps);
    static void get_all_tracking(istream & fin, GTP *gtp, vector< vector<Express*> >& exps);

    static void get_all_fpkm_from_one_sample_gtf(istream & fin, vector< FPKMInfo* > & fpkm_infos);

    static void split_tracking_by_chrom(string & tracking_file, string & out_dir, vector<string> * out_chroms=NULL);
    static void split_tracking_by_chrom(string & tracking_file, string & out_dir, vector<string> & out_chroms);

    static void split_gtf_by_chrom(string & gtf_file, string & out_dir, vector<string> * out_chroms=NULL);
    static void split_gtf_by_chrom(string & gtf_file, string & out_dir, vector<string> & out_chroms);
};
bool less_fpkminfo_gene_id(Cufflink::FPKMInfo* a, Cufflink::FPKMInfo* b);
bool less_fpkminfo_pos(Cufflink::FPKMInfo* a, Cufflink::FPKMInfo* b);
*/

class CZLDB {
//  static void split_by_pos(string &db_dir, string &out_dir, Int split_len=0, string * out_pos_file=NULL);
};

/// class IOMergeSort
/// @brief template class for IO merge sort
/// Type T must be have a '=';
template<typename T>
class IOMergeSort {
private:
    int (*read_fun)(istream &, T&);
    int (*write_fun)(ostream &, const T&);
    bool (*less_fun)(const T&, const T&);

public:
    class A {
    public:
        T data;
        int id;
        A()
        {
            id = -1;
        }
        A(const A & a)
        {
            data = a.data;
            id   = a.id;
        }
        bool operator<(const A & a)
        {
            return less_fun(this->data, a.data);
        }
        bool operator>(const A & a)
        {
            return less_fun(a.data, this->data);
        }
        istream & operator >> (istream & s)
        {
            read_fun(s, data);
            return s;
        }
        ostream & operator << (ostream & s)
        {
            write_fun(s, data);
            return s;
        }

        A & operator = (const A & a)
        {
            data = a.data;
            id   = a.id;
            return (*this);
        }
    };

    IOMergeSort()
    {
        this->read_fun  = NULL;
        this->write_fun = NULL;
        this->less_fun  = NULL;
    }

    IOMergeSort(int (*read_fun)(istream &, T&), int (*write_fun)(ostream &, const T&), bool (*less_fun)(const T&, const T&))
    {
        this->read_fun  = read_fun;
        this->write_fun = write_fun;
        this->less_fun  = less_fun;
    }

    class Less {
    public:
        Less(bool (*less)(const T&, const T&)) {
            less_fun = less;
        }
        bool operator() (const A & a, const A & b)
        {
            return less_fun(a.data, b.data);
        }

        bool (*less_fun)(const T&, const T&);
    };

    class Great {
    public:
        Great(bool (*less)(const T&, const T&)) {
            less_fun = less;
        }
        bool operator() (const A & a, const A & b)
        {
            return less_fun(b.data, a.data);
        }

        bool (*less_fun)(const T&, const T&);
    };

    int run(const string & in, const string & out, const string & tmp_dir, size_t slice_size, size_t merge_size, short is_io_binary=0)
    {
        Less less(less_fun);
        Great great(less_fun);
        std::ios_base::openmode mode;
        if (is_io_binary) {
            mode |= std::ios::binary;
        }
        ifstream fin(in.c_str(), mode);
        if ( fin.fail() ) {
            cerr << "Can't open FILE " << in << endl;
            return 1;
        }
        int m=0;
        size_t n=0;
        vector<A> ts;
        while (!fin.eof()) {
            if (n>=slice_size) {
                sort(ts.begin(), ts.end(), less);
                stringstream ss;
                ss << tmp_dir << "/0_" << m;
                string file = ss.str();
                ofstream fout( file.c_str(), mode );
                for (size_t i=0; i<n; i++) write_fun(fout, ts[i].data);
                fout.close();
                ts.clear();
                n=0;
                m++;
            }
            A t;
            read_fun(fin, t.data);
            if (fin.fail()) break;
            ts.push_back(t);
            n++;
        }
        if (ts.size()>0) {
            sort(ts.begin(), ts.end(), less);
            stringstream ss;
            ss << tmp_dir << "/0_" << m;
            ofstream fout( ss.str().c_str(), mode);
            for (size_t i=0; i<n; i++) write_fun(fout, ts[i].data);
            fout.close();
            ts.clear();
            n=0;
            m++;
        }
        fin.close();

        /// merge
        // {{{
        int run=0;
        int buf_size = slice_size/(merge_size+1)+1;
    //  char *in_buf[merge_size];
        while (m>1) {
            int m1 = 0;
            for (int i1=0; i1<m; i1+=merge_size) {
                stringstream ss;
                ss << tmp_dir << "/" << run+1 << "_" << m1 ;
                string out_name1 = ss.str();

                int i2 = i1+merge_size<m ? i1+merge_size : m;
                if (i2==i1+1) {
                    stringstream ss1;
                    ss1 << tmp_dir << "/" << run << "_" << i1;
                    string in_name1 = ss1.str();
                    boost::filesystem::rename(in_name1, out_name1);
                    m1++;
                    break;
                }
                int mm = i2-i1;
                int write_n=0;
                ifstream fin1[mm];
                for (int i=i1; i<i2; i++) {
                    stringstream ss1;
                    ss1 << tmp_dir << "/" << run << "_" << i;
                    string in_name1 = ss1.str();
                    int j=i-i1;
                    fin1[j].open(in_name1.c_str(), mode);
                //  fin1[j].rdbuf()->pubsetbuf(in_buf[j], N);
                }
                ofstream fout1(out_name1.c_str(), ios::binary);
            //  fout1.rdbuf()->pubsetbuf(out_buf, N);

                vector<A> heap;
                A t;
                for (int i=0; i<mm; i++) {
                    if ( !fin1[i].eof() ) {
                        A t;
                        read_fun(fin1[i], t.data);
                        if ( !fin1[i].fail() ) {
                            t.id = i;
                            heap.push_back(t);
                        }
                    }
                }
                make_heap(heap.begin(), heap.end(), great);
                while (!heap.empty()) {
                    t = heap[0];
                    int i = t.id;
                    pop_heap(heap.begin(), heap.end(), great);
                    heap.pop_back();
                    write_fun(fout1, t.data);
                    write_n++;
                    if ( !fin1[i].eof() ) {
                        A t;
                        read_fun(fin1[i], t.data);
                        if ( !fin1[i].fail() ) {
                            t.id = i;
                            heap.push_back(t);
                            push_heap(heap.begin(), heap.end(), great);
                        }
                    }
                }
                fout1.close();
                for (int i=0; i<mm; i++) {
                    fin1[i].close();
                    stringstream ss1;
                    ss1 << tmp_dir << "/" << run << "_" << i+i1;
                    string in_name1 = ss1.str();
                    boost::filesystem::remove(in_name1);
                }
                m1++;
            }
            m = m1;
            run++;
        }
        {
            stringstream ss;
            ss << tmp_dir << "/" << run << "_0";
            string tmp=ss.str();
            if (boost::filesystem::exists(out)) {
                boost::filesystem::remove(out);
            }
            if (boost::filesystem::exists(tmp)) {
                boost::filesystem::rename(tmp, out);
            }
        }
        // }}}
        return 0;
    }

    ~IOMergeSort()
    {
    }

};

/* 
 * check if a fname is .gz
 */
bool is_gz(const char *fname);
bool is_gz(const string & fname);

/*
 * check if big-endian
 */
bool is_be()
{
	uint8_t a[4] = {0x01, 0x0, 0x0, 0x0};
	uint32_t b;
	memcpy(&b, a, sizeof(uint32_t));
	return b!=1;
}

/*
 * return the number of address bits on the current machine
 */
int machine_bit()
{
	return sizeof(void*)*8;
}

};

#endif
