#ifndef CZL_BIO_SEQ_H
#define CZL_BIO_SEQ_H

#include "czl_common.h"

namespace czl_bio {

    class BioSeq;

    static map<string, char> MAP_CODON_TO_AA1_GENERAL;
    static const int N_AA0 = 20;
    static const int N_AA = 28;
    static const int N_NT0 = 4;
    static const int N_NT = 18;
    static const int N_CN0 = 64;
    static const string AA3[N_AA] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASX", "XLE", "GLX", "XXX", "***", "---", "???", "UUU"};
    static const char AA1[N_AA+1] = "ARNDCQEGHILKMFPSTWYVBJZX*-?U";
    static const char NT1[N_NT+1] = "TCAGYRMKSWHBVDN?-";
    static const string NT4[N_NT] = {"T","C","A","G", "TC","AG","CA","TG","CG","TA", "TCA","TCG","CAG","TAG", "TCAG", "?", "-"};
    static const string CODON[N_CN0] = {
        "TTT", "TTC", "TTA", "TTG",
        "TCT", "TCC", "TCA", "TCG",
        "TAT", "TAC", "TAA", "TAG",
        "TGT", "TGC", "TGA", "TGG",
        "CTT", "CTC", "CTA", "CTG",
        "CCT", "CCC", "CCA", "CCG",
        "CAT", "CAC", "CAA", "CAG",
        "CGT", "CGC", "CGA", "CGG",
        "ATT", "ATC", "ATA", "ATG",
        "ACT", "ACC", "ACA", "ACG",
        "AAT", "AAC", "AAA", "AAG",
        "AGT", "AGC", "AGA", "AGG",
        "GTT", "GTC", "GTA", "GTG",
        "GCT", "GCC", "GCA", "GCG",
        "GAT", "GAC", "GAA", "GAG",
        "GGT", "GGC", "GGA", "GGG"
    };
    static const char CODON_AA1_GENERAL[N_CN0+1] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    static map<char, int> MAP_AA1_TO_ID;
    static map<char, int> MAP_NT1_TO_ID;
    static map<string, int> MAP_CODON_TO_ID;

    void seq_init_codon_table();
    char seq_to_complement1(char nt1);
    Int seq_to_complement(string & seq);
    Int seq_to_complement_copy(string & seq, string & tag);
    char seq_codon_to_aa1(string & seq);
	char seq_codon_to_id(string & codon);
    Int seq_translate(string & seq, string & tag);
    Int seq_to_rev_complement(string & seq);
    Int seq_to_rev_complement_copy(string & seq, string & tag);
	Int seq_aa1_to_id(string & aa1);
	Int seq_aa1_to_id(char aa1);
	Int seq_nt1_to_id(string & nt1);
	Int seq_nt1_to_id(char nt1);
	void seq_load_aa_score_mat(string & file, vector< vector<float> >& out_mat);
	void seq_load_aa_score_mat(string & file, Int* out_n, float ***out_mat);

    class BioSeq: virtual public BioObj {
    protected:
        string seq_m;

    public:
        BioSeq() {};
        BioSeq(string & id, string & symbol, string & seq): BioObj(id, symbol), seq_m(seq) {}
        BioSeq(string & seq): seq_m(seq) {}
        virtual ~BioSeq() {};
        
        virtual Int get_seq_len() { return seq_m.size(); }
        virtual string get_seq(Int begin=0, Int len=-1);
        virtual Int get_seq(string & out, Int begin=0, Int len=-1);

        virtual void set_seq(const char *seq) { seq_m.assign(seq); }
        virtual void set_seq(const string & seq) { seq_m.assign(seq); }
		virtual void set(BioSeq & seq);
		virtual void set(const BioSeq & seq);

        virtual void app_seq(const char *seq) { seq_m += seq; }
        virtual void app_seq(const string & seq) { seq_m += seq; }

        BioSeq & operator = (BioSeq & seq);
        BioSeq & operator = (const BioSeq & seq);
    };

    namespace io {
        class Fasta;
    }

    class SeqRegion: virtual public ObjRegion {
    protected:
        short strand_m;

    public:
        SeqRegion(BioSeq *bio_seq=NULL, short strand=0): ObjRegion(bio_seq), strand_m(strand) {}
		SeqRegion( SeqRegion & seq_region );
		SeqRegion(const SeqRegion & seq_region);
		virtual ~SeqRegion() {}

        BioSeq* get_seq_pt() { return dynamic_cast<BioSeq*>(obj_pt_m); }
        short get_strand() { return strand_m; }
        string get_seq();
        Int get_seq(string & seq);
        Int get_seq_len(Int i);
        Int get_seq_len();

		SeqRegion & set(SeqRegion & seq_region);
		SeqRegion & set(const SeqRegion & seq_region);
        void set_strand(short strand) { strand_m = strand>0 ? 1:(strand<0? -1:0); }
        void set_seq_pt(BioSeq *seq_pt) { obj_pt_m = seq_pt; }

		SeqRegion & operator =(SeqRegion & seq_region);
		SeqRegion & operator =(const SeqRegion & seq_region);
    };
};

#endif
