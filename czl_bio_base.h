/* czl_bio_base.h
 * @brief basic class, basic definition, tables
 */
#ifndef CZL_BIO_BASE_H
#define CZL_BIO_BASE_H

#include "czl_common.h"

#ifdef PTHREAD
#include <pthread.h>
#endif

using namespace std;

namespace czl_bio {

    static const uint16_t N_NT0 = 4;
    static const uint16_t N_NT = 16;
//  static const char NT1[N_NT+1] = "TCAGYRMKSWHBVDN?-";
//  static const string NT4[N_NT] = {"T","C","A","G", "TC","AG","CA","TG","CG","TA", "TCA","TCG","CAG","TAG", "TCAG", "?", "-"};
    /// 8-bit id: bit 0-3: A,C,G,T; 4: 0 for capital, 1 for lower case (mask), 5-7: not use
    static const char BYTE_TO_NT1[N_NT*2+1] = "-ACMGRSVTWYHKDBN-acmgrsvtwyhkdbn";
    static const string BYTE_TO_NT4[N_NT*2] = {"-","A","C","AC","G", "AG","CG","ACG","T","AT","CT", "ACT", "GT","AGT","CGT", "ACGT", "-","a","c","ag","g", "ag","cg","acg","t","at","ct", "act", "gt","agt","cgt", "acgt"};
    static map<char,uint8_t> NT1_TO_BYTE = CreateMap<char,uint8_t>('-', 0x00)('A',0x01)('C',0x02)('M',0x03)('G',0x04)('R',0x05)('S',0x06)('V',0x07)('T',0x08)('W',0x09)('Y',0x0A)('H',0x0B)('K',0x0C)('D',0x0D)('B',0x0E)('N',0x0F)('?',0x0F)('-', 0x10)('a',0x11)('c',0x12)('m',0x13)('g',0x14)('r',0x15)('s',0x16)('v',0x17)('t',0x18)('w',0x19)('y',0x1a)('h',0x1b)('k',0x1c)('d',0x1d)('b',0x1e)('n',0x1f);
    static map<char,string> NT1_TO_NT4 = CreateMap<char,string>('A',"A")('C',"C")('G',"G")('T',"T")('R',"AG")('Y',"CT")('M',"AC")('K',"GT")('W',"AT")('S',"CG")('V',"ACG")('H',"ACT")('D',"AGT")('B',"CGT")('N',"ACGT")('?',"ACGT")('-',"-");
    static map<string,char> NT4_TO_NT1 = CreateMap<string,char>("A",'A')("C",'C')("G",'G')("T",'T')("AG",'R')("CT",'Y')("AC",'M')("GT",'K')("AT",'W')("CG",'S')("ACG",'V')("ACT",'H')("AGT",'D')("CGT",'B')("ACGT",'N')("ACGT",'?')("-",'-');
    static map<char,char> NT1_TO_COMPL= CreateMap<char,char>('A','T')('T','A')('C','G')('G','C')('R','Y')('Y','R')('K','M')('M','K')('W','W')('S','S')('H','D')('D','H')('V','B')('B','T')('a','t')('t','a')('c','g')('g','c')('r','y')('y','r')('k','m')('m','k')('w','w')('s','s')('h','d')('d','h')('v','b')('b','t');
	static const uint8_t BYTE_TO_COMPL[256]={0x00, 0x08, 0x04, 0x0C, 0x02, 0x0A, 0x06, 0x0E, 0x01, 0x09, 0x05, 0x0D, 0x03, 0x0B, 0x07, 0x0F,
		0x10, 0x18, 0x14, 0x1C, 0x12, 0x1A, 0x16, 0x1E, 0x11, 0x19, 0x15, 0x1D, 0x13, 0x1B, 0x17, 0x1F};

	/// amino acid
	/// Positive: R H K
	/// Negative: D E
	/// Polar unchared: S T N Q
	/// special: C U(selenocysteine) G P
	/// Hydrophobic: A B I L M F Y W
	/// 
	/// uncharged or chareged
	/// positive or negative
	/// aliphatic or aromatic
	/// cyclic or not
	/// sulphur or not
	/// OH or not
	/// small or large
    static const int N_AA0 = 21;
    static const int N_AA = 28;
	static map<char,string> AA1_TO_AA3=CreateMap<char,string>('A',"ALA")('B',"ASX")('C',"CYS")('D',"ASP")('E',"GLU")('F',"PHE")('G',"GLY")('H',"HIS")('I',"ILE")('J',"XLE")('K',"LYS")('L',"LEU")('M',"MET")('N',"ASN")('P',"PRO")('Q',"GLN")('R',"ARG")('S',"SER")('T',"THR")('U',"SEC")('V',"VAL")('W',"TRP")('X',"XXX")('Y',"TYR")('Z',"GLX")('*',"***")('-',"---")('?',"???");
	static map<string,char> AA3_TO_AA1=CreateMap<string,char>("ALA",'A')("ASX",'B')("CYS",'C')("ASP",'D')("GLU",'E')("PHE",'F')("GLY",'G')("HIS",'H')("ILE",'I')("XLE",'J')("LYS",'K')("LEU",'L')("MET",'M')("ASN",'N')("PRO",'P')("GLN",'Q')("ARG",'R')("SER",'S')("THR",'T')("SEC",'U')("VAL",'V')("TRP",'W')("XXX",'X')("TYR",'Y')("GLX",'Z')("***",'*')("---",'-')("???",'?');
	static map<char,string> AA1_TO_NAME=CreateMap<char,string>('A',"Alanine")('B',"AsparticAcidOrAsparagine")('C',"Cysteine")('D',"AsparticAcid")('E',"GlutamicAcid")('F',"Phenylalanine")('G',"Glycine")('H',"Histidine")('I',"Isoleucine")('J',"LeucineOrIsoleucine")('K',"Lysine")('L',"Leucine")('M',"Methionine")('N',"Asparagine")('P',"Proline")('Q',"Glutamine")('R',"ARG")('S',"Serine")('T',"Threonine")('U',"Selenocysteine")('V',"Valine")('W',"Tryptophan")('X',"Unknown")('Y',"Tyrosine")('Z',"GlutaminAcidOrGlutamine")('*',"Stop")('-',"Gap")('?',"Unknown");
    static map<string, char> GENERAL_CODON_TABLE = CreateMap<string,char>
        ("TTT",'F') ("TTC",'F') ("TTA",'L') ("TTG",'L')
        ("TCT",'S') ("TCC",'S') ("TCA",'S') ("TCG",'S')
        ("TAT",'Y') ("TAC",'Y') ("TAA",'*') ("TAG",'*')
        ("TGT",'C') ("TGC",'C') ("TGA",'*') ("TGG",'W')
        ("CTT",'L') ("CTC",'L') ("CTA",'L') ("CTG",'L')
        ("CCT",'P') ("CCC",'P') ("CCA",'P') ("CCG",'P')
        ("CAT",'H') ("CAC",'H') ("CAA",'Q') ("CAG",'Q')
        ("CGT",'R') ("CGC",'R') ("CGA",'R') ("CGG",'R')
        ("ATT",'I') ("ATC",'I') ("ATA",'I') ("ATG",'M')
        ("ACT",'T') ("ACC",'T') ("ACA",'T') ("ACG",'T')
        ("AAT",'N') ("AAC",'N') ("AAA",'K') ("AAG",'K')
        ("AGT",'S') ("AGC",'S') ("AGA",'R') ("AGG",'R')
        ("GTT",'V') ("GTC",'V') ("GTA",'V') ("GTG",'V')
        ("GCT",'A') ("GCC",'A') ("GCA",'A') ("GCG",'A')
        ("GAT",'D') ("GAC",'D') ("GAA",'E') ("GAG",'E')
        ("GGT",'G') ("GGC",'G') ("GGA",'G') ("GGG",'G')
		("GCN",'A') ("CGN",'R') ("MGR",'R') ("AAY",'N')
		("GAY",'D') ("TGY",'C') ("CAR",'Q') ("GAR",'E')
		("GGN",'G') ("CAY",'H') ("ATH",'I') ("YTR",'L')
		("CTN",'L') ("AAR",'K') ("TTY",'F') ("CCN",'P')
		("TCN",'S') ("AGY",'S') ("ACN",'T') ("TAY",'Y')
		("GTN",'V') ("TAR",'*') ("TRA",'*');
    static map<string, char> MT_CODON_TABLE = CreateMap<string,char>
        ("TTT",'F') ("TTC",'F') ("TTA",'L') ("TTG",'L') ("TTY",'F')("YTR",'L')
        ("TCT",'S') ("TCC",'S') ("TCA",'S') ("TCG",'S') ("TCN",'S')
        ("TAT",'Y') ("TAC",'Y') ("TAA",'*') ("TAG",'*') ("TAY",'Y')("TAR",'*')
        ("TGT",'C') ("TGC",'C') ("TGA",'W') ("TGG",'W') ("TGY",'C')("TGR",'W')
        ("CTT",'L') ("CTC",'L') ("CTA",'L') ("CTG",'L') ("CTN",'L')
        ("CCT",'P') ("CCC",'P') ("CCA",'P') ("CCG",'P') ("CCN",'P')
        ("CAT",'H') ("CAC",'H') ("CAA",'Q') ("CAG",'Q') ("CAY",'H')("CAR",'Q')
        ("CGT",'R') ("CGC",'R') ("CGA",'R') ("CGG",'R') ("CGN",'R')
        ("ATT",'I') ("ATC",'I') ("ATA",'M') ("ATG",'M') ("ATY",'I')("ATR",'M')
        ("ACT",'T') ("ACC",'T') ("ACA",'T') ("ACG",'T') ("ACN",'T')
        ("AAT",'N') ("AAC",'N') ("AAA",'K') ("AAG",'K') ("AAY",'N')("AAR",'K')
        ("AGT",'S') ("AGC",'S') ("AGA",'*') ("AGG",'*') ("AGY",'S')("AGR",'*')
        ("GTT",'V') ("GTC",'V') ("GTA",'V') ("GTG",'V') ("GTN",'V')
        ("GCT",'A') ("GCC",'A') ("GCA",'A') ("GCG",'A') ("GCN",'A')
        ("GAT",'D') ("GAC",'D') ("GAA",'E') ("GAG",'E') ("GAY",'D')("GAR",'E')
        ("GGT",'G') ("GGC",'G') ("GGA",'G') ("GGG",'G') ("GGN",'G') ;
    static const int N_CN0 = 64;
//  static const string AA3[N_AA] = {"ALA","ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASX", "XLE", "GLX", "XXX", "***", "---", "???", "UUU"};
//  static const char AA1[N_AA+1] = "ARNDCQEGHILKMFPSTWYVBJZX*-?U";
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

	/// byte to other
    char byte_to_nt1(uint8_t byte);
	string byte_to_nt4(uint8_t byte);
	/// bit4 (byte's lower 4 bits) to other
	char bit4_to_nt1(uint8_t bit4);
	string bit4_to_nt4(uint8_t bit4);
	/// bit2
	char bit2_to_nt1(uint8_t bit2);
	///
	string nt1_to_nt4(char nt1);
	uint8_t nt1_to_byte(char nt1);
	uint8_t nt1_to_bit4(char nt1);
	char nt4_to_nt1(const string & nt4);
	uint8_t nt4_to_byte(const string & nt4);
	//
	template<typename Iter>
	int bit4_to_nt1(const Iter & begin, const int num, string & out)
	{
		if (num<=0) return 0;
		Iter it = begin;
		int nb=sizeof(*begin)*8;
		int nb1=0;
		int k=0;
		int l=4;
		typename iterator_traits<Iter>::value_type a;
		uint8_t b4 = 0;
		a = (*it);
		for (int i=0; i<num; i++) {
			b4 = (a&0xf);
			a >>= l;
			k+=4;
			if (k==nb) {
				it++;
				a=(*it);
				k=0;
			} if (k>nb) {
				it++;
				a=(*it);
				k-=nb;
				b4 |= ((a<<(l-k))&0xf);
				a >>= k;
			}
			out[i] = bit4_to_nt1(b4);
		}

		return 0;
	}

	template<typename T>
	int nt1_to_bit4(const string::iterator & begin, const string::iterator & end, T *out_begin)
	{
		if (begin==end) return 0;
		T a=0;
		int nb=sizeof(a)*8;
		int l=4;
		int k=0;
		for (string::iterator it=begin; it!=end; it++) {
			uint8_t b = nt1_to_bit4((*it));
			a |= (b<<k);
			if (k+l==nb) {
				(*out_begin++) = a;
				a = 0;
				k = 0;
			} else if (k+l>nb) {
				(*out_begin++) = a;
				a = 0;
				a |= b>>(nb-k);
				k = k+l-nb;
			} else {
				k += l;
			}
		}
		if (k>0) (*out_begin++) = a;
		return 0;
	}

	template<typename Iter>
	int nt1_to_bit4(const string::iterator & begin, const string::iterator & end, Iter & out_begin)
	{
		if (begin==end) return 0;
		typename Iter::value_type a = 0;
		int nb=sizeof(a)*8;
		int l=4;
		int k=0;
		for (string::iterator it=begin; it!=end; it++) {
			uint8_t b = nt1_to_bit4((*it));
			a |= (b<<k);
			if (k+l==nb) {
				(*out_begin++) = a;
				a = 0;
				k = 0;
			} else if (k+l>nb) {
				(*out_begin++) = a;
				a = 0;
				a |= b>>(nb-k);
				k = k+l-nb;
			} else {
				k += l;
			}
		}
		if (k>0) (*out_begin++) = a;
		return 0;
	}

	template<typename Iter>
	void print_bit4_as_nt1(const Iter & begin, const int num)
	{
		if (num<=0) return ;
		Iter it = begin;
		int nb=sizeof(*begin)*8;
		int nb1=0;
		int k=0;
		int l=4;
		typename iterator_traits<Iter>::value_type a;
		uint8_t b4 = 0;
		a = (*it);
		for (int i=0; i<num; i++) {
			b4 |= (a&0xf);
			a >>= l;
			k+=4;
			if (k>nb) {
				it++;
				a=(*it);
				k-=nb;
				b4 |= (a<<(l-k));
				a >>= k;
			}
			cout << bit4_to_nt1(b4);
		}
	}

	///

	char complement_nt1(char nt1);
	uint8_t complement_byte(uint8_t byte);
	int complement_copy(const string & seq, string & out_seq);
	int complement(string & seq);
	int rev_complement_copy(const string & seq, string & out_seq);
	int rev_complement(string & seq);

	string aa1_to_aa3(char aa1);
	string aa1_to_name(char aa1);
	char aa3_to_aa1(const string & aa3);

	char codon_to_aa1(char const codon[3], const map<string,char> & codon_table=GENERAL_CODON_TABLE);
	char codon_to_aa1(const string & codon, const map<string,char> & codon_table=GENERAL_CODON_TABLE);
	int translate(const string & codon, string & prot, const map<string,char> & codon_table=GENERAL_CODON_TABLE);

	void load_aa_score_mat(string & file, vector< vector<float> >& out_mat);
	void load_aa_score_mat(string & file, int* out_n, float ***out_mat);
    /*
    void init_codon_table();
    char complement1(char nt1);
    Int complement(string & seq);
    Int complement_copy(string & seq, string & tag);
    char seq_codon_to_aa1(string & seq);
	char seq_codon_to_id(string & codon);
    Int seq_translate(string & seq, string & tag);
    Int seq_to_rev_complement(string & seq);
    Int seq_to_rev_complement_copy(string & seq, string & tag);
	Int seq_aa1_to_id(string & aa1);
	Int seq_aa1_to_id(char aa1);
	Int seq_nt1_to_id(string & nt1);
	Int seq_nt1_to_id(char nt1);
    */


/** class BioObj
 * @brief    Basic Object for all biological element
 * @details  This class is used for the basic classes of all other element.
 *           It includes an ID, symbol and a bitset flag, 
 *           you can set the flag by get_flag() first, 
 *           which will return the REF of the bitset flag
 *
 * @author   Zelin
 */
class BioObj {
protected:
	const static size_t flag_size_m=32;
	void *data_m;
	uint64_t flag_m;
public:
	BioObj();
	BioObj(const BioObj & obj);
	virtual ~BioObj();

    virtual bool is_bit_set(size_t bit);
	virtual bool is_flag_set(uint64_t flag);
	virtual bool is_any_flag_set(uint64_t flag);
	virtual bool is_all_flag_set(uint64_t flag);
	virtual void * get_data();

	virtual void set(const BioObj & obj);
	virtual void set(BioObj & obj);
	virtual void set(BioObj * obj);
	virtual int set_bit(size_t bit);
	virtual void set_flag(uint64_t flag);
	virtual void set_flag();
	virtual int unset_bit(size_t bit);
	virtual void unset_flag(size_t bit);
	virtual void unset_flag();
	virtual void set_data(void * data);

	BioObj & copy(const BioObj & obj);

	BioObj & operator =(BioObj & obj);
};


/*
 * class BioRegion
 * @brief base class of any biology regions
 * @comment make sure sub-region is sorted before calling get_begin() 
 *          and get_end()
 */
/*
template<typename T>
class BioRegion : public BioObj {
protected:
	BioObj *obj_m;
	vector< pair<T,T> > reg_m; // regions
public:
	BioRegion(BioObj *obj=NULL): BioObj(), obj_m(obj) {}

	BioRegion(BioRegion& objreg): BioObj(objreg), obj_m(objreg.obj_m), reg_m(objreg.reg_m) {};
	BioRegion(const BioRegion& objreg): BioObj(objreg), obj_m(objreg.obj_m), reg_m(objreg.reg_m) {};
	virtual ~BioRegion() {}

	virtual BioObj* get_obj() { return obj_m; }
    vector< pair<T,T> > & get_reg_ref() { return reg_m; }
	T get_begin(T i) { return reg_m[i].first; }
	T get_begin() { return reg_m[0].first; }
	T get_end(size_t i) { return reg_m[i].second; }
	T get_end() { return reg_m[reg_m.size()-1].second; }
	T get_len(size_t i) { return reg_m[i].second-reg_m[i].first; }
	T get_len()
	{
		T len();
		for (size_t i=0; i<reg_m.size(); i++) {
			len+=reg_m[i].second-reg_m[i].first;
		}
		return len;
	}

	T get_span_len()
	{
		return get_end()-get_begin();
	}

	BioRegion& set(BioRegion& obj_region)
    {
        obj_m = obj_region.obj_m;
        reg_m = obj_region.reg_m;
        return (*this);
    }

	BioRegion& set(const BioRegion& obj_region)
    {
        obj_m = obj_region.obj_m;
        reg_m = obj_region.reg_m;
        return (*this);
    }

	virtual void set_obj(BioObj *obj) { obj_m = obj; }

	void push_back(T begin, T end)
    {
        reg_m.push_back(pair<T,T>(begin, end));
    }

    void sort(bool (*less)(const T&, const T&))
    {
        sort(reg_m.begin(), reg_m.end(), less);
    }
};
*/

/*
class BioSeq: virtual public BioObj {
protected:
    string seq_m;

public:
    BioSeq() {};
    BioSeq(string & seq): BioObj(id, symbol), seq_m(seq) {}
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
*/

};

#endif
