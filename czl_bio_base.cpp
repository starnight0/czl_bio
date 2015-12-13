#include "czl_bio_base.hpp"

namespace czl_bio {

char byte_to_nt1(uint8_t byte)
{
    return BYTE_TO_NT1[byte];
}

string byte_to_nt4(uint8_t byte)
{
    return BYTE_TO_NT4[byte];
}

char bit4_to_nt1(uint8_t bit4)
{
    return BYTE_TO_NT1[bit4&0xf];
}
string bit4_to_nt4(uint8_t bit4)
{
    return BYTE_TO_NT4[bit4&0xf];
}

/// bit2: bit0: is pyrimidine; bit1: is strong (triple bond)
char bit2_to_nt1(uint8_t bit2)
{
	switch(bit2&0x3) {
	case 0x0: return 'A';
	case 0x1: return 'T';
	case 0x2: return 'G';
	case 0x3: return 'C';
	default: return 'N';
	}
	return 'N';
}

string nt1_to_nt4(char nt1)
{
    if (NT1_TO_NT4.find(nt1)!=NT1_TO_NT4.end()) return NT1_TO_NT4[nt1];
	else return "";
}

uint8_t nt1_to_byte(char nt1)
{
    if (NT1_TO_BYTE.find(nt1)!=NT1_TO_BYTE.end()) return NT1_TO_BYTE[nt1];
	else return 0xf;
}

uint8_t nt1_to_bit4(char nt1)
{
    if (NT1_TO_BYTE.find(nt1)!=NT1_TO_BYTE.end()) return NT1_TO_BYTE[nt1]&0xf;
	else return 0xf;
}

char nt4_to_nt1(const string & nt4)
{
    if (NT4_TO_NT1.find(nt4)!=NT4_TO_NT1.end()) return NT4_TO_NT1[nt4];
	else return 'N';
}

uint8_t nt4_to_byte(const string & nt4)
{
    if (NT4_TO_NT1.find(nt4)!=NT4_TO_NT1.end()) return NT1_TO_BYTE[NT4_TO_NT1[nt4]];
	else return 0xf;
}

char complement_nt1(char nt1)
{
    if (NT1_TO_COMPL.find(nt1)!=NT1_TO_COMPL.end()) return NT1_TO_COMPL[nt1];
	else return 'N';
}

uint8_t complement_byte(uint8_t byte)
{
    return BYTE_TO_COMPL[byte];
}

int complement_nt1_copy(const string & seq, string & out_seq)
{
	out_seq.resize(seq.size());
    for (int i=0; i<seq.size(); i++) {
        if (NT1_TO_COMPL.find(seq[i])==NT1_TO_COMPL.end()) {
            out_seq[i] = 'N';
        } else {
            out_seq[i] = NT1_TO_COMPL[seq[i]];
        }
    }
    return 0;
}

int complement_nt1(string & seq)
{
    for (int i=0; i<seq.size(); i++) {
        if (NT1_TO_COMPL.find(seq[i])==NT1_TO_COMPL.end()) {
            seq[i] = 'N';
        } else {
            seq[i] = NT1_TO_COMPL[seq[i]];
        }
    }
    return 0;
}

int rev_complement_nt1_copy(const string & seq, string & out_seq)
{
	out_seq.resize(seq.size());
    for (int i=0; i<seq.size(); i++) {
        if (NT1_TO_COMPL.find(seq[i])==NT1_TO_COMPL.end()) {
            out_seq[i] = 'N';
        } else {
            out_seq[i] = NT1_TO_COMPL[seq[i]];
        }
    }
    reverse(out_seq.begin(),out_seq.end());
    return 0;
}

int rev_complement_nt1(string & seq)
{
	size_t n = seq.size();
    for (int i=0; i<n/2; i++) {
		char c1, c2;
        if (NT1_TO_COMPL.find(seq[i])==NT1_TO_COMPL.end()) {
			c1 = 'N';
        } else {
			c1 = NT1_TO_COMPL[seq[i]];
        }
		int j = n-i-1;
        if (NT1_TO_COMPL.find(seq[j])==NT1_TO_COMPL.end()) {
			c2 = 'N';
        } else {
			c2 = NT1_TO_COMPL[seq[j]];
        }
		seq[i] = c2;
		seq[j] = c1;
    }
	if (n%2==1) {
		int i=n/2;
        if (NT1_TO_COMPL.find(seq[i])==NT1_TO_COMPL.end()) {
			seq[i] = 'N';
        } else {
			seq[i] = NT1_TO_COMPL[seq[i]];
        }
	}
    return 0;
}

string aa1_to_aa3(char aa1)
{
    if (AA1_TO_AA3.find(aa1)==AA1_TO_AA3.end()) return "";
    else return AA1_TO_AA3[aa1];
}

string aa1_to_name(char aa1)
{
    if (AA1_TO_NAME.find(aa1)==AA1_TO_NAME.end()) return "";
    else return AA1_TO_NAME[aa1];
}

char aa3_to_aa1(const string & aa3)
{
    if (AA3_TO_AA1.find(aa3)==AA3_TO_AA1.end()) return '\0';
    else return AA3_TO_AA1[aa3];
}

char codon_to_aa1(char const codon[3], map<string,char> & codon_table)
{
    uint8_t is_mask=0;
    string s(codon);
    for (int j=0; j<3; j++) {
        is_mask |= (s[j]&0x20);
        s[j]&=~0x20;
    }
    char a;
    if (codon_table.find(s)==codon_table.end()) {
        a = '\0';
    } else {
        a = (codon_table[s] | is_mask);
    }
    return a;
}

char codon_to_aa1(const string & codon, map<string,char> & codon_table)
{
    if (codon.size()<3) return '\0';
    uint8_t is_mask=0;
    string s = codon;
    for (int j=0; j<3; j++) {
        is_mask |= (s[j]&0x20);
        s[j]&=~0x20;
    }
    char a;
    if (codon_table.find(s)==codon_table.end()) {
        a = '\0';
    } else {
        a = (codon_table[s] | is_mask);
    }
    return a;
}

int translate(const string & codon, string & prot, map<string,char> & codon_table)
{
	int k=0;
	prot.resize(codon.size()/3);
    for (int i=0; i<codon.size()-2; i+=3, k++) {
        uint8_t is_mask=0;
        string s = codon.substr(i, 3);
        for (int j=0; j<3; j++) {
            is_mask |= (s[j]&0x20);
            s[j]&=~0x20;
        }
        if (codon_table.find(s)!=codon_table.end()) {
            char a = (codon_table[s] | is_mask);
            prot[k] = a;
        }
    }
    return 0;
}

/// read a protein alignment score matrix from file
void load_aa_score_mat(string & file, vector< vector<float> >& out_mat)
// {{{
{
    ifstream fin(file.c_str());
    vector<char> aa1s;
    int r=0;
    int n=0;
    while (!fin.eof()) {
        string line;
        getline(fin, line);
        if (line.empty()) continue;
        if (line[0]=='#') continue;

        vector<string> tab;
        if (r==0) {
            size_t i=0, j, k=0;
            while (1) {
                j = line.find("\t", i);
                if (j>i+1) {
                    char aa1 = line[i];
                    if (aa1 > n && AA1_TO_AA3.find(aa1)!=AA1_TO_AA3.end()) n = aa1;
                }
                if (j==string::npos) break;
                i=j+1;
            }
            r=1;
        } else {
            char aa1 = line[0];
            if (aa1 > n && AA1_TO_AA3.find(aa1)!=AA1_TO_AA3.end()) n = aa1;
        }
    }
    n++;
    out_mat.resize(n);
//  for (Int i=0; i<n; i++) out_mat[i].resize(n);
    fin.close();

    r=0;
    fin.open(file.c_str());
    while (!fin.eof()) {
        string line;
        getline(fin, line);
        if (line.empty()) continue;
        if (line[0]=='#') continue;

        vector<string> tab;
   //   boost::split(tab, line, boost::is_any_of("\t"));
        if (r==0) {
            size_t i=0, j, k=0;
            while (1) {
                j = line.find("\t", i);
                if (k>0) {
                    if (j==i+1) {
                        aa1s.push_back(-1);
                    } else if (j>i+1) {
                        char aa1 = line[i];
                        if (AA1_TO_AA3.find(aa1)!=AA1_TO_AA3.end()) {
                            aa1s.push_back(aa1);
                            out_mat[aa1].resize(n);
                        } else {
                            aa1s.push_back(-1);
                        }
                    }
                }
                if (j==string::npos) break;
                i=j+1;
            }
            r=1;
        } else {
        //  Int id1 = seq_aa1_to_id(tab[0]);
            size_t i=0, j, k=0;
            char aa1;
            while (1) {
                j = line.find("\t", i);
                if (k==0) { aa1 = line[i]; }
                else {
                    if (aa1s[k-1]>0) {
                        char aa2 = aa1s[k-1];
                        out_mat[aa1][aa2] = atof(line.substr(i,j-i).c_str());
                    }
                }
                if (j==string::npos) break;
                i=j+1;
                k++;
            //  Int id2 = aa1s[i-1];
            //  out_mat[id1][id2] = atof(tab[i].c_str());
            }
        }
    }
    fin.close();
}
// }}}

void load_aa_score_mat(string & file, int *out_n, float ***out_mat)
// {{{
{
    ifstream fin(file.c_str());
//  vector<Int> aa1s;
    vector<char> aa1s;
    int r=0;
    int n=0;
    while (!fin.eof()) {
        string line;
        getline(fin, line);
        if (line.empty()) continue;
        if (line[0]=='#') continue;

        vector<string> tab;
    //  boost::split(tab, line, boost::is_any_of("\t"));
        if (r==0) {
        //  for (size_t i=1; i<tab.size(); i++) {
        //      if (!tab[i].empty()) {
        //      //  Int id = seq_aa1_to_id(tab[i]);
        //      //  if (id > n) n = id;
        //          char aa = tab[i];
        //          if (aa > n) n = aa;
        //      }
        //  }
            size_t i=0, j, k=0;
            while (1) {
                j = line.find("\t", i);
                if (j>i+1) {
                    char aa1 = line[i];
                    if (aa1 > n && AA1_TO_AA3.find(aa1)!=AA1_TO_AA3.end()) n = aa1;
                }
                if (j==string::npos) break;
                i=j+1;
            }
            r=1;
        } else {
        //  Int id = seq_aa1_to_id(tab[0]);
        //  if (id > n) n = id;
            char aa1 = line[0];
            if (aa1 > n && AA1_TO_AA3.find(aa1)!=AA1_TO_AA3.end()) n = aa1;
        }
    }
    fin.close();

    n++;
    (*out_n) = n;
    (*out_mat) = new float*[n];

    r=0;
    fin.open(file.c_str());
    while (!fin.eof()) {
        string line;
        getline(fin, line);
        if (line.empty()) continue;
        if (line[0]=='#') continue;

        vector<string> tab;
   //   boost::split(tab, line, boost::is_any_of("\t"));
        if (r==0) {
        //  for (size_t i=1; i<tab.size(); i++) {
        //      if (tab[i].empty()) aa1s.push_back(-1);
        //  //  else aa1s.push_back(seq_aa1_to_id(tab[i]));
        //      else {
        //          aa1s.push_back(tab[i]);
        //          (*out_mat)[i] = new float[n];
        //      }
        //  }
            size_t i=0, j, k=0;
            while (1) {
                j = line.find("\t", i);
                if (k>0) {
                    if (j==i+1) {
                        aa1s.push_back(-1);
                    } else {
                        char aa1 = line[i];
                        if (AA1_TO_AA3.find(aa1)!=AA1_TO_AA3.end()) {
                            aa1s.push_back(aa1);
                            (*out_mat)[aa1] = new float[n];
                        } else {
                            aa1s.push_back(-1);
                        }
                    }
                }
                if (j==string::npos) break;
                i=j+1;
            }
            r=1;
        } else {
        //  Int id1 = seq_aa1_to_id(tab[0]);
        //  char aa1 = tab[0];
        //  for (size_t i=1; i<tab.size(); i++) {
        //      if (aa1s[i-1]<0) continue;
        //  //  Int id2 = aa1s[i-1];
        //  //  (*out_mat)[id1][id2] = atof(tab[i].c_str());
        //      aa2 = aa1s[i-1];
        //      (*out_mat)[aa1][aa2] = atof(tab[i].c_str());
        //  }
            size_t i=0, j, k=0;
            char aa1;
            while (1) {
                j = line.find("\t", i);
                if (k==0) {
                    aa1 = line[i];
                } else {
                    if (aa1s[k-1]>0) {
                        char aa2 = aa1s[k-1];
                        (*out_mat)[aa1][aa2] = atof(line.substr(i,j-i).c_str());
                    }
                }
                if (j==string::npos) break;
                i=j+1;
                k++;
            }
        }
    }
    fin.close();
}
// }}}
///



/*
 * class BioObj
 */
// {{{
BioObj::BioObj()
{
    data_m = NULL;
    flag_m = 0;
}

BioObj::BioObj(const BioObj & obj)
{
    data_m = obj.data_m;
    flag_m = obj.flag_m;
}

BioObj::~BioObj()
{
}

bool BioObj::is_bit_set(size_t bit)
{
    return (bool)(flag_m & (0x1<<bit));
}

bool BioObj::is_flag_set(uint64_t flag)
{
	return (bool)(flag_m & flag);
}

bool BioObj::is_any_flag_set(uint64_t flag)
{
    return (bool)(flag_m & flag);
}

bool BioObj::is_all_flag_set(uint64_t flag)
{
    return ((flag_m & flag) == flag);
}

void * BioObj::get_data()
{
    return data_m;
}

void BioObj::set(const BioObj & obj)
{
    flag_m = obj.flag_m;
    data_m = obj.data_m;
}

void BioObj::set(BioObj & obj)
{
    flag_m = obj.flag_m;
    data_m = obj.data_m;
}

void BioObj::set(BioObj * obj)
{
    if (obj!=NULL) {
        flag_m = obj->flag_m;
        data_m = obj->data_m;
    }
}

int BioObj::set_bit(size_t bit)
{
    if (bit < flag_size_m) {
        flag_m |= (0x1<<bit);
        return 0;
    } else {
        return 1;
    }
}

void BioObj::set_flag(uint64_t flag)
{
    flag_m |= flag;
}

void BioObj::set_flag()
{
    flag_m = 0xffffffffffffffff;
}

int BioObj::unset_bit(size_t bit)
{
    if (bit < flag_size_m) {
        flag_m &= ~(0x1<<bit);
        return 0;
    } else {
        return 1;
    }
}

void BioObj::unset_flag(uint64_t flag)
{
    flag_m &= ~flag;
}

void BioObj::unset_flag()
{
    flag_m = 0x0000000000000000;
}

void BioObj::set_data(void * data)
{
    data_m = data;
}

BioObj & BioObj::copy(const BioObj & obj)
{
    data_m = obj.data_m;
    flag_m = obj.flag_m;
    return (*this);
}

BioObj & BioObj::operator =(BioObj & obj)
{
    data_m = obj.data_m;
    flag_m = obj.flag_m;
    return (*this);
}
// }}}

};
