#include "czl_bio_annot.h"
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"

namespace czl_bio {

/**
 * Class Gene : inherit from SeqObj
 */
// {{{
Gene::Gene(): canon_tscript_m(NULL)
{
}

Gene::~Gene()
{
}

Chrom* Gene::get_chrom()
{
    return static_cast<Chrom*>(get_seq_pt());
}

Tscript* Gene::get_tscript(Int i)
{
    if (i<0 || i>=tscript_m.size()) {
        stringstream ss;
        ss << "In Gene::get_tscript, 1st parameter " << i << " NOT in [0," << tscript_m.size() << ")";
        Msg::error(ss);
    }
    return tscript_m[i];
}

Tscript* Gene::get_canon_tscript()
{
    return canon_tscript_m;
}

Int Gene::get_tscript_num(short is_include_null)
{
    if (is_include_null) {
        return tscript_m.size();
    } else {
        Int n=0;
        for (vector<Tscript*>::iterator it=tscript_m.begin(); it!=tscript_m.end(); ++it) {
            if (*it != NULL) n++;
        }
        return n;
    }
}

void Gene::set_chrom(Chrom* chrom)
{
    set_seq_pt(chrom);
}

void Gene::set_tscript(Int i, Tscript* tscript)
{
    if (i<0 || i>=tscript_m.size()) {
        stringstream ss;
        ss << "In Gene::set_tscript, 1st parameter " << i << " NOT in [0," << tscript_m.size() << ")";
        Msg::error(ss);
    }
    tscript_m[i] = tscript;
}

void Gene::set_canon_tscript(Tscript* t)
{
    canon_tscript_m = t;
}

void Gene::ins_tscript(Tscript* t)
{
    tscript_m.push_back(t);
}

/// Set the transcript as NULL by index
Tscript * Gene::rm_tscript(Int i)
{
    Tscript *p = NULL;
    if (i<0 || i>=tscript_m.size()) {
        stringstream ss;
        ss << "In Gene::rm_tscript, 1st parameter " << i << " NOT in [0," << tscript_m.size() << ")";
        Msg::error(ss);
    } else {
        p = tscript_m[i];
        tscript_m[i] = NULL;
    }
    return p;
}

// void Gene::find_tscript(Tscript* t)

void Gene::compact_vec()
{
    vector<Tscript*>::iterator it0=tscript_m.begin();
    Int n=0;
    for (vector<Tscript*>::iterator it=tscript_m.begin(); it!=tscript_m.end(); ++it) {
        if (*it != NULL) {
            if (it0!=it) {
                (*it0) = (*it);
                *it = NULL;
            }
            it0++;
            n++;
        }
    }
    tscript_m.resize(n);
}

void Gene::sort_tscript( bool(*cmp_fh)(Tscript*, Tscript*) )
{
    sort(tscript_m.begin(), tscript_m.end(), cmp_fh);
}

bool less_gene_pos(Gene *t1, Gene *t2)
{
    if (t1->get_seq_pt()->get_id() < t2->get_seq_pt()->get_id()) return true;
    else if (t1->get_seq_pt()->get_id() > t2->get_seq_pt()->get_id()) return false;

    if (t1->get_begin() < t2->get_begin()) return true;
    if (t1->get_begin() > t2->get_begin()) return false;

    if (t1->get_end() < t2->get_end()) return true;
    if (t1->get_end() > t2->get_end()) return false;

    if (t1->get_strand() < t2->get_strand()) return true;
    if (t1->get_strand() > t2->get_strand()) return false;

	return false;
}

bool less_gene_pos(Gene &t1, Gene &t2)
{
    if (t1.get_seq_pt()->get_id() < t2.get_seq_pt()->get_id()) return true;
    else if (t1.get_seq_pt()->get_id() > t2.get_seq_pt()->get_id()) return false;

    if (t1.get_begin() < t2.get_begin()) return true;
    if (t1.get_begin() > t2.get_begin()) return false;

    if (t1.get_end() < t2.get_end()) return true;
    if (t1.get_end() > t2.get_end()) return false;

    if (t1.get_strand() < t2.get_strand()) return true;
    if (t1.get_strand() > t2.get_strand()) return false;

	return false;
}

bool less_gene_chromsymbol_pos(Gene *t1, Gene *t2)
{
    if (t1->get_seq_pt()->get_symbol() < t2->get_seq_pt()->get_symbol()) return true;
    else if (t1->get_seq_pt()->get_symbol() > t2->get_seq_pt()->get_symbol()) return false;

    if (t1->get_begin() < t2->get_begin()) return true;
    if (t1->get_begin() > t2->get_begin()) return false;

    if (t1->get_end() < t2->get_end()) return true;
    if (t1->get_end() > t2->get_end()) return false;

    if (t1->get_strand() < t2->get_strand()) return true;
    if (t1->get_strand() > t2->get_strand()) return false;

	return false;
}


// }}}

/**
 * Class Tscript: inherit from BioObj
 */
// {{{
Tscript::Tscript(): canon_pep_m(NULL)
{
    is_et_sorted_m=false;
    n__null_et_m=0;
    n__null_pep_m=0;
}

Tscript::~Tscript()
{
}

Chrom* Tscript::get_chrom()
{
    return chrom_m;
}

Gene* Tscript::get_gene()
{
    return gene_m;
}

Pep* Tscript::get_pep(Int i)
{
    if (i<0 || i>=pep_m.size()) {
        stringstream ss;
        ss << "In Tscript::get_pep, 1st parameter " << i << " NOT in [0," << pep_m.size() << ")";
        Msg::error(ss);
    }
    return pep_m[i];
}

Pep* Tscript::get_canon_pep()
{
    return canon_pep_m;
}

Exon* Tscript::get_exon(Int i)
{
    if (i<0 || i>=et_m.size()) {
        stringstream ss;
        ss << "In Tscript::get_exon, 1st parameter " << i << " NOT in [0," << et_m.size() << ")";
        Msg::error(ss);
    }
    if (et_m[i]==NULL) return NULL;
    else return et_m[i]->get_exon();
}

ExonTscriptPair* Tscript::get_et(Int i)
{
    if (i<0 || i>=et_m.size()) {
        stringstream ss;
        ss << "In Tscript::get_et, 1st parameter " << i << " NOT in [0," << et_m.size() << ")";
        Msg::error(ss);
    }
    return et_m[i];
}

Int Tscript::get_pep_num(short is_include_null)
{
    if (is_include_null) {
        return pep_m.size();
    } else {
        Int n=0;
        for (vector<Pep*>::iterator it=pep_m.begin(); it!=pep_m.end(); ++it) {
            if (*it != NULL) n++;
        }
        return n;
    }
}

Int Tscript::get_et_num(short is_include_null)
{
    if (is_include_null) {
        return et_m.size();
    } else {
        Int n=0;
        for (vector<ExonTscriptPair*>::iterator it=et_m.begin(); it!=et_m.end(); ++it) {
            if (*it != NULL) n++;
        }
        return n;
    }
}

Int Tscript::get_exon_num(short is_include_null)
{
    return get_et_num(is_include_null);
}

Int Tscript::get_begin()
{
    if ( et_m.size()<=0 ) return -1;
    if (n__null_et_m>0) { compact_vec(); }
    if (!is_et_sorted_m) {
        sort(et_m.begin(), et_m.end(), less_et_exon_pos);
        is_et_sorted_m=true;
    }
    return et_m.front()->get_exon()->get_begin();
}

Int Tscript::get_end()
{
    if ( et_m.size()<=0 ) return -1;
    if (n__null_et_m>0) { compact_vec(); }
    if (!is_et_sorted_m) {
        sort(et_m.begin(), et_m.end(), less_et_exon_pos);
        is_et_sorted_m=true;
    }
    return et_m.back()->get_exon()->get_end();
}

Int Tscript::get_span_len()
{
    if ( et_m.size()<=0 ) return -1;
    if (n__null_et_m>0) { compact_vec(); }
    if (!is_et_sorted_m) {
        sort(et_m.begin(), et_m.end(), less_et_exon_pos);
        is_et_sorted_m=true;
    }
    return et_m.back()->get_exon()->get_end() - et_m.front()->get_exon()->get_begin();
}

Int Tscript::get_seq_len()
{
    Int len = 0;
    if ( et_m.size()<=0 ) return 0;
    if (n__null_et_m>0) { compact_vec(); }
    for (Int i=0; i<et_m.size(); i++) {
        len += et_m[i]->get_exon()->get_len();
    }
    return len;
}

void Tscript::set_chrom(Chrom* chrom)
{
    chrom_m = chrom;
}

void Tscript::set_gene(Gene* gene)
{
    gene_m = gene;
}

void Tscript::set_pep(Int i, Pep* pep)
{
    if (i<0 || i>=pep_m.size()) {
        stringstream ss;
        ss << "In Tscript::set_pep " << i << " not in [0," << pep_m.size() << ")";
        Msg::error(ss);
    }
    if (pep==NULL) {
        stringstream ss;
        ss << "In Tscript::set_pep, 2nd PARAM pep should not be 'NULL'";
        Msg::error(ss);
    }
    pep_m[i] = pep;
    if (pep_m[i]==NULL && pep!=NULL) n__null_pep_m++;
}

void Tscript::set_canon_pep(Pep * p)
{
    canon_pep_m = p;
}

void Tscript::set_et(Int i, ExonTscriptPair* et)
{
    if (i<0 || i>=et_m.size()) {
        stringstream ss;
        ss << "In Tscript::set_cp " << i << " not in [0," << et_m.size() << ")";
        Msg::error(ss);
    //    throw bpp::OutOfRangeException("In Tscript::set_cp", i, 0, et_m.size());
    }
    if (et==NULL) {
        stringstream ss;
        ss << "In Tscript::set_et, 2nd PARAM et should not be 'NULL'";
        Msg::error(ss);
    }
    if (et_m[i]==NULL && et!=NULL) n__null_et_m++;
    et_m[i] = et;
}

void Tscript::ins_pep(Pep* pep)
{
    pep_m.push_back(pep);
}

void Tscript::ins_et(ExonTscriptPair* et)
{
    et_m.push_back(et);
}

Pep* Tscript::rm_pep(Int i)
{
    Pep *p = NULL;
    if (i<0 || i>=pep_m.size()) {
        return NULL;
    } else {
        p = pep_m[i];
        pep_m[i]=NULL;
    }
    if (p!=NULL) { n__null_pep_m++; }
    return p;
}

ExonTscriptPair* Tscript::rm_et(Int i)
{
    ExonTscriptPair *p = NULL;
    if (i<0 || i>=et_m.size()) {
        return NULL;
    } else {
        p = et_m[i];
        et_m[i]=NULL;
    }
    if (p!=NULL) { n__null_et_m++; }
    return p;
}

void Tscript::compact_vec()
{
    Int n=0;
    if (n__null_pep_m>0) { 
        vector<Pep*>::iterator it1=pep_m.begin();
        for (vector<Pep*>::iterator it=pep_m.begin(); it!=pep_m.end(); ++it) {
            if (*it != NULL) {
                if (it1!=it) {
                    (*it1) = (*it);
                    *it = NULL;
                }
                n++;
                it1++;
            }
        }
        pep_m.resize(n);
        n__null_pep_m=0;
    }

    if (n__null_et_m>0) { 
        n=0;
        vector<ExonTscriptPair*>::iterator it2=et_m.begin();
        for (vector<ExonTscriptPair*>::iterator it=et_m.begin(); it!=et_m.end(); ++it) {
            if (*it != NULL) {
                if (it2!=it) {
                    (*it2) = (*it);
                    *it = NULL;
                }
                it2++;
                n++;
            }
        }
        et_m.resize(n);
        n__null_et_m=0;
    }
}

void Tscript::sort_et() 
{
    if ( n__null_et_m>0 ) compact_vec();
    sort(et_m.begin(), et_m.end(), less_et_exon_pos );
    is_et_sorted_m=true;
}

void Tscript::sort_et( bool(*cmp_fh)(ExonTscriptPair*, ExonTscriptPair*) )
{
    if ( n__null_et_m>0 ) compact_vec();
    sort(et_m.begin(), et_m.end(), cmp_fh);
}

bool less_tscript_pos(Tscript *t1, Tscript *t2)
{
    if (t1->get_chrom()->get_id() < t2->get_chrom()->get_id()) return true;
    else if (t1->get_chrom()->get_id() > t2->get_chrom()->get_id()) return false;

    if (t1->get_begin() < t2->get_begin()) return true;
    if (t1->get_begin() > t2->get_begin()) return false;

    if (t1->get_end() < t2->get_end()) {
        return true;
    } else {
        return false;
    }
}

bool less_tscript_pos(Tscript &t1, Tscript &t2)
{
    return less_tscript_pos(&t1, &t2);
}
// }}}

/**
 * Class Pep: inherit from BioObj
 */
// {{{
Chrom* Pep::get_chrom()
{
    return tscript_m->get_chrom();
}

short Pep::get_strand()
{
    return tscript_m->get_strand();
}

CDS* Pep::get_cds(Int i)
{
    if (i<0 || i>=cp_m.size()) {
        stringstream ss;
        ss << "In Pep::get_cds, 1st parameter " << i << " NOT in [0," << cp_m.size() << ")";
        Msg::error(ss);
    }
    if (cp_m[i] == NULL) return NULL;
    else return cp_m[i]->get_cds();
}

CDSPepPair* Pep::get_cp(Int i)
{
    if (i<0 || i>=cp_m.size()) {
        stringstream ss;
        ss << "In Pep::get_cp, 1st parameter " << i << " NOT in [0," << cp_m.size() << ")";
        Msg::error(ss);
    }
    return cp_m[i];
}

Int Pep::get_cp_num(short is_include_null)
{
    if (is_include_null) {
        return cp_m.size();
    } else {
        Int n=0;
        for (vector<CDSPepPair*>::iterator it=cp_m.begin(); it!=cp_m.end(); ++it) {
            if (*it != NULL) n++;
        }
        return n;
    }
}

Int Pep::get_cds_num(short is_include_null)
{
    return get_cp_num(is_include_null);
}

Int Pep::get_begin()
{
    if ( cp_m.size()<=0 ) return -1;
    if (n__null_cp_m>0) { compact_vec(); }
    if (!is_cp_sorted_m) {
        sort(cp_m.begin(), cp_m.end(), less_cp_cds_pos);
        is_cp_sorted_m=true;
    }
    return cp_m.front()->get_cds()->get_begin();
}

Int Pep::get_end()
{
    if ( cp_m.size()<=0 ) return -1;
    if (n__null_cp_m>0) { compact_vec(); }
    if (!is_cp_sorted_m) {
        sort(cp_m.begin(), cp_m.end(), less_cp_cds_pos);
        is_cp_sorted_m=true;
    }
    return cp_m.back()->get_cds()->get_end();
}

Int Pep::get_nt_len()
{
    return nt_len_m;
}

Int Pep::get_aa_len()
{
    return get_nt_len()/3;
}

void Pep::set_chrom(Chrom* chrom)
{
    chrom_m = chrom;
}

void Pep::set_tscript(Tscript* tscript)
{
    tscript_m = tscript;
}

void Pep::set_cp(Int i, CDSPepPair* cp)
{
    if (i<0 || i>=cp_m.size()) {
        stringstream ss;
        ss << "In Pep::rm_cp " << i << " is not in [0," << cp_m.size() << ")";
        Msg::error(ss);
    }
    if ( cp==NULL ) {
        stringstream ss;
        ss << "In Pep::rm_cp, 2nd parameter 'cp' should not be NULL";
        Msg::error(ss);
    }
    if ( cp_m[i]==NULL ) n__null_cp_m--;
    cp_m[i] = cp;
}

void Pep::set_nt_len(Int nt_len)
{
    nt_len_m = nt_len;
}

void Pep::ins_cp(CDSPepPair* cp)
{
    cp_m.push_back(cp);
}

void Pep::cal_len_from_cp()
{
    nt_len_m = 0;
    if ( n__null_cp_m > 0 ) compact_vec();
    if ( cp_m.size()<=0 ) return;
    for (Int i=0; i<cp_m.size(); i++) {
        nt_len_m += cp_m[i]->get_cds()->get_len();
    }
}

CDSPepPair* Pep::rm_cp(Int i)
{
    CDSPepPair *p = NULL;
    if (i<0 || i>=cp_m.size()) {
        stringstream ss;
        ss << "In Pep::rm_cp " << i << " is not in [0," << cp_m.size() << ")";
        Msg::error(ss);
    } else {
        p = cp_m[i];
        if (p!=NULL) n__null_cp_m++;
        cp_m[i]=NULL;
    }
    return p;
}

void Pep::compact_vec()
{
    if (n__null_cp_m<=0) return;
    vector<CDSPepPair*>::iterator it1=cp_m.begin();
    Int i=0;
    for (vector<CDSPepPair*>::iterator it=cp_m.begin(); it!=cp_m.end(); ++it) {
        if (*it != NULL) {
            if (it1!=it) {
                (*it1) = (*it);
                *it = NULL;
            }
            it1++;
            i++;
        }
    }
    cp_m.resize(i);
    n__null_cp_m = 0;
}

void Pep::sort_cp()
{
    if (n__null_cp_m>0) compact_vec();
    sort(cp_m.begin(), cp_m.end(), less_cp_cds_pos);
    is_cp_sorted_m = true;
}

void Pep::sort_cp( bool(*cmp_fh)(CDSPepPair*, CDSPepPair*) )
{
    if (n__null_cp_m>0) compact_vec();
    sort(cp_m.begin(), cp_m.end(), cmp_fh);
}

// }}}

/**
 * Class Exon: inherit from SeqRegion
 */
// {{{
ExonTscriptPair* Exon::get_et(Int i)
{
    if (i<0 || i>=et_m.size()) {
        stringstream ss;
        ss << "In Exon::get_et, 1st parameter " << i << "NOT in [0," << et_m.size() << ")";
        Msg::error(ss);
    }
    return et_m[i];
}

Tscript* Exon::get_tscript(Int i)
{
    ExonTscriptPair *et = get_et(i);
    if (et==NULL) return NULL;
    else return et->get_tscript();
}

Int Exon::get_et_num(short is_include_null)
{
    if (is_include_null) {
        return et_m.size();
    } else {
        Int n=0;
        for (vector<ExonTscriptPair*>::iterator it=et_m.begin(); it!=et_m.end(); ++it) {
            if (*it != NULL) n++;
        }
        return n;
    }
}

Int Exon::get_tscript_num(short is_include_null)
{
    return get_et_num(is_include_null);
}

CDS* Exon::get_cds(Int i)
{
    if (i<0 || i>=cds_m.size()) {
        stringstream ss;
        ss << "In Exon::get_cds, 1st parameter " << i << "NOT in [0," << cds_m.size() << ")";
        Msg::error(ss);
    }
    return cds_m[i];
}

Int Exon::get_cds_num(short is_include_null)
{
    if (is_include_null) {
        return cds_m.size();
    } else {
        Int n=0;
        for (vector<CDS*>::iterator it=cds_m.begin(); it!=cds_m.end(); ++it) {
            if (*it != NULL) n++;
        }
        return n;
    }
}

Chrom * Exon::get_chrom()
{
    return (Chrom*)get_seq_pt();
}

void Exon::set_et(Int i, ExonTscriptPair* et)
{
    if (i<0 || i>=et_m.size()) {
        stringstream ss;
        ss << "In Exon::set_et, 1st parameter " << i << "NOT in [0," << et_m.size() << ")";
        Msg::error(ss);
    }
    et_m[i] = et;
}

void Exon::set_cds(Int i, CDS* cds)
{
    if (i<0 || i>=cds_m.size()) {
        stringstream ss;
        ss << "In Exon::set_cds, 1st parameter " << i << "NOT in [0," << cds_m.size() << ")";
        Msg::error(ss);
    }
    cds_m[i] = cds;
}

void Exon::ins_et(ExonTscriptPair* et)
{
    et_m.push_back(et);
}

void Exon::ins_cds(CDS* cds)
{
    cds_m.push_back(cds);
}

ExonTscriptPair* Exon::rm_et(Int i)
{
    ExonTscriptPair *p = NULL;
    if (i<0 || i>=et_m.size()) {
        stringstream ss;
        ss << "In Exon::rm_et, 1st parameter " << i << "NOT in [0," << et_m.size() << ")";
        Msg::error(ss);
    } else {
        p = et_m[i];
        et_m[i]=NULL;
    }
    return p;
}

CDS* Exon::rm_cds(Int i)
{
    CDS *p = NULL;
    if (i<0 || i>=cds_m.size()) {
        stringstream ss;
        ss << "In Exon::rm_cds, 1st parameter " << i << "NOT in [0," << cds_m.size() << ")";
        Msg::error(ss);
    } else {
        p = cds_m[i];
        cds_m[i]=NULL;
    }
    return p;
}

void Exon::compact_vec()
{
    vector<ExonTscriptPair*>::iterator it0=et_m.begin();
    Int i=0;
    for (vector<ExonTscriptPair*>::iterator it=et_m.begin(); it!=et_m.end(); ++it) {
        if (*it != NULL) {
            if (it0!=it) {
                (*it0) = (*it);
                *it = NULL;
            }
            it0++;
            i++;
        }
    }
    et_m.resize(i);

    i=0;
    vector<CDS*>::iterator it1=cds_m.begin();
    for (vector<CDS*>::iterator it=cds_m.begin(); it!=cds_m.end(); ++it) {
        if (*it != NULL) {
            if (it1!=it) {
                (*it1) = (*it);
                *it = NULL;
            }
            it1++;
            i++;
        }
    }
    cds_m.resize(i);
}

bool less_exon_pos(Exon *t1, Exon *t2)
{
    if (t1->get_seq_pt()->get_id() < t2->get_seq_pt()->get_id()) return true;
    else if (t1->get_seq_pt()->get_id() > t2->get_seq_pt()->get_id()) return false;

    if (t1->get_begin() < t2->get_begin()) return true;
    if (t1->get_begin() > t2->get_begin()) return false;

    if (t1->get_end() < t2->get_end()) {
        return true;
    } else {
        return false;
    }
}
// }}}

/**
 * Class CDS: inherit from SeqRegion
 */
// {{{
CDSPepPair* CDS::get_cp(Int i)
{
    if (i<0 || i>=cp_m.size()) {
        stringstream ss;
        ss << "In In CDS::get_cp, 1st parameter " << i << "NOT in [0," << cp_m.size() << ")";
        Msg::error(ss);
    }
    return cp_m[i];
}

Pep* CDS::get_pep(Int i)
{
    CDSPepPair *cp = get_cp(i);
    if (cp==NULL) return NULL;
    else return cp->get_pep();
}

Int CDS::get_cp_num(short is_include_null)
{
    if (is_include_null) {
        return cp_m.size();
    } else {
        Int n=0;
        for (vector<CDSPepPair*>::iterator it=cp_m.begin(); it!=cp_m.end(); ++it) {
            if (*it != NULL) n++;
        }
        return n;
    }
}

Int CDS::get_pep_num(short is_include_null)
{
    return get_cp_num(is_include_null);
}

//void CDS::set_cp(Int i, CDSPepPair* cp)
//{
//    if (i<0 || i>=cp_m.size()) {
//        throw bpp::OutOfRangeException("In CDS::set_cp", i, 0, cp_m.size());
//    }
//    cp_m[i] = cp;
//}
Chrom * CDS::get_chrom()
{
    return (Chrom*)get_seq_pt();
}


void CDS::ins_cp(CDSPepPair* cp)
{
    cp_m.push_back(cp);
}

CDSPepPair* CDS::rm_cp(Int i)
{
    CDSPepPair *p = NULL;
    if (i<0 || i>=cp_m.size()) {
        stringstream ss;
        ss << "In CDSPepPair::rm_cp, 1st parameter " << i << "NOT in [0," << cp_m.size() << ")";
        Msg::error(ss);
    } else {
        p = cp_m[i];
        cp_m[i]=NULL;
    }
    return p;
}

void CDS::compact_vec()
{
    Int n=0, i;
    for (i=0; i<cp_m.size(); i++) {
        if (cp_m[i]!=NULL) {
            if (i!=n) {
                cp_m[n] = cp_m[i];
                n++;
            }
        }
    }
    cp_m.resize(n);
}

bool less_cds_pos(CDS *t1, CDS *t2)
{
    if (t1->get_seq_pt()->get_id() < t2->get_seq_pt()->get_id()) return true;
    else if (t1->get_seq_pt()->get_id() > t2->get_seq_pt()->get_id()) return false;

    if (t1->get_begin() < t2->get_begin()) return true;
    if (t1->get_begin() > t2->get_begin()) return false;

    if (t1->get_end() < t2->get_end()) {
        return true;
    } else {
        return false;
    }
}
// }}}

/**
 * Class ExonTscriptPair
 */
// {{{
/**
 * @brief   less function for sorting by CDS position
 * @return  true if (begin1<begin2 OR begin1=begin2 AND end1<end2) 
 */
bool less_et_exon_pos(ExonTscriptPair *et1, ExonTscriptPair *et2)
{
    if (et1->get_exon()->get_begin() < et2->get_exon()->get_begin() ) {
        return true;
    } else if (et1->get_exon()->get_begin() == et2->get_exon()->get_begin() && et1->get_exon()->get_end() < et2->get_exon()->get_end() ) {
        return true;
    } else return false;
}

/**
 * @brief   less function for sorting by transcript position
 * @return  true if (begin1<begin2 OR begin1=begin2 AND end1<end2) 
 * @warning run compact_vec before sorting
 */
bool less_et_tscript_pos(ExonTscriptPair *et1, ExonTscriptPair *et2)
{
    if (et1->get_tscript()->get_begin() < et2->get_tscript()->get_begin() ) {
        return true;
    } else if (et1->get_tscript()->get_begin() == et2->get_tscript()->get_begin() && et1->get_tscript()->get_end() < et2->get_tscript()->get_end() ) {
        return true;
    } else return false;
}

// }}}

/**
 * Class CDSPepPair
 */
// {{{
/**
 * @brief   less function for sorting by CDS position
 * @return  true if (begin1<begin2 OR begin1=begin2 AND end1<end2) 
 */
bool less_cp_cds_pos(CDSPepPair *cp1, CDSPepPair *cp2)
{
    if (cp1->get_cds()->get_begin() < cp2->get_cds()->get_begin() ) {
        return true;
    } else if (cp1->get_cds()->get_begin() == cp2->get_cds()->get_begin() && cp1->get_cds()->get_end() < cp2->get_cds()->get_end() ) {
        return true;
    } else return false;
}

/**
 * @brief   less function for sorting by transcript position
 * @return  true if (begin1<begin2 OR begin1=begin2 AND end1<end2) 
 * @warning run compact_vec before sorting
 */
bool less_cp_pep_pos(CDSPepPair *cp1, CDSPepPair *cp2)
{
    if (cp1->get_pep()->get_begin() < cp2->get_pep()->get_begin() ) {
        return true;
    } else if (cp1->get_pep()->get_begin() == cp2->get_pep()->get_begin() && cp1->get_pep()->get_end() < cp2->get_pep()->get_end() ) {
        return true;
    } else return false;
}

// }}}

/**
 * Class GTP: inherit from BioObj
 */
GTP::GTP()
{
}

GTP::~GTP()
{
}

vector<Gene*> & GTP::get_gene()
{
    return gene_m;
}

Gene* GTP::get_gene(Int i)
{
    if (i<0 || i>=gene_m.size()) {
        stringstream ss;
        ss << "In GTP::get_gene, 1st parameter " << i << "should in [0," << gene_m.size() << ")";
        Msg::error(ss);
        return NULL;
    } else {
        return gene_m[i];
    }
}

vector<Tscript*> & GTP::get_tscript()
{
    return tscript_m;
}

Tscript* GTP::get_tscript(Int i)
{
    if (i<0 || i>=tscript_m.size()) {
        stringstream ss;
        ss << "In GTP::get_tscript, 1st parameter " << i << "should in [0," << tscript_m.size() << ")";
        Msg::error(ss);
        return NULL;
    } else {
        return tscript_m[i];
    }
}

vector<Pep*> & GTP::get_pep()
{
    return pep_m;
}

Pep* GTP::get_pep(Int i)
{
    if (i<0 || i>=pep_m.size()) {
        stringstream ss;
        ss << "In GTP::get_pep, 1st parameter " << i << "should in [0," << pep_m.size() << ")";
        Msg::error(ss);
        return NULL;
    } else {
        return pep_m[i];
    }
}

vector<Exon*> & GTP::get_exon()
{
    return exon_m;
}

Exon* GTP::get_exon(Int i)
{
    if (i<0 || i>=exon_m.size()) {
        stringstream ss;
        ss << "In GTP::get_exon, 1st parameter " << i << "should in [0," << exon_m.size() << ")";
        Msg::error(ss);
        return NULL;
    } else {
        return exon_m[i];
    }
}

vector<CDS*> & GTP::get_cds()
{
    return cds_m;
}

CDS* GTP::get_cds(Int i)
{
    if (i<0 || i>=cds_m.size()) {
        stringstream ss;
        ss << "In GTP::get_cds, 1st parameter " << i << "should in [0," << cds_m.size() << ")";
        Msg::error(ss);
        return NULL;
    } else {
        return cds_m[i];
    }
}

//vector<Variant*> & GTP::get_variant()
//{
//    return variant_m;
//}

Int GTP::get_chrom_num()
{
    return chrom_m.size();
}

Int GTP::get_gene_num()
{
    return gene_m.size();
}

Int GTP::get_tscript_num()
{
    return tscript_m.size();
}

Int GTP::get_pep_num()
{
    return pep_m.size();
}

Int GTP::get_exon_num()
{
    return exon_m.size();
}

Int GTP::get_cds_num()
{
    return cds_m.size();
}


Chrom* GTP::get_chrom_by_symbol(string & symbol)
{
    if ( map_chrom_symbol_to_obj_m.find(symbol) == map_chrom_symbol_to_obj_m.end() ) {
        return NULL;
    } else {
        return map_chrom_symbol_to_obj_m[symbol];
    }
}

Chrom* GTP::get_chrom_by_id(string & id)
{
    if ( map_chrom_id_to_obj_m.find(id) == map_chrom_id_to_obj_m.end() ) {
        return NULL;
    } else {
        return map_chrom_id_to_obj_m[id];
    }
}

Gene* GTP::get_gene_by_id(string & id)
{
    if ( map_gene_id_to_obj_m.find(id) == map_gene_id_to_obj_m.end() ) {
        return NULL;
    } else {
        return map_gene_id_to_obj_m[id];
    }
}

Tscript* GTP::get_tscript_by_id(string & id)
{
    if ( map_tscript_id_to_obj_m.find(id) == map_tscript_id_to_obj_m.end() ) {
        return NULL;
    } else {
        return map_tscript_id_to_obj_m[id];
    }
}

Pep* GTP::get_pep_by_id(string & id)
{
    if ( map_pep_id_to_obj_m.find(id) == map_pep_id_to_obj_m.end() ) {
        return NULL;
    } else {
        return map_pep_id_to_obj_m[id];
    }
}

Variant* GTP::get_variant_by_id_range_mut(string& id, Int begin, Int end, string & mut_seq)
{
    stringstream ss;
    ss << id << "\t" << begin << "\t" << end << "\t" << mut_seq;
    string s = ss.str();
    if ( map_variant_to_obj_m.find(s) == map_variant_to_obj_m.end() ) {
        return NULL;
    } else {
        return map_variant_to_obj_m[s];
    }
}

void GTP::set_chrom_by_id(string & id, Chrom* chrom)
{
    map_chrom_id_to_obj_m[id] = chrom;
}

void GTP::set_gene_by_id(string & id, Gene* g)
{
    map_gene_id_to_obj_m[id] = g;
}

void GTP::set_tscript_by_id(string & id, Tscript* t)
{
    map_tscript_id_to_obj_m[id] = t;
}

void GTP::set_pep_by_id(string & id, Pep *p)
{
    map_pep_id_to_obj_m[id] = p;
}

void GTP::set_variant_by_id_range_mut(string& id, Int begin, Int end, string & mut_seq, Variant *v)
{
    stringstream ss;
    ss << id << "\t" << begin << "\t" << end << "\t" << mut_seq;
    string s = ss.str();
    map_variant_to_obj_m[s] = v;
}

void GTP::read_from_db(string dir)
{
    Int i, j, k, l, m, n;

    /// read genes from <dir>/gene.txt
    string chrom_file = dir + "/seq_region.txt";
    string gene_file = dir + "/gene.txt";
    string tscript_file = dir + "/tscript.txt";
    string pep_file = dir + "/pep.txt";
    string exon_file = dir + "/exon.txt";
    string cds_file = dir + "/cds.txt";
    string exon_to_tscript_file = dir + "/exon_to_tscript.txt";
    string cds_to_pep_file = dir + "/cds_to_pep.txt";
    string files[8] = {chrom_file, gene_file, tscript_file, pep_file, exon_file, cds_file, exon_to_tscript_file, cds_to_pep_file};
    stringstream ss;
    for (i=0; i<6; i++) {
        if ( !boost::filesystem::exists(files[i]) ) {
            ss << "Can't open FILE " << files[i];
            Msg::error(ss);
            return;
        }
    }

    string str;
    string value;
    ifstream fin;

    /// read chromosome/scaffolds
    // {{{
    fin.open(chrom_file.c_str());
    while ( !fin.eof() ) {
        getline(fin, str);
        if ( str.empty() ) { continue; }
        vector<string> tab;
        boost::split( tab, str, boost::is_any_of("\t") );
        BOOST_FOREACH(string & s, tab) {
            boost::trim(s);
        }
        /// id  symbol  coord_system_id  length
        Chrom * chrom = new Chrom;
        chrom->set_id(tab[0]);
        if (map_chrom_id_to_obj_m.find(chrom->get_id()) != map_chrom_id_to_obj_m.end()) {
            cout << "Warnning: Chrom " << chrom->get_id() << " exists." << endl;
            delete chrom;
            continue;
        }
        chrom->set_symbol(tab[1]);
        chrom->set_coord_system_id( atoi(tab[2].c_str()) );
//      chrom->length = atol(tab[3].c_str());

        chrom_m.push_back(chrom);
        map_chrom_id_to_obj_m[chrom->get_id()] = chrom;
        map_chrom_symbol_to_obj_m[chrom->get_symbol()] = chrom;
    }
    fin.close();
    // }}}

    /// read gene
    // {{{
    fin.open(gene_file.c_str());
    vector<string> gene_canon_tscript;
    while ( !fin.eof() ) {
        getline(fin, str);
        if ( str.empty() ) { continue; }
        vector<string> tab;
        boost::split( tab, str, boost::is_any_of("\t") );
        BOOST_FOREACH(string & s, tab) {
            boost::trim(s);
        }
        /// gene header in db file
        /// part 1, 1--5 : id  symbol  name  synonym  desc
        Gene* gene = new Gene;
        gene->set_id(tab[0]);
        if (map_gene_id_to_obj_m.find(gene->get_id()) != map_gene_id_to_obj_m.end()) {
            string msg = "Gene " + gene->get_id() + " exists.";
            Msg::warn(msg);
            delete gene;
            continue;
        }
        gene->set_symbol(tab[1]);

        /// 6--10 : chrom  start  stop  strand  seq_length
        Chrom* chrom_p;
        value = tab[5];
        if ( map_chrom_id_to_obj_m.find(value) != map_chrom_id_to_obj_m.end() ) {
            chrom_p = map_chrom_id_to_obj_m[value];
            gene->set_seq_pt(chrom_p);
        }
        Int begin = atoi(tab[6].c_str()) -1;
        Int end = atoi(tab[7].c_str()) ;
        gene->ins_range(begin, end);
        Int strand = atoi(tab[8].c_str()) ;
        gene->set_strand(strand);

        /// 11--E : cannonical_tscript  biostatus ...
        string canon_tscript = tab[10];
        gene_canon_tscript.push_back(canon_tscript);

        boost::to_upper(tab[11]);
        string biostatus = tab[11];
        if ( biostatus == "NOVEL" ) {
            gene->set_flag(NOVEL_BIT);
        }

        string version = tab[12]; ///< ignore version

        gene_m.push_back(gene);
        map_gene_id_to_obj_m[gene->get_id()] = gene;
    }
    fin.close();
    // }}}

    /// read tscript
    // {{{
    fin.open(tscript_file.c_str());
    vector<string> tscript_canon_pep;
    while ( !fin.eof() ) {
        getline(fin, str);
        if ( str.empty() ) { continue; }
        vector<string> tab;
        boost::split( tab, str, boost::is_any_of("\t") );
        BOOST_FOREACH(string & s, tab) {
            boost::trim(s);
        }
        Tscript* tscript = new Tscript;
        /// part 1, 1--5 : id  symbol  name  synonym  desc
        tscript->set_id(tab[0]);
        if (map_tscript_id_to_obj_m.find(tscript->get_id()) != map_tscript_id_to_obj_m.end()) {
            string msg = "Tscript " + tscript->get_id() + " exists.";
            Msg::warn(msg);
            delete tscript;
            continue;
        }
        tscript->set_symbol(tab[1]);
        ///

        /// part 2, 6--10: chrom  start  stop  strand  seq_len
        Chrom* chrom_p;
        if ( map_chrom_id_to_obj_m.find(tab[5]) != map_chrom_id_to_obj_m.end() ) {
            chrom_p = map_chrom_id_to_obj_m[tab[5]];
            tscript->set_chrom(chrom_p);
        }
        Int begin = atoi(tab[6].c_str()) -1;  ///< 1-base in file, 0-base in program, ignore
        Int end = atoi(tab[7].c_str()) ;  ///< ignore
        Int strand = atoi(tab[8].c_str()) ;
        tscript->set_strand(strand);
        Int seq_len = atoi(tab[9].c_str()) ;  ///< ignore
//      tscript->set_seq_len(seq_len);
        ///
        
        /// part 3: parent_gene_id  start_in_gene  exon_num  canon_pep_id
        string gene_id = tab[10];
        if ( map_gene_id_to_obj_m.find(gene_id) != map_gene_id_to_obj_m.end() ) {
            tscript->set_gene( map_gene_id_to_obj_m[gene_id] );
        }

        value = tab[11]; ///< start_in_gene, ignored!
        value = tab[12]; ///< exon_num, ignored!

        string canon_pep = tab[13];
        tscript_canon_pep.push_back(canon_pep);
        ///

        string biostatus = tab[14]; ///< ignore biostatus!
        string version = tab[15]; ///< ignore version!

        tscript_m.push_back(tscript);
        map_tscript_id_to_obj_m[tscript->get_id()] = tscript;
    }
    fin.close();
    // }}}

    /// read pep
    // {{{
    fin.open(pep_file.c_str());
    vector<string> pep_canon_pep;
    while ( !fin.eof() ) {
        getline(fin, str);
        if ( str.empty() ) { continue; }
        vector<string> tab;
        boost::split( tab, str, boost::is_any_of("\t") );
        BOOST_FOREACH(string & s, tab) {
            boost::trim(s);
        }
        /// part 1: id  symbol  name  synonym  desc
        Pep* pep = new Pep;
        pep->set_id(tab[0]);
        if (map_pep_id_to_obj_m.find(pep->get_id()) != map_pep_id_to_obj_m.end()) {
            string msg = "Pep " + pep->get_id() + " exists.";
            Msg::warn(msg);
            delete pep;
            continue;
        }
        pep->set_symbol(tab[1]);
        string name    = tab[2];  ///< ignore
        string synonym = tab[3];  ///< ignore
        string desc    = tab[4];  ///< ignore
        ///

        /// part 2: start  stop  seq_len(aa_len)
        Int begin = atol(tab[5].c_str()) -1;  ///< ignore
        Int end = atol(tab[6].c_str()) ;  ///< ignore
        Int aa_seq_len = atol(tab[7].c_str());  ///< ignore aa_seq_len
        ///

        /// part 3: canon_transcript  start_in_tscript  cds_num  cds_len
        string tscript_id = tab[8];
        if ( map_tscript_id_to_obj_m.find(tscript_id) != map_tscript_id_to_obj_m.end() ) {
            Tscript *p = map_tscript_id_to_obj_m[tscript_id];
            pep->set_tscript( p );
            pep->set_chrom(p->get_chrom());
        }

        value = tab[9]; ///< start_in_tscript, ignored!

        value = tab[10]; ///< cds_num, ignored!

        Int nt_len = atoi(tab[11].c_str()) ; ///< cds_len, ignored!
        pep->set_nt_len(nt_len);
        ///

        string version = tab[12]; ///< ignored!

        pep_m.push_back(pep);
        map_pep_id_to_obj_m[pep->get_id()] = pep;
    }
    fin.close();
    // }}}

    /// read exon
    // {{{
    fin.open(exon_file.c_str());
    map<string, Exon*> exon_id_to_obj;
    while ( !fin.eof() ) {
        getline(fin, str);
        if ( str.empty() ) { continue; }
        vector<string> tab;
        boost::split( tab, str, boost::is_any_of("\t") );
        BOOST_FOREACH(string & s, tab) {
            boost::trim(s);
        }
        /// id  start  stop  seq_len  is_constitutive
        Exon* exon = new Exon;
        exon->set_id(tab[0]);
        if (exon_id_to_obj.find(exon->get_id()) != exon_id_to_obj.end()) {
            string msg = "Exon " + exon->get_id() + " exists.";
            Msg::warn(msg);
            delete exon;
            continue;
        }

        Int begin = atoi(tab[1].c_str()) -1;
        Int end   = atoi(tab[2].c_str()) ;
        exon->ins_range(begin, end);

        Int seq_len = atoi(tab[3].c_str()); ///< ignore seq_len

        Int is_constitutive = atoi(tab[4].c_str());
        if (is_constitutive) exon->set_flag(CONSTITUTIVE_BIT);

        exon_m.push_back(exon);
        exon_id_to_obj[exon->get_id()] = exon;
    }
    fin.close();
    // }}}

    /// read cds
    // {{{
    fin.open(cds_file.c_str());
    map<string, CDS*> cds_id_to_obj;
    while ( !fin.eof() ) {
        getline(fin, str);
        if ( str.empty() ) { continue; }
        vector<string> tab;
        boost::split( tab, str, boost::is_any_of("\t") );
        BOOST_FOREACH(string & s, tab) {
            boost::trim(s);
        }
        /// id  start  stop  seq_len  exon_id
        CDS* cds = new CDS;
        cds->set_id(tab[0]);
        if (cds_id_to_obj.find(cds->get_id()) != cds_id_to_obj.end()) {
            string msg = "CDS " + cds->get_id() + " exists.";
            Msg::warn(msg);
            delete cds;
            continue;
        }

        Int begin = atoi(tab[1].c_str()) -1;
        Int end   = atoi(tab[2].c_str()) ;
        cds->ins_range(begin, end);

        Int seq_len = atoi(tab[3].c_str());  ///< ignore seeq_len

        string exon_id = tab[4];
        if ( exon_id_to_obj.find(exon_id) != exon_id_to_obj.end() ) {
            Exon *p = exon_id_to_obj[exon_id];
            cds->set_exon( p );
        }

        cds_m.push_back(cds);
        cds_id_to_obj[cds->get_id()] = cds;
    }
    fin.close();
    // }}}

    /// set gene.tscript
    // {{{
    for (i=0; i<tscript_m.size(); i++) {
        tscript_m[i]->get_gene()->ins_tscript(tscript_m[i]);
    }
    // }}}

    /// set tscript.pep
    // {{{
    for (i=0; i<pep_m.size(); i++) {
        pep_m[i]->get_tscript()->ins_pep(pep_m[i]);
    }
    // }}}

    /// set exon.cds
    // {{{
    for (i=0; i<cds_m.size(); i++) {
        cds_m[i]->get_exon()->ins_cds(cds_m[i]);
    }
    // }}}
    
    /// read exon_to_tscript
    // {{{
    fin.open(exon_to_tscript_file.c_str());
    while ( !fin.eof() ) {
        getline(fin, str);
        if ( str.empty() ) { continue; }
        vector<string> tab;
        boost::split( tab, str, boost::is_any_of("\t") );
        BOOST_FOREACH(string & s, tab) {
            boost::trim(s);
        }
        /// exon_id  tscript_id  start_in_tscript  
        Exon* exon;
        string exon_id = tab[0];
        if (exon_id_to_obj.find(exon_id) == exon_id_to_obj.end()) {
            string msg = "Exon " + exon_id + "not exists.";
            Msg::warn(msg);
            continue;
        } else {
            exon = exon_id_to_obj[exon_id];
        }

        Tscript* tscript;
        string tscript_id = tab[1];
        if (map_tscript_id_to_obj_m.find(tscript_id) == map_tscript_id_to_obj_m.end()) {
            string msg = "Tscript " + tscript_id + "not exists.";
            Msg::warn(msg);
            continue;
        } else {
            tscript = map_tscript_id_to_obj_m[tscript_id];
        }
        exon->set_seq_pt(tscript->get_chrom());

        Int begin_in_tscript = atoi(tab[2].c_str()) - 1;

        Int rank = atoi(tab[3].c_str()) - 1;

        Int start_phase = atoi(tab[4].c_str());
        Int stop_phase = atoi(tab[5].c_str());

        ExonTscriptPair *et = new ExonTscriptPair(exon, tscript, begin_in_tscript, rank, start_phase, stop_phase);

        et_m.push_back(et);
    }
    fin.close();
    // }}}
    ///

    /// set tscript.et, exon.et
    for (i=0; i<et_m.size(); i++) {
        et_m[i]->get_exon()->ins_et(et_m[i]);
        et_m[i]->get_tscript()->ins_et(et_m[i]);
    }
    ///

    /// read cds_to_pep
    // {{{
    fin.open(cds_to_pep_file.c_str());
    while ( !fin.eof() ) {
        getline(fin, str);
        if ( str.empty() ) { continue; }
        vector<string> tab;
        boost::split( tab, str, boost::is_any_of("\t") );
        BOOST_FOREACH(string & s, tab) {
            boost::trim(s);
        }
        /// cds_id  pep_id  start_in_pep  
        CDS* cds;
        string cds_id = tab[0];
        if (cds_id_to_obj.find(cds_id) == cds_id_to_obj.end()) {
            string msg = "CDS " + cds_id + "not exists.";
            Msg::warn(msg);
            continue;
        } else {
            cds = cds_id_to_obj[cds_id];
        }

        Pep* pep;
        string pep_id = tab[1];
        if (map_pep_id_to_obj_m.find(pep_id) == map_pep_id_to_obj_m.end()) {
            string msg = "Pep " + pep_id + "not exists.";
            Msg::warn(msg);
            continue;
        } else {
            pep = map_pep_id_to_obj_m[pep_id];
        }
        cds->set_seq_pt(pep->get_chrom());

        Int begin_in_tscript = atoi(tab[2].c_str()) - 1;

        Int rank = atoi(tab[3].c_str()) - 1;

        Int start_phase = atoi(tab[4].c_str());
        Int stop_phase = atoi(tab[5].c_str());

        CDSPepPair *cp = new CDSPepPair(cds, pep, begin_in_tscript, rank, start_phase, stop_phase);

        cp_m.push_back(cp);
    }
    fin.close();
    // }}}
    ///
    
    /// set pep.et, cds.et
    for (i=0; i<cp_m.size(); i++) {
        cp_m[i]->get_cds()->ins_cp(cp_m[i]);
        cp_m[i]->get_pep()->ins_cp(cp_m[i]);
    }
    ///

    /// set bioflag:CANONICAL, CODING 
    // {{{
    for (i=0; i<gene_m.size(); i++) {
        Tscript *ts = gene_m[i]->get_canon_tscript();
        if (ts!=NULL) {
            ts->set_flag(CANONICAL_BIT);
        }
    }
    for (i=0; i<tscript_m.size(); i++) {
        Pep *p = tscript_m[i]->get_canon_pep();
        if (p!=NULL) {
            p->set_flag(CANONICAL_BIT);
            tscript_m[i]->set_flag(CODING_BIT);
            Gene *g = tscript_m[i]->get_gene();
            if (g!=NULL) {
                g->set_flag(CODING_BIT);
            }
        }
    }
    // }}}
    ///
}

void GTP::sort_gene( bool(*cmp_fh)(Gene*, Gene*) )
{
    sort(gene_m.begin(), gene_m.end(), cmp_fh);
}

void GTP::sort_tscript( bool(*cmp_fh)(Tscript*, Tscript*) )
{
    sort(tscript_m.begin(), tscript_m.end(), cmp_fh);
}

void GTP::sort_exon( bool(*cmp_fh)(Exon*, Exon*) )
{
    sort(exon_m.begin(), exon_m.end(), cmp_fh);
}

void GTP::sort_cds( bool(*cmp_fh)(CDS*, CDS*) )
{
    sort(cds_m.begin(), cds_m.end(), cmp_fh);
}

Chrom* GTP::search_chrom(bool (*is_fh)(Chrom*, Chrom*), Chrom* c)
{
    vector<Chrom*> cv(1, c);
    vector<Chrom*>::iterator it = search(chrom_m.begin(), chrom_m.end(), cv.begin(), cv.end(), is_fh);
    if (it==chrom_m.end()) return NULL;
    return (*it);
}

Gene* GTP::search_gene(bool (*is_fh)(Gene*, Gene*), Gene* g)
{
    vector<Gene*> gv(1, g);
    vector<Gene*>::iterator it = search(gene_m.begin(), gene_m.end(), gv.begin(), gv.end(), is_fh);
    if (it==gene_m.end()) return NULL;
    return (*it);
}

void GTP::search_gene(bool (*is_fh)(Gene*, Gene*), Gene* g, vector<Gene*> & out_genes)
{
    vector<Gene*> gv(1, g);
    vector<Gene*>::iterator it = search(gene_m.begin(), gene_m.end(), gv.begin(), gv.end(), is_fh);
    if (it==gene_m.end()) return;

    vector<Gene*>::iterator it0 = it;
    out_genes.push_back((*it));
    ++it0;
    while (it0!=gene_m.end()) {
        Gene *pg = (*it0);
        if (g->get_chrom()->get_id()!=pg->get_chrom()->get_id()) break;
        if (g->get_begin() < pg->get_begin()) {
            break;
        } else {
            if (g->get_end() < pg->get_end()) {
                out_genes.push_back(pg);
            }
            ++it0;
        }
    }
}

Tscript* GTP::search_tscript(bool (*is_fh)(Tscript*, Tscript*), Tscript* t)
{
    vector<Tscript*> tv(1, t);
    vector<Tscript*>::iterator it = search(tscript_m.begin(), tscript_m.end(), tv.begin(), tv.end(), is_fh);
    if (it==tscript_m.end()) return NULL;
    return (*it);
}

Int GTP::get_codon_by_pos(string & chr_oid, Int pos, vector<string> & codons, vector<Int> & phases, vector<CDS*> & out_cdss, vector<Pep*> & out_peps, vector<Gene*> &out_genes)
{
    Chrom c;
    c.set_id(chr_oid);
    Chrom *chr_obj = search_chrom(is_same_chrom, &c);
    return get_codon_by_pos(chr_obj, pos, codons, phases, out_cdss, out_peps, out_genes);
}

Int GTP::get_codon_by_pos(Chrom *chr_obj, Int pos, vector<string> & codons, vector<Int> & phases, vector<CDS*> & out_cdss, vector<Pep*> & out_peps, vector<Gene*> &out_genes)
{
    Gene g;
    g.set_chrom(chr_obj);
    g.ins_range(pos, pos+1);
    search_gene(is_same_gene, &g, out_genes);
    if (out_genes.size()<=0) return 0;
    map<string, Int> map_codon_to_phase;
    Int i;

    for (Int gi=0; gi<out_genes.size(); gi++) {
        Gene *gene_obj = out_genes[gi];
        for (i=0; i<gene_obj->get_tscript_num(); i++) {
            Tscript *pt = gene_obj->get_tscript(i);
            for (Int i_p=0; i_p<pt->get_pep_num(); i_p++) {
                Pep *pp = pt->get_pep(i_p);
                for (Int i_cp=0; i_cp<pp->get_cp_num(); i_cp++) {
                    CDSPepPair *pcp = pp->get_cp(i_cp);
                    CDS * pc = pcp->get_cds();
                    Int b_phase = pcp->get_begin_phase();
                    Int e_phase = pcp->get_end_phase();
                    if ( pos>=pc->get_begin() && pos<pc->get_end()) {
                        string codon;
                        Int phase;
                        phase = ( (pos - pc->get_begin()) + b_phase )%3;
                        Int begin = pos - phase;
                        if ( begin < pc->get_begin() ) {
                            Int l = phase;
                            CDSPepPair *prev_pcp;
                            i_cp--;
                            while (i_cp>=0) { 
                                prev_pcp = pp->get_cp(i_cp);
                                if (prev_pcp!=NULL) break;
                                i_cp--;
                            }
                            if (i_cp>=0) {
                                CDS * prev_pc = prev_pcp->get_cds();
                                codon = chr_obj->get_seq(prev_pc->get_end()-l, l);
                                codon += chr_obj->get_seq(pc->get_begin(), 3-l);
                            }
                        } else if ( begin+3 >= pc->get_end() ) {
                            Int l = e_phase;
                            CDSPepPair *next_pcp;
                            i_cp++;
                            while (i_cp<pp->get_cp_num()) {
                                next_pcp = pp->get_cp(i_cp);
                                   if (next_pcp!=NULL) break;
                                i_cp++;
                            }
                            if (i_cp<pp->get_cp_num()) {
                                CDS * next_pc = next_pcp->get_cds();
                                codon = chr_obj->get_seq(pc->get_end()-l, l);
                                codon += chr_obj->get_seq(next_pc->get_begin(), 3-l);
                            }
                        } else {
                            codon = chr_obj->get_seq(begin, 3);
                        }
                        if (!codon.empty() && map_codon_to_phase.find(codon)==map_codon_to_phase.end() ) {
                            map_codon_to_phase[codon] = phase;
                            codons.push_back(codon);
                            phases.push_back(phase);
                            out_cdss.push_back(pc);
                            out_peps.push_back(pp);
                        } else {
                        }
                        break;
                    }
                }
            }
        }
    }

    return 0;
}

Int GTP::get_codon_by_pos(Gene *gene_obj, Int pos, vector<string> & codons, vector<Int> & phases, vector<CDS*> & out_cdss, vector<Pep*> & out_peps)
{
	Chrom * chr_obj = gene_obj->get_chrom();
	Int i;
    map<string, Int> map_codon_to_phase;
	for (i=0; i<gene_obj->get_tscript_num(); i++) {
		Tscript *pt = gene_obj->get_tscript(i);
		for (Int i_p=0; i_p<pt->get_pep_num(); i_p++) {
			Pep *pp = pt->get_pep(i_p);
			for (Int i_cp=0; i_cp<pp->get_cp_num(); i_cp++) {
				CDSPepPair *pcp = pp->get_cp(i_cp);
				CDS * pc = pcp->get_cds();
				Int b_phase = pcp->get_begin_phase();
				Int e_phase = pcp->get_end_phase();
				if ( pos>=pc->get_begin() && pos<pc->get_end()) {
					string codon;
					Int phase;
					phase = ( (pos - pc->get_begin()) + b_phase )%3;
					Int begin = pos - phase;
					if ( begin < pc->get_begin() ) {
						Int l = phase;
						CDSPepPair *prev_pcp;
						i_cp--;
						while (i_cp>=0) { 
							prev_pcp = pp->get_cp(i_cp);
							if (prev_pcp!=NULL) break;
							i_cp--;
						}
						if (i_cp>=0) {
							CDS * prev_pc = prev_pcp->get_cds();
							codon = chr_obj->get_seq(prev_pc->get_end()-l, l);
							codon += chr_obj->get_seq(pc->get_begin(), 3-l);
						}
					} else if ( begin+3 >= pc->get_end() ) {
						Int l = e_phase;
						CDSPepPair *next_pcp;
						i_cp++;
						while (i_cp<pp->get_cp_num()) {
							next_pcp = pp->get_cp(i_cp);
							   if (next_pcp!=NULL) break;
							i_cp++;
						}
						if (i_cp<pp->get_cp_num()) {
							CDS * next_pc = next_pcp->get_cds();
							codon = chr_obj->get_seq(pc->get_end()-l, l);
							codon += chr_obj->get_seq(next_pc->get_begin(), 3-l);
						}
					} else {
						codon = chr_obj->get_seq(begin, 3);
					}
					if (!codon.empty() && map_codon_to_phase.find(codon)==map_codon_to_phase.end() ) {
						map_codon_to_phase[codon] = phase;
						codons.push_back(codon);
						phases.push_back(phase);
						out_cdss.push_back(pc);
						out_peps.push_back(pp);
					} else {
					}
					break;
				}
			}
		}
	}
}

bool is_same_chrom(Chrom* g1, Chrom* g2)
{
    if (g1->get_id() == g2->get_id()) return true;
    else return false;
}

bool is_same_gene(Gene* g1, Gene* g2)
{
    if (g2->get_chrom()->get_id()==g1->get_chrom()->get_id() && g2->get_begin()>=g1->get_begin() && g2->get_end()<=g1->get_end()) return true;
    else return false;
}

bool is_same_cp_1(CDSPepPair* cp1, CDSPepPair* cp2)
{
    if (cp2->get_cds()->get_begin() >= cp1->get_cds()->get_begin() && cp2->get_cds()->get_end() <= cp1->get_cds()->get_end() ) return true;
    else return false;
}

/**
 * @brief class Variant
 */
Variant::Variant(Variant & variant)
{
    set(variant);
}

Variant::Variant(const Variant & variant )
{
    set(variant);
}

bool Variant::is_shift_frame(Int i)
{
    if (i<0 || i>=mut_seq_m.size()) {
        stringstream ss;
        ss << "In Variant::is_shift_frame, 1st parameter " << i << "should in [0," << mut_seq_m.size() << ")";
        Msg::error(ss);
        return false;
    } else {
        Int o_len = src_m.get_seq_len(0);
        Int m_len = mut_seq_m[i].get_seq_len();
        return ( (m_len-o_len)%3 == 0 );
    }
}

Int Variant::get_forw_read_num(Int i)
{
    if (i<0 || i>=forw_read_num_m.size()) {
        stringstream ss;
        ss << "In Variant::get_forw_read_num, 1st parameter " << i << "should in [0," << forw_read_num_m.size() << ")";
        Msg::error(ss);
        return 0;
    } else {
        return forw_read_num_m[i];
    }
}

Int Variant::get_rev_read_num(Int i)
{
    if (i<0 || i>=rev_read_num_m.size()) {
        stringstream ss;
        ss << "In Variant::get_rev_read_num, 1st parameter " << i << "should in [0," << rev_read_num_m.size() << ")";
        Msg::error(ss);
        return 0;
    } else {
        return rev_read_num_m[i];
    }
}

Int Variant::get_read_num(Int i)
{
    if (i<0 || i>=forw_read_num_m.size()) {
        stringstream ss;
        ss << "In Variant::get_read_num, 1st parameter " << i << "should in [0," << forw_read_num_m.size() << ")";
        Msg::error(ss);
        return 0;
    } else {
        return forw_read_num_m[i] + rev_read_num_m[i];
    }
}

Gene* Variant::get_gene(Int i) 
{
    if ( i<0 || i>=gene_m.size() ) {
        stringstream ss;
        ss << "In Variant::get_gene, 1st param : " << i << " should in [0," << gene_m.size() << ")";
        Msg::error(ss);
        return NULL;
    } else {
        return gene_m[i];
    }
}

Exon* Variant::get_exon(Int i) 
{
    if ( i<0 || i>=exon_m.size() ) {
        stringstream ss;
        ss << "In Variant::get_exon, 1st param : " << i << " should in [0," << exon_m.size() << ")";
        Msg::error(ss);
        return NULL;
    } else {
        return exon_m[i];
    }
}

CDS* Variant::get_cds(Int i) 
{
    if ( i<0 || i>=cds_m.size() ) {
        stringstream ss;
        ss << "In Variant::get_cds, 1st param : " << i << " should in [0," << cds_m.size() << ")";
        Msg::error(ss);
        return NULL;
    } else {
        return cds_m[i];
    }
}

BioSeq & Variant::get_mut_seq(Int i)
{ 
    if (i<0 || i>=mut_seq_m.size()) {
        stringstream ss;
        ss << "In get_mut_seq, 1st param : " << i << " should in [0," << mut_seq_m.size();
        Msg::error(ss);
    }
    return mut_seq_m[i];
}

VariantInfo * Variant::get_info(Int i)
{
    if (i<0 || i>=info_m.size()) {
        stringstream ss;
        ss << "In get_info, 1st param : " << i << " should in [0," << info_m.size();
        Msg::error(ss);
    } else {
        return info_m[i];
    }
}

Variant & Variant::set( Variant & variant )
{
    src_m      = variant.src_m;
    mut_seq_m  = variant.mut_seq_m;
    gene_m     = variant.gene_m;
    qual_m     = variant.qual_m;
    map_qual_m = variant.map_qual_m;
    genotype_m[0] = variant.genotype_m[0];
    genotype_m[1] = variant.genotype_m[1];
    forw_read_num_m = variant.forw_read_num_m;
    rev_read_num_m  = variant.rev_read_num_m;
    total_forw_read_num_m = variant.total_forw_read_num_m;
    total_rev_read_num_m  = variant.total_rev_read_num_m;
    variant_type_m = variant.variant_type_m;
    return *this;
}

Variant & Variant::set( const Variant & variant )
{
    src_m      = variant.src_m;
    mut_seq_m  = variant.mut_seq_m;
    gene_m     = variant.gene_m;
    qual_m     = variant.qual_m;
    map_qual_m = variant.map_qual_m;
    genotype_m[0] = variant.genotype_m[0];
    genotype_m[1] = variant.genotype_m[1];
    forw_read_num_m = variant.forw_read_num_m;
    rev_read_num_m  = variant.rev_read_num_m;
    total_forw_read_num_m = variant.total_forw_read_num_m;
    total_rev_read_num_m  = variant.total_rev_read_num_m;
    variant_type_m = variant.variant_type_m;
    return *this;
}

void Variant::set_forw_read_num(Int i, Int forw_read_num)
{
    if (i<0 || i>=forw_read_num_m.size()) {
        stringstream ss;
        ss << "In Variant::set_read_num, 1st parameter " << i << "should in [0," << forw_read_num_m.size() << ")";
        Msg::error(ss);
        return ;
    } else {
        total_forw_read_num_m -= forw_read_num_m[i];
        forw_read_num_m[i] = forw_read_num;
        total_forw_read_num_m += forw_read_num;
    }
}

void Variant::set_rev_read_num(Int i, Int rev_read_num)
{
    if (i<0 || i>=rev_read_num_m.size()) {
        stringstream ss;
        ss << "In Variant::set_read_num, 1st parameter " << i << "should in [0," << rev_read_num_m.size() << ")";
        Msg::error(ss);
        return ;
    } else {
        total_rev_read_num_m -= rev_read_num_m[i];
        rev_read_num_m[i]  = rev_read_num;
        total_rev_read_num_m += rev_read_num;
    }
}

void Variant::set_read_num(Int i, Int forw_read_num, Int rev_read_num)
{
    if (i<0 || i>=forw_read_num_m.size()) {
        stringstream ss;
        ss << "In Variant::set_read_num, 1st parameter " << i << "should in [0," << forw_read_num_m.size() << ")";
        Msg::error(ss);
        return ;
    } else {
        total_forw_read_num_m -= forw_read_num_m[i];
        forw_read_num_m[i] = forw_read_num;
        total_forw_read_num_m += forw_read_num;

        total_rev_read_num_m -= rev_read_num_m[i];
        rev_read_num_m[i]  = rev_read_num;
        total_rev_read_num_m += rev_read_num;
    }
}

void Variant::ins_read_num(Int forw_read_num, Int rev_read_num)
{
    forw_read_num_m.push_back(forw_read_num);
    total_forw_read_num_m += forw_read_num;
    rev_read_num_m.push_back(rev_read_num);
    total_rev_read_num_m += rev_read_num;
}

void Variant::clear_info()
{
    info_m.clear();
}

Variant & Variant::operator =(Variant & variant)
{
    return set(variant);
}

Variant & Variant::operator =(const Variant & variant)
{
    return set(variant);
}

bool less_variant_pos(Variant *v1, Variant *v2)
{
    SeqRegion *sr1 = v1->get_seq_region_pt();
    SeqRegion *sr2 = v2->get_seq_region_pt();
    string id1 = sr1->get_seq_pt()->get_id();
    string id2 = sr2->get_seq_pt()->get_id();
    if ( id1 < id2 ) {
        return true;
    } else if ( id1 == id2 ) {
        if ( sr1->get_begin() < sr2->get_begin() ) {
            return true;
        } else if ( sr1->get_begin() == sr2->get_begin() ) {
            if ( sr1->get_end() < sr2->get_end() ) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool less_variant_range(Variant *v1, Variant *v2)
{
    SeqRegion *sr1 = v1->get_seq_region_pt();
    SeqRegion *sr2 = v2->get_seq_region_pt();
    if ( sr1->get_begin() < sr2->get_begin() ) {
        return true;
    } else if ( sr1->get_begin() == sr2->get_begin() ) {
        if ( sr1->get_end() < sr2->get_end() ) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

/**
 * @brief class VCFVariantInfo
 */
// {{{

Int VariantInfo::get_forw_read_num(Int i)
{
    if (i<0 || i>=forw_read_num_m.size()) {
        stringstream ss;
        ss << "In VariantInfo::get_forw_read_num, 1st parameter " << i << "should in [0," << forw_read_num_m.size() << ")";
        Msg::error(ss);
        return 0;
    } else {
        return forw_read_num_m[i];
    }
}

Int VariantInfo::get_rev_read_num(Int i)
{
    if (i<0 || i>=rev_read_num_m.size()) {
        stringstream ss;
        ss << "In VariantInfo::get_rev_read_num, 1st parameter " << i << "should in [0," << rev_read_num_m.size() << ")";
        Msg::error(ss);
        return 0;
    } else {
        return rev_read_num_m[i];
    }
}

Int VariantInfo::get_read_num(Int i)
{
    if (i<0 || i>=forw_read_num_m.size()) {
        stringstream ss;
        ss << "In VariantInfo::get_read_num, 1st parameter " << i << "should in [0," << forw_read_num_m.size() << ")";
        Msg::error(ss);
        return 0;
    } else {
        return forw_read_num_m[i] + rev_read_num_m[i];
    }
}

string VariantInfo::get_genotype_str()
{
    stringstream out;
    out << genotype_m[0] << "/" << genotype_m[1];
    return out.str();
}

void VariantInfo::set_forw_read_num(Int i, Int forw_read_num)
{
    if (i<0 || i>=forw_read_num_m.size()) {
        stringstream ss;
        ss << "In VariantInfo::set_read_num, 1st parameter " << i << "should in [0," << forw_read_num_m.size() << ")";
        Msg::error(ss);
        return ;
    } else {
        total_forw_read_num_m -= forw_read_num_m[i];
        forw_read_num_m[i] = forw_read_num;
        total_forw_read_num_m += forw_read_num;
    }
}

void VariantInfo::set_rev_read_num(Int i, Int rev_read_num)
{
    if (i<0 || i>=rev_read_num_m.size()) {
        stringstream ss;
        ss << "In VariantInfo::set_read_num, 1st parameter " << i << "should in [0," << rev_read_num_m.size() << ")";
        Msg::error(ss);
        return ;
    } else {
        total_rev_read_num_m -= rev_read_num_m[i];
        rev_read_num_m[i]  = rev_read_num;
        total_rev_read_num_m += rev_read_num;
    }
}

void VariantInfo::set_read_num(Int i, Int forw_read_num, Int rev_read_num)
{
    if (i<0 || i>=forw_read_num_m.size()) {
        stringstream ss;
        ss << "In VariantInfo::set_read_num, 1st parameter " << i << "should in [0," << forw_read_num_m.size() << ")";
        Msg::error(ss);
        return ;
    } else {
        total_forw_read_num_m -= forw_read_num_m[i];
        forw_read_num_m[i] = forw_read_num;
        total_forw_read_num_m += forw_read_num;

        total_rev_read_num_m -= rev_read_num_m[i];
        rev_read_num_m[i]  = rev_read_num;
        total_rev_read_num_m += rev_read_num;
    }
}

void VariantInfo::ins_read_num(Int forw_read_num, Int rev_read_num)
{
    forw_read_num_m.push_back(forw_read_num);
    total_forw_read_num_m += forw_read_num;
    rev_read_num_m.push_back(rev_read_num);
    total_rev_read_num_m += rev_read_num;
}

// }}}

/**
 * @brief Express
 */
// {{{
Express::Express(BioSeq *bio_seq, short strand): SeqRegion(bio_seq, strand)
{
    for (Int i=0; i<3; i++) fpkm_m[i]=0;
}

Express::Express(Express & express)
{
}

Express::Express(const Express & express)
{
}

Express::~Express()
{
}

double Express::get_fpkm()
{
    return fpkm_m[0];
}

double Express::get_fpkm_low()
{
    return fpkm_m[1];
}

double Express::get_fpkm_high()
{
    return fpkm_m[2];
}

void Express::set(Express & express)
{
    SeqRegion::set(express);
    for (Int i=0; i<3; i++) fpkm_m[i] = express.fpkm_m[i];
}

void Express::set_fpkm(double fpkm)
{
    fpkm_m[0] = fpkm;
}

void Express::set_fpkm_conf(double low, double high)
{
    fpkm_m[1] = low;
    fpkm_m[2] = high;
}

void Express::set_fpkm(double fpkm, double low, double high)
{
    fpkm_m[0] = fpkm;
    fpkm_m[1] = low;
    fpkm_m[2] = high;
}
// }}}

};
