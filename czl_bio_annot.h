#ifndef ANNOT_HH
#define ANNOT_HH

// #include <boost/interval.hpp>
#include "czl_common.h"
#include "czl_bio_base.h"

/**
 * @brief  classes for gene-transcript-peptide struture
 * @author Zelin Chen
 */

namespace czl_bio {

enum BioStatus {
    CODING=0x0001,  ///< coding genes
    NOVEL=0x0002,  ///< novel, function not known
    PSEUDO=0x0004,  ///< pseudo gene
    PARTIAL=0x0008,  ///< partial gene
    CONSTITUTIVE=0x0010,  ///< constitutive element (exon)
    CANONICAL=0x0020  ///< canonical element (transcript, peptide)
};

enum BioStatusBit {
    CODING_BIT = 0,  ///< coding genes
    NOVEL_BIT  = 1,  ///< novel, function not known
    PSEUDO_BIT = 2,  ///< pseudo gene
    PARTIAL_BIT= 3,  ///< partial gene
    CONSTITUTIVE_BIT= 4,  ///< constitutive element (exon)
    CANONICAL_BIT   = 5  ///< canonical element (transcript, peptide)
};

enum GenomeRegion {
    GR_CDS = 1,
    GR_UTR = 2,
    GR_INTRON     = 3,
    GR_INTER_GENE = 4,
    GR_CONS       = 5,
};

enum VariantType {
    VT_NULL = 0,
    VT_SNP  = 1,
    VT_INDEL = 2,
    VT_SNP_INDEL = 3,
};

class Chrom;
class Gene;
class Tscript;
class Pep;
class Exon;
class CDS;
class CDSPepPair;
class ExonTscriptPair;
class Variant;
class VariantInfo;
class Express;

/**
 *
 */
class Chrom: public BioSeq {
protected:
    Int coord_system_id_m;
public:
    Chrom(): coord_system_id_m(-1) {}
    virtual ~Chrom() {};

    Int get_coord_system_id() {return coord_system_id_m;}

    Int set_coord_system_id(Int coord_system_id) { coord_system_id_m = coord_system_id; }
};
bool is_same_chrom(Chrom* g1, Chrom* g2);

/**
 * Class Gene : inherit from SeqObj
 */
class Gene: virtual public BioObj, virtual public SeqRegion {
protected:
    vector<Tscript*> tscript_m;
    Tscript * canon_tscript_m;
public:
    Gene();
    virtual ~Gene();

    Chrom* get_chrom();
    Tscript* get_tscript(Int i);
    Tscript* get_canon_tscript();
    Int get_tscript_num(short is_include_null=1);
	void * get_data() { return BioObj::get_data();}

    void set_chrom(Chrom* chrom);
    void set_tscript(Int i, Tscript* tscript);
    void set_canon_tscript(Tscript* t);
	void set_data(void * data) { BioObj::set_data(data);}

    void ins_tscript(Tscript* t);

    Tscript* rm_tscript(Int i);
//    void find_tscript(Tscript* t);

    void compact_vec();
    void sort_tscript( bool (*cmp_fh)(Tscript*, Tscript*) );

    friend bool less_gene_pos(Gene *t1, Gene *t2);
    friend bool less_gene_pos(Gene &t1, Gene &t2);
};
bool less_gene_pos(Gene *t1, Gene *t2);
bool less_gene_pos(Gene &t1, Gene &t2);
// bool less_gene_range(Gene *t1, Gene *t2);
// bool less_gene_range(Gene &t1, Gene &t2);
bool less_gene_chromsymbol_pos(Gene *t1, Gene *t2);
bool is_same_gene(Gene* g1, Gene* g2);

/**
 * Class Tscript : inherit from SeqObj
 */
class Tscript: public BioObj {
protected:
    Chrom *chrom_m; ///< point to Chrom
    Gene* gene_m;
    vector<Pep*> pep_m;
    Int n__null_pep_m;
    vector<ExonTscriptPair*> et_m;
    Int n__null_et_m;
    bool is_et_sorted_m;
    Pep* canon_pep_m;
    short strand_m;
public:
    Tscript();
    virtual ~Tscript();
    Chrom* get_chrom();
    Gene* get_gene();
    Pep* get_pep(Int i);
    Pep* get_canon_pep();
    Exon* get_exon(Int i);
    ExonTscriptPair* get_et(Int i);
    Int get_pep_num(short is_include_null=1);
    Int get_et_num(short is_include_null=1);
    Int get_exon_num(short is_include_null=1);
    Int get_begin();
    Int get_begin_in_gene() { return get_begin() - gene_m->get_begin(); }
    Int get_end();
    Int get_end_in_gene() { return get_end() - gene_m->get_begin(); }
    Int get_span_len();
    Int get_seq_len();
    Int get_strand() { return strand_m; }

    void set_chrom(Chrom* chrom);
    void set_gene(Gene* gene);
    void set_pep(Int i, Pep* pep);
    void set_canon_pep(Pep* pep);
    void set_et(Int i, ExonTscriptPair* et);
    void set_strand(short strand) { strand_m = strand; }

    void ins_gene(Gene* gene);
    void ins_pep(Pep* pep);
    void ins_et(ExonTscriptPair* et);

    Gene* rm_gene(Int i);
    Pep* rm_pep(Int i);
    ExonTscriptPair* rm_et(Int i);

    void compact_vec();
    void sort_et();
    void sort_et( bool(*cmp_fh)(ExonTscriptPair*, ExonTscriptPair*) );

};
bool less_tscript_pos(Tscript *t1, Tscript *t2);
bool less_tscript_pos(Tscript &t1, Tscript &t2);

/**
 * Class Pep : inherit from SeqObj
 */
class Pep: public BioObj {
protected:
    Chrom *chrom_m;
    Tscript* tscript_m;
    vector<CDSPepPair*> cp_m;
    Int n__null_cp_m;
    bool is_cp_sorted_m;
    Int nt_len_m;
public:
    Pep(): n__null_cp_m(0), is_cp_sorted_m(false) {}
    virtual ~Pep() {}

    Chrom* get_chrom();
    short get_strand();
    Tscript* get_tscript() { return tscript_m; }
    CDS* get_cds(Int i);
    CDSPepPair* get_cp(Int i);
    Int get_cp_num(short is_include_null=1);
    Int get_cds_num(short is_include_null=1);
    Int get_begin();
    Int get_end();
    Int get_nt_len();
    Int get_aa_len();
    Int get_seq_len() { return get_aa_len(); }
    Int get_begin_in_tscript() { return get_begin() - tscript_m->get_begin(); }
    Int get_end_in_tscript() { return get_end() - tscript_m->get_begin(); }

    void set_chrom(Chrom* chrom);
    void set_tscript(Tscript* tscript);
    void set_cp(Int i, CDSPepPair* cp);
    void set_nt_len(Int len);
//  Int set_cds_len(Int l) { SeqRegion::set_seq_len(l); }

    void ins_tscript(Tscript* tscript);
    void ins_cp(CDSPepPair* cp);

    void cal_len_from_cp();

    CDSPepPair* rm_cp(Int i);

    void compact_vec();
//    void sort_cds( bool(*cmp_fh)(CDS*, CDS*) );
    void sort_cp();
    void sort_cp( bool(*cmp_fh)(CDSPepPair*, CDSPepPair*) );
};

/**
 * Class Exon: inherit from SeqObj
 */
class Exon: public BioObj, public SeqRegion {
protected:
    vector<ExonTscriptPair*> et_m;
    vector<CDS*> cds_m;
public:
    Tscript* get_tscript(Int i);
    ExonTscriptPair* get_et(Int i);
    Int get_tscript_num(short is_include_null=1);
    Int get_et_num(short is_include_null=1);
    CDS* get_cds(Int i);
    Int get_cds_num(short is_include_null=1);
	Chrom * get_chrom();

    void set_et(Int i, ExonTscriptPair* et);
    void set_cds(Int i, CDS* cds);

    void ins_et(ExonTscriptPair* et);
    void ins_cds(CDS* cds);

    ExonTscriptPair* rm_et(Int i);
    CDS* rm_cds(Int i);

    void compact_vec();
};
bool less_exon_pos(Exon *t1, Exon *t2);

/**
 * Class CDS: inherit from SeqObj
 */
class CDS: public BioObj, public SeqRegion {
protected:
    vector<CDSPepPair*> cp_m;
    Exon* exon_m;
public:
    Pep* get_pep(Int i);
    CDSPepPair* get_cp(Int i);
    Int get_pep_num(short is_include_null=1);
    Int get_cp_num(short is_include_null=1);
    Exon* get_exon() { return exon_m; }
//  Exon* get_exon(Int i);
//  Int get_exon_num(short is_include_null=1);

//  void set_pep(Int i, Pep* pep);
    void set_cp(Int i, CDSPepPair* cp);
//  void set_exon(Int i, Exon* exon);
    void set_exon(Exon* exon) { exon_m = exon; }
	Chrom * get_chrom();

//  void ins_pep(Pep* pep);
    void ins_cp(CDSPepPair* cp);
//  void ins_exon(Exon* exon);

//  Pep* rm_pep(Int i);
    CDSPepPair* rm_cp(Int i);
//  Exon* rm_exon(Int i);

    void compact_vec();
};
bool less_cds_pos(CDS *t1, CDS *t2);

class ExonTscriptPair: public BioObjPair {
protected:
    Int begin_in_tscript_m;
    Int exon_rank_m;
    Int exon_begin_phase_m;
    Int exon_end_phase_m;
public:
    ExonTscriptPair(Exon* exon=NULL, Tscript * tscript=NULL, Int begin_in_tscript = -1, Int rank=-1, Int begin_phase=3, Int end_phase=3): BioObjPair(exon, tscript), begin_in_tscript_m(begin_in_tscript),  exon_rank_m(rank), exon_begin_phase_m(begin_phase), exon_end_phase_m(end_phase) {}

    Exon* get_exon() { return static_cast<Exon*>(get_obj(0)); }
    Tscript* get_tscript() { return static_cast<Tscript*>(get_obj(1)); }
    Int get_begin_in_tscript() { return begin_in_tscript_m; }
    Int get_rank() {return exon_rank_m;}
    Int get_begin_phase() {return exon_begin_phase_m;}
    Int get_end_phase() {return exon_end_phase_m;}

    void set_exon(Exon* exon) { set_obj(0, exon);}
    void set_tscript(Tscript* tscript) { set_obj(1, tscript);}
    void set_begin_in_tscript(Int begin_in_tscript) { begin_in_tscript_m = begin_in_tscript; }
    void set_rank(Int rank) {exon_rank_m = rank;}
    void set_phase(Int begin_phase, Int end_phase) {
        exon_begin_phase_m = begin_phase; 
        exon_end_phase_m = end_phase; 
    }

    friend bool less_et_exon_pos(ExonTscriptPair *et1, ExonTscriptPair *et2);
    friend bool less_et_tscript_pos(ExonTscriptPair *et1, ExonTscriptPair *et2);
};
bool less_et_exon_pos(ExonTscriptPair *et1, ExonTscriptPair *et2);
bool less_et_tscript_pos(ExonTscriptPair *et1, ExonTscriptPair *et2);

class CDSPepPair: BioObjPair {
protected:
    Int begin_in_tscript_m;
    Int cds_rank_m;
    Int cds_begin_phase_m;
    Int cds_end_phase_m;
public:
    CDSPepPair(CDS* cds=NULL, Pep * pep=NULL, Int begin_in_tscript=-1, Int rank=-1, Int begin_phase=3, Int end_phase=3): BioObjPair(cds, pep), begin_in_tscript_m(begin_in_tscript), cds_rank_m(rank), cds_begin_phase_m(begin_phase), cds_end_phase_m(end_phase) {}

    void set_cds(CDS* cds) { set_obj(0, cds);}
    void set_pep(Tscript* pep) { set_obj(1, pep);}
    void set_begin_in_tscript(Int begin_in_tscript) {begin_in_tscript_m = begin_in_tscript;}
    void set_rank(Int rank) {cds_rank_m = rank;}
    void set_phase(Int begin_phase, Int end_phase) {
        cds_begin_phase_m = begin_phase; 
        cds_end_phase_m = end_phase; 
    }

    CDS* get_cds() { return static_cast<CDS*>(get_obj(0)); }
    Pep* get_pep() { return static_cast<Pep*>(get_obj(1)); }
    Int get_begin_in_tscript() { return begin_in_tscript_m;}
    Int get_rank() {return cds_rank_m;}
    Int get_begin_phase() {return cds_begin_phase_m;}
    Int get_end_phase() {return cds_end_phase_m;}

    friend bool less_cp_cds_pos(CDSPepPair *cp1, CDSPepPair *cp2);
    friend bool is_same_cp_1(CDSPepPair* cp1, CDSPepPair* cp2);
};
bool less_cp_cds_pos(CDSPepPair *cp1, CDSPepPair *cp2);
bool is_same_cp_1(CDSPepPair* cp1, CDSPepPair* cp2);

class GTP: public BioObj {
protected:
    vector<Chrom*> chrom_m; 
    vector<Gene*> gene_m; 
    vector<Tscript*> tscript_m; 
    vector<Pep*> pep_m; 
    vector<Exon*> exon_m; 
    vector<CDS*> cds_m; 
    vector<ExonTscriptPair*> et_m;
    vector<CDSPepPair*> cp_m;
    vector<Variant*> variant_m;

    map<string, Variant*> map_variant_to_obj_m;
    map<string, Chrom*> map_chrom_id_to_obj_m;
    map<string, Chrom*> map_chrom_symbol_to_obj_m;
    map<string, Gene*> map_gene_id_to_obj_m;
    map<string, Tscript*> map_tscript_id_to_obj_m;
    map<string, Pep*> map_pep_id_to_obj_m;
public:
    GTP();
    virtual ~GTP();

    vector<Gene*> & get_gene();
    Gene* get_gene(Int i);
    vector<Tscript*> & get_tscript();
    Tscript* get_tscript(Int i);
    Pep* get_pep(Int i);
    vector<Pep*> & get_pep();
    vector<Exon*> & get_exon();
    Exon* get_exon(Int i);
    vector<CDS*> & get_cds();
    CDS* get_cds(Int i);
//  vector<Variant*> & get_variant();

    Int get_chrom_num();
    Int get_gene_num();
    Int get_tscript_num();
    Int get_pep_num();
    Int get_exon_num();
    Int get_cds_num();

    Chrom* get_chrom_by_id(string & id);
    Chrom* get_chrom_by_symbol(string & symbol);
    Gene* get_gene_by_id(string & id);
    Tscript* get_tscript_by_id(string & id);
    Pep* get_pep_by_id(string & id);
    Variant* get_variant_by_id_range_mut(string& id, Int begin, Int end, string & mut_seq);

    void set_chrom_by_id(string & id, Chrom* chrom);
    void set_gene_by_id(string & id, Gene* g);
    void set_tscript_by_id(string & id, Tscript* t);
    void set_pep_by_id(string & id, Pep *p);
    void set_variant_by_id_range_mut(string& id, Int begin, Int end, string & mut_seq, Variant *v);

    void read_file_seq_from_fasta(string file);
    void read_from_db(string dir);
    void sort_gene( bool (*cmp_fh)(Gene*, Gene*) );
    void sort_tscript( bool (*cmp_fh)(Tscript*, Tscript*) );
	void sort_exon( bool(*cmp_fh)(Exon*, Exon*) );
    void sort_cds( bool (*cmp_fh)(CDS*, CDS*) );

    Chrom* search_chrom(bool (*is_fh)(Chrom*, Chrom*), Chrom* c);
    Gene* search_gene(bool (*is_fh)(Gene*, Gene*), Gene* g);
	void search_gene(bool (*is_fh)(Gene*, Gene*), Gene* g, vector<Gene*> & out_genes);
    Tscript* search_tscript(bool (*is_fh)(Tscript*, Tscript*), Tscript* g);
	Int get_codon_by_pos(string & chr_oid, Int pos, vector<string> & codons, vector<Int> & phases, vector<CDS*> & out_cdss, vector<Pep*> & out_peps, vector<Gene*> &out_genes);
	Int get_codon_by_pos(Chrom *chr_obj, Int pos, vector<string> & codons, vector<Int> & phases, vector<CDS*> & out_cdss, vector<Pep*> & out_peps, vector<Gene*> &out_genes);
	Int get_codon_by_pos(Gene *gene_obj, Int pos, vector<string> & codons, vector<Int> & phases, vector<CDS*> & out_cdss, vector<Pep*> & out_peps);
};

class Variant: public BioObj {
protected:
    SeqRegion src_m;
    vector<BioSeq> mut_seq_m; ///< mutated sequence, first one is ref sequence
    vector<Gene*> gene_m;
    vector<Exon*> exon_m;
    vector<CDS*> cds_m;
    float qual_m;
    float map_qual_m;  ///< VCF:MQ
//  float phred_qual__is_same_m; ///< VCF:FQ
    Int genotype_m[2];
    vector<Int> forw_read_num_m;  ///<  forward mapping read number for REF and each Mut
    vector<Int> rev_read_num_m;  ///<  reverse mapping read number for REF and each Mut
    Int total_forw_read_num_m;
    Int total_rev_read_num_m;
    vector<VariantInfo*> info_m;
    VariantType variant_type_m;

public:
    Variant(Variant & variant);
    Variant(const Variant & variant);
    Variant(): total_forw_read_num_m(0), total_rev_read_num_m(0), variant_type_m(VT_NULL) {}
    virtual ~Variant() {}

    float get_qual() {return qual_m; }
    Int get_total_read_num() {return total_forw_read_num_m + total_rev_read_num_m; }
    Int get_total_forw_read_num() {return total_forw_read_num_m;}
    Int get_total_rev_read_num() {return total_rev_read_num_m;}
    Int get_forw_read_num(Int i); 
    Int get_rev_read_num(Int i); 
    Int get_read_num(Int i); 
    SeqRegion & get_seq_region_ref() { return src_m; }
    SeqRegion * get_seq_region_pt() { return &src_m; }
    Int get_begin() { return src_m.get_begin(); }
    Int get_end() { return src_m.get_end(); }
    vector<Gene*> & get_gene() { return gene_m; }
    Gene* get_gene(Int i);
    Int get_gene_num() {return gene_m.size();}
    vector<Exon*> & get_exon() { return exon_m; }
    Exon* get_exon(Int i);
    Int get_exon_num() {return exon_m.size();}
    vector<CDS*> & get_cds() { return cds_m; }
    CDS* get_cds(Int i);
    Int get_cds_num() {return cds_m.size();}
    vector<BioSeq> & get_mut_seq() { return mut_seq_m; }
    BioSeq & get_mut_seq(Int i);
    string & get_mut_seq_str(Int i, string & out) { out = get_mut_seq(i).get_seq(); return out;}
    float get_map_qual() { return map_qual_m; }
    Int get_genotype(Int i) { 
        if ( i<0 || i>=2 ) {
            stringstream ss;
            ss << "In get_mut_seq, 1st param : " << i << " should in [0,2)";
            Msg::error(ss);
            return -1;
        } else {
            return genotype_m[i]; 
        }
    }
    Int get_site_num() { return mut_seq_m.size(); }
    VariantType get_type() { return variant_type_m; }
    VariantInfo * get_info(Int i);
    Int get_info_num() { return info_m.size(); }


    Variant & set( Variant & variant );
    Variant & set( const Variant & variant );
    void set_qual(float qual) { qual_m = qual; }
    void set_forw_read_num(Int i, Int forw_read_num);
    void set_rev_read_num(Int i, Int forw_read_num);
    void set_read_num(Int i, Int forw_read_num, Int rev_read_num); 
    void set_map_qual(float map_qual) { map_qual_m = map_qual; }
    void set_genotype(Int t1, Int t2) { genotype_m[0]=t1; genotype_m[1]=t2; }
    void set_type(VariantType type) { variant_type_m = type; }

    void ins_gene(Gene* gene) { gene_m.push_back(gene); }
    void ins_exon(Exon* exon) { exon_m.push_back(exon); }
    void ins_cds(CDS* cds) { cds_m.push_back(cds); }
    void ins_mut_seq(string & seq) { 
        BioSeq bs(seq); 
        mut_seq_m.push_back(bs);
    }

    void ins_read_num(Int forw_read_num, Int rev_read_num); 
    void ins_mut_seq(BioSeq & seq) { mut_seq_m.push_back(seq); }
    bool is_shift_frame(Int i);
    void ins_info(VariantInfo *info) { info_m.push_back(info); }
//  void set_orgin_seq(bpp::Sequence *seq) { SeqObj::set_seq(seq); }
//  void set_mut_seq(bpp::Sequence *seq) { mut_seq = seq; }

	void clear_info();

    Variant & operator =(Variant & variant);
    Variant & operator =(const Variant & variant);
};
bool less_variant_pos(Variant *v1, Variant *v2);
bool less_variant_range(Variant *v1, Variant *v2);

/**
 * class IVariantInfo for individual variant information
 */
class VariantInfo: public BioObj {
protected:
    Variant * variant_m;
    Int genotype_m[2];
    float map_qual_m;
//  float qual_m;
//  float phred_qual__is_same_m; ///< VCF:FQ
    vector<Int> forw_read_num_m;  ///<  forward mapping read number for REF and each Mut
    vector<Int> rev_read_num_m;  ///<  reverse mapping read number for REF and each Mut
    Int total_forw_read_num_m;
    Int total_rev_read_num_m;
public:
    VariantInfo(Variant *variant=NULL): total_forw_read_num_m(0), total_rev_read_num_m(0) {variant_m = variant;}
    virtual ~VariantInfo() {};

    Variant* get_variant() { return variant_m; }
    float get_map_qual() { return map_qual_m; }
    Int get_total_read_num() {return total_forw_read_num_m + total_rev_read_num_m; }
    Int get_total_forw_read_num() {return total_forw_read_num_m;}
    Int get_total_rev_read_num() {return total_rev_read_num_m;}
    Int get_forw_read_num(Int i); 
    Int get_rev_read_num(Int i); 
    Int get_read_num(Int i); 
    Int get_genotype(Int i) { 
        if ( i<0 || i>=2 ) {
            stringstream ss;
            ss << "In get_mut_seq, 1st param : " << i << " should in [0,2)";
            Msg::error(ss);
            return -1;
        } else {
            return genotype_m[i]; 
        }
    }
	string get_genotype_str();

    void set_variant(Variant* variant) { variant_m = variant; }
    void set_forw_read_num(Int i, Int forw_read_num);
    void set_rev_read_num(Int i, Int forw_read_num);
    void set_read_num(Int i, Int forw_read_num, Int rev_read_num); 
    void set_map_qual(float map_qual) { map_qual_m = map_qual; }
    void set_genotype(Int t1, Int t2) { genotype_m[0]=t1; genotype_m[1]=t2; }

    void ins_read_num(Int forw_read_num, Int rev_read_num); 

};

class Express: public SeqRegion {
protected:
	double fpkm_m[3]; ///< 0:FPKM and 1,2:confident interval
public:
	Express(BioSeq *seq_obj=NULL, short strand=0);
	Express(Express & express);
	Express(const Express & express);
	virtual ~Express();

	double get_fpkm();
	double get_fpkm_low();
	double get_fpkm_high();

	void set(Express & express);
	void set_fpkm(double fpkm);
	void set_fpkm_conf(double low, double high);
	void set_fpkm(double fpkm, double low, double high);
};


}; // end of namespace czl_bio

#endif
