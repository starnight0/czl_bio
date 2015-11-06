#ifndef CZL_VCF_H
#define CZL_VCF_H
/* @file czl_vcf.h
 * @brief manipulate VCF record, read/write VCF files
 */
#include "common.hh"

class Variant {
}

/* @brief class VarInfo
class VarInfo {
public:
    int32_t total_allele_t;  // e.g. from VCF INFO 'AN'
    vector<int32_t> allele_count_v_t; // e.g. from VCF 'AC' for REF and each ALT
    int32_t total_depth_t;  // total read depth
    vector<int32_t> allele_depth_v_t; // e.g. from VCF 'DP'(read depth before filter)
    float qual_t;   // quality, e.g from VCF 'QUAL'
    float rsm_map_qual; // mapping quality, e.g. from VCF MQ: RMS mapping quality
//  vector<int32_t> mle_allele_count_v;  // MLEAC for REF and each ALT
//  vector<int32_t> mle_allele_freq_v;   // MLEAE for REF and each ALT
//  float strand_bias;     // FS
//  int map_qual0_total;   // MQ0
//  float qual_by_depth; // QD: quality depth
//  float base_qual_rank_sum; // BaseQRankSum
//  float read_pos_rank_sum;  // ReadPosRankSum
};
 */

/* class VarSample
 */
class VarSample {
public: 
    int32_t total_depth_t;
    vector<int32_t> allele_depth_v_t;  // AD for REF and each ALT
//  int32_t a_read_depth;     // DP, approximate read depth, after filter
    float geno_qual;           // genotype quality, e.g. from VCF 'GQ'
    vector<int8_t> genotype_v; // e.g. from VCF 'GT' (genotype)
    vector<float> phred_likelihood_v;  // PL
    int32_t flag;
    void* data_t;

    VarSample();
    VarSample(const VarSample & vs);
    VarSample & operator=(const VarSample & vs);
};

/* class Var
 */
class Variant {
public:
    int32_t id_t;
    int32_t seq_id_t;
    int64_t pos_t;
    vector<string> allele_v_t;  // allele sequence, first is ref
    int32_t total_allele_t;  // e.g. from VCF INFO 'AN'
    vector<int32_t> allele_count_v_t; // e.g. from VCF 'AC' for REF and each ALT
    int32_t total_depth_t;  // total read depth
    vector<int32_t> allele_depth_v_t; // e.g. from VCF 'DP'(read depth before filter)
    float qual_t;   // quality, e.g from VCF 'QUAL'
    float rsm_map_qual_t; // mapping quality, e.g. from VCF MQ: RMS mapping quality

    vector<VarSample*> sample_v_t;

    int32_t flag_t;
    void* data_t;

    Variant();
    Variant(const Variant & v);
    Variant & operator =(const Variant & v);
};

#endif
