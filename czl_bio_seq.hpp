#ifndef CZL_BIO_SEQ_HPP
#define CZL_BIO_SEQ_HPP

#include "seq/bit2_seq.hpp"

namespace czl_bio {

int detect_base_qual_type(vector<string> const & qual_v);

int trim_seq_by_qual(string const & qual,
           int win_size, int step_size, 
           int qual_thres, float win_frac_thres, int win_min,
           int cut_pos[2]);

int trim_seq_by_qual(string const & qual, int begin, int length,
           int win_size, int step_size, 
           int qual_thres, float win_frac_thres, int win_min,
           int cut_pos[2]);

int trim_seq_by_qual_m2(string const & qual, int begin, int length,
            int qual_thres, float frac_thres, int size_thres,
            int cut_pos[2]);

int trim_seq_by_qual_m3(string const & qual,
           int win_size, int step_size, 
           int qual_thres, float win_frac_thres,
           int cut_pos[2]);

int trim_seq_by_qual_m3(string const & qual, int begin, int length,
           int win_size, int step_size, 
           int qual_thres, float win_frac_thres,
           int cut_pos[2]);

/**
 * @brief trim sequence by quality, method 4: devide the sequence into 
 *        alternative good runs (base quality >= 'qual_thres') and bad runs 
 *        (base quality < 'qual_thres'). Foreach pair of good runs, such as 
 *        run i from b_i to e_i, run j from b_j to e_j, compute the fraction
 *        of good base pairs. if the fraction >= 'frac_thres', store it. At last
 *        get the best run pairs with fraction >= 'frac_thres' and return the 
 *        cutting position of the run pairs. Best can be the longest sequence or  *        largest number of good base pairs.
 * @param  qual  quality of sequence
 * @param  qual_thres  good quality thresthold (greater than it is good)
 * @param  frac_thres  return region must have at least 'frac_thres' of good 
 *                     base pairs 
 * @param  cut_pos  region range returned
 * @return  0 if success, not 0 if other.
 * @remark  This method is exactly when you want to look for the longest 
 *          trimmed sequence with fraction of good bps > 'frac_thres', but slow,
 *          O(n^2), n is the size of good runs (<= half of the sequence length).
 * @sa  int trim_seq_by_qual_m4(string const &, size_t, size_t, int, float,
 *      int [2])
 */
int trim_seq_by_qual_m4(string const & qual, 
           int qual_thres, float frac_thres,
           int cut_pos[4]);
/**
 * @brief trim sequence by quality, method 4: devide the sequence into 
 *        alternative good runs (base quality >= 'qual_thres') and bad runs 
 *        (base quality < 'qual_thres'). Foreach pair of good runs, such as 
 *        run i from b_i to e_i, run j from b_j to e_j, compute the fraction
 *        of good base pairs. if the fraction >= 'frac_thres', store it. At last
 *        get the best run pairs with fraction >= 'frac_thres' and return the 
 *        cutting position of the run pairs. Best can be the longest sequence or  *        largest number of good base pairs.
 * @param  qual_thres  good quality thresthold (greater than it is good)
 *
 * @param  qual  quality of sequence
 * @param  begin  begin position of the sequence for trimming 
 * @param  length length of the sequence (from 'begin') for trimming 
 * @param  qual_thres  good quality thresthold (greater than it is good)
 * @param  frac_thres  return region must have at least 'frac_thres' of good 
 *                     base pairs 
 * @param  cut_pos  region range returned, 0-1: longest sequence position, 
 *                  2-3: trimmed position of sequence with maximal nubmer of 
 *                  good base pairs.
 * @return  0 if success, not 0 if other.
 * @remark  This method is exactly when you want to look for the longest 
 *          trimmed sequence with fraction of good bps > 'frac_thres', but slow,
 *          O(n^2), n is the size of good runs (<= half of the sequence length).
 * @sa  int trim_seq_by_qual_m4(string const &, size_t, size_t, int, float,
 *      int [2])
 */
int trim_seq_by_qual_m4(string const & qual, int begin, int length,
           int qual_thres, float frac_thres,
           int cut_pos[4]);

int trim_seq_by_qual_m4(string const & qual1, string const & qual2,
           int qual_thres,
           float frac_thres, float frac_thres1, float frac_thres2,
           int cut_pos[8]);
int trim_seq_by_qual_m4(string const & qual1, int begin1, int length1,
           string const & qual2, int begin2, int length2,
           int qual_thres,
           float frac_thres, float frac_thres1, float frac_thres2,
           int cut_pos[8]);

short is_seq_low_qual(string const & qual,
           int qual_thres, float frac_thres, int min_len);

int get_high_qual_len(string const & qual, int qual_thres);

int get_high_qual_len(string const & qual, int begin, int len, int qual_thres);

};

#endif
