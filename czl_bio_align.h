/* czl_bio_align.h
 * @brief class for represent alignment object
 * @author ZELIN CHEN
 * @date 2015-08-15
 */
#ifndef CZL_BIO_ALIGN_H
#define CZL_BIO_ALIGN_H

#include "czl_common.h"
#include "czl_bio_base.h"

using namespace czl_bio;
namespace czl_bio {
    template <typename Pos_T, typename ID_T> class AlignBlock;
    class AlignStat;

    /* Alignment
     *
     * It's represented by a tree with treenode as alignment blocks
     * 
     * REMARK!
     * - Changing number of sequence is expensive
     */
    // {{{
    template <typename Pos_T, typename ID_T>
    class Align: public AlignBlock<Pos_T, ID_T> 
    {
    public:
        Align(): AlignBlock<Pos_T, ID_T>(), consensus_i_m(0)
        {
        }
        Align(Align & src) : AlignBlock<Pos_T, ID_T>(src)
        {
            seq_id_vm = src.seq_id_vm;
            consensus_i_m = src.consensus_i_m;
        }

        Align & operator =(Align & src)
        {
            (*this) = src;
            seq_id_vm = src.seq_id_vm;
            consensus_i_m = src.consensus_i_m;
            return *this;
        }

        vector<ID_T> & get_all_seq_id_ref()
        {
            return seq_id_vm;
        }
        void set_all_seq_id(vector<ID_T> & seq_id_v)
        {
            seq_id_vm = seq_id_v;
        }

        size_t get_consensus_i()
        {
            return consensus_i_m;
        }
        void set_consensus_i(size_t i)
        {
#ifndef NDEBUG
            assert(i<seq_id_vm.size());
#endif
            consensus_i_m = i;
        }
        size_t get_consensus_seq_id()
        {
            return seq_id_vm[consensus_i_m];
        }
        ID_T get_seq_id(size_t i)
        {
#ifndef NDEBUG
            assert(i<seq_id_vm.size());
#endif
            return seq_id_vm[i];
        }
        void set_seq_id(size_t i, ID_T id)
        {
#ifndef NDEBUG
            assert(i<seq_id_vm.size());
#endif
            seq_id_vm[i] = id;
        }
        
        void set_seq_n(size_t n)
        {
            seq_id_vm.resize(n);
            AlignBlock<Pos_T, ID_T>::set_seq_n(n);
        }

    private:
        size_t consensus_i_m;
        vector<ID_T> seq_id_vm;
    };
    // }}}

    /* Alignment Block class 
     *
     * An alignemnt block contains an array of aligned sequence range
     * and statistics, and a list of child sub-block
     */
    // {{{
    template <typename Pos_T, typename ID_T>
    class AlignBlock: public BioObj {
    public:
        typedef Range<Pos_T> Range_T;

        AlignBlock(): BioObj(), parent_m(NULL)
        {
        }

        AlignBlock(size_t n): BioObj(), parent_m(NULL)
        {
			if (n>0) {
				range_vm.resize(n);
				stat_vm.resize(n);
			}
        }

        AlignBlock(const AlignBlock<Pos_T,ID_T> & src): BioObj(src)
        {
            range_vm = src.range_vm;
            stat_vm = src.stat_vm;
        }

        AlignBlock & operator = (const AlignBlock<Pos_T,ID_T> & src)
        {
            (*this) = src;
            range_vm = src.range_vm;
            stat_vm = src.stat_vm;
        }

        size_t get_seq_n()
        {
            return range_vm.size();
        }

        void set_seq_n(size_t n)
        {
            if (n>0) {
                range_vm.resize(n);
                stat_vm.resize(n);
                if (child_vm.size() > 0) {
                    for (size_t i=0; i<child_vm.size(); i++) {
                        child_vm[i]->set_seq_n(n);
                    }
                }
            }
        }

        Range_T & get_range_ref(size_t i)
        {
            return range_vm[i];
        }

        void set_range(size_t i, const Range_T & range)
        {
            range_vm[i] = range;
        }

        AlignStat & get_stat_ref(size_t i)
        {
            return stat_vm[i];
        }
        
        void set_stat(size_t i, const AlignStat & stat)
        {
#ifndef NDEBUG
            assert(i<stat_vm.size());
#endif
            stat_vm[i] = stat;
        }
        
        AlignBlock<Pos_T,ID_T>* get_parent()
        {
            return parent_m;
        }
        void set_parent(AlignBlock* p)
        {
            parent_m[0] = p;
        }

        vector<AlignBlock<Pos_T,ID_T>*> & get_all_child_ref()
        {
            return child_vm;
        }
        void set_all_child(const vector<AlignBlock*> & child_v)
        {
            return child_vm = child_v;
        }

        AlignBlock* get_child(size_t i)
        {
#ifndef NDEBUG
            assert(i<child_vm.size());
#endif
            return child_vm[i];
        }
        void push_back_child(AlignBlock* p)
        {
            child_vm.push_back(p);
        }
        void set_child(size_t i, AlignBlock* p)
        {
#ifndef NDEBUG
            assert(i<child_vm.size());
#endif
            child_vm[i] = p;
        }
        void set_child_n(size_t n)
        {
            child_vm.resize(n);
        }

        void cal_stat_from_bottom()
        {
            for (size_t i=0; i<stat_vm.size(); i++) {
            //  stat_vm[i].clear();
            }
            if (child_vm.size()>0) {
                for (size_t i=0; i<child_vm.size(); i++) {
                    if (child_vm[i]==NULL) continue;
                    child_vm[i]->cal_from_bottom();
                    stat_vm[i] += child_vm[i]->stat_vm[i];
                }
            }
        }

        void clear()
        {
            clear_help();
            if (parent_m != NULL) {
                for (size_t i=0; i<parent_m->child_vm.size(); i++) {
                    if (parent_m->child_vm[i]==this) parent_m->child_vm[i]=NULL;
                }
            }
        }

        ~AlignBlock()
        {
        }

    protected:
        void clear_help()
        {
            for (size_t i=0; i<child_vm.size(); i++) {
                child_vm[i].clear_help();
            }
            for (size_t i=0; i<range_vm.size(); i++) {
                range_vm[i].clear();
            //  stat_vm[i].clear();
            }
            for (size_t i=0; i<child_vm.size(); i++) {
                if (child_vm[i]==NULL) continue;
                delete child_vm[i];
            }
            child_vm.clear();
        }

        vector<Range_T> range_vm;
        vector<AlignStat> stat_vm;
        vector<AlignBlock<Pos_T,ID_T>*> child_vm;
        AlignBlock<Pos_T,ID_T>* parent_m;
    };
    // }}}

    /* class AlignStat
     * statistics of an alignment 
     */
    // {{{
    class AlignStat {
    public:
        AlignStat() 
        {
            len_m      = 0;
            match_m    = 0;
            mismatch_m = 0;
            gap_n_m    = 0;
            gap_len_m  = 0;
            rep_match_m= 0;
            N_n_m      = 0;
            score_m    = 0.0;
        }

        AlignStat(const AlignStat & src) 
        {
            len_m      = src.len_m;
            match_m    = src.match_m;
            mismatch_m = src.mismatch_m;
            gap_n_m= src.gap_n_m;
            gap_len_m  = src.gap_len_m;
            rep_match_m= src.rep_match_m;
            N_n_m      = src.N_n_m;
            score_m    = src.score_m;
        }

        AlignStat & operator =(const AlignStat & src) 
        {
            len_m      = src.len_m;
            match_m    = src.match_m;
            mismatch_m = src.mismatch_m;
            gap_n_m= src.gap_n_m;
            gap_len_m  = src.gap_len_m;
            rep_match_m= src.rep_match_m;
            N_n_m      = src.N_n_m;
            score_m    = src.score_m;
        }

        void clear()
        {
            len_m      = 0;
            match_m    = 0;
            mismatch_m = 0;
            gap_n_m    = 0;
            gap_len_m  = 0;
            rep_match_m= 0;
            N_n_m      = 0;
            score_m    = 0.0;
        }

        bool operator ==(const AlignStat & src) {
            return (
                len_m      == src.len_m
                &&
                match_m    == src.match_m
                &&
                mismatch_m == src.mismatch_m
                &&
                gap_n_m== src.gap_n_m
                &&
                gap_len_m  == src.gap_len_m
                &&
                rep_match_m == src.rep_match_m
                &&
                N_n_m      == src.N_n_m
                &&
                score_m    == src.score_m
            );
        }

        size_t get_len()
        {
            return len_m;    
        }
        void set_len(size_t len)
        {
            len_m = len;    
        }

        size_t get_match()
        {
            return match_m;    
        }
        void set_match(size_t match)
        {
            match_m = match;    
        }

        size_t get_mismatch()
        {
            return mismatch_m;    
        }
        void set_mismatch(size_t mismatch)
        {
            mismatch_m = mismatch;    
        }

        size_t get_gap_n()
        {
            return gap_n_m;    
        }
        void set_gap_n(size_t gap_n)
        {
            gap_n_m = gap_n;
        }

        size_t get_gap_len()
        {
            return gap_len_m;    
        }
        void set_gap_len(size_t gap_len)
        {
            gap_len_m = gap_len;
        }

        double get_score()
        {
            return score_m;    
        }
        void set_score(double score)
        {
            score_m = score;
        }

        size_t get_rep_match()
        {
            return rep_match_m;    
        }
        void set_rep_match(size_t rep_match)
        {
            rep_match_m = rep_match;    
        }

        size_t get_N_n()
        {
            return N_n_m;    
        }
        void set_N_n(size_t N_n)
        {
            N_n_m = N_n;    
        }

        double get_identity()
        {
            return (double)match_m/len_m;
        }

        double get_coverage(size_t seq_len)
        {
            return (double)(match_m+mismatch_m)/seq_len;
        }

        AlignStat & operator +=(const AlignStat & a)
        {
            len_m      += a.len_m;
            match_m    += a.match_m;
            mismatch_m += a.mismatch_m;
            gap_n_m    += a.gap_n_m;
            gap_len_m  += a.gap_len_m;
            rep_match_m += a.rep_match_m;
            N_n_m      += a.N_n_m;
            score_m    += a.score_m;
        }

        AlignStat operator +(AlignStat & a)
        {
            AlignStat out;
            out.len_m      = len_m      + a.len_m;
            out.match_m    = match_m    + a.match_m;
            out.mismatch_m = mismatch_m + a.mismatch_m;
            out.gap_n_m    = gap_n_m    + a.gap_n_m;
            out.gap_len_m  = gap_len_m  + a.gap_len_m;
            out.rep_match_m = rep_match_m + a.rep_match_m;
            out.N_n_m      = N_n_m      + a.N_n_m;
            out.score_m    = score_m    + a.score_m;
            return out;
        }

        ~AlignStat()
        {}
    protected:
        size_t len_m;
        size_t match_m;
        size_t mismatch_m;
        size_t gap_n_m;
        size_t gap_len_m;
        size_t rep_match_m;
        size_t N_n_m;
        double score_m;
    };
    // }}}

    /* class IO_BLAST
     */
    // {{{
    template <typename Pos_T=int64_t, typename ID_T=int64_t>
    class IOBLAST {
        typedef Align<Pos_T, ID_T> Align_T;
    public:
        static int get_one(istream & in, Align<Pos_T, ID_T> & align, map<string, ID_T> & name_to_id, vector<string> & name_v)
        // {{{
        {
            string line;
            if (in.good()) {
                getline(in, line);
                return get_one(line, align, name_to_id, name_v);
            } else {
                return -1;
            }
        }
        // }}}

        static int get_one(string & line, Align<Pos_T, ID_T> & align, map<string, ID_T> & name_to_id, vector<string> & name_v) 
        // {{{
        {
            int i=0;
            int n=0;
            int j;
            string qname, tname;
            size_t align_len, match, mismatch, gap_open;
            ID_T tbegin, tend, qbegin, qend;
            float iden_perc;
            double bit, evalue;
            ID_T qid, tid;
            // qname
            n++;
            i = Token::get_next(line, i, "\t", qname);
            if (i==string::npos) {
                return n;
            }
            if (name_to_id.find(qname)==name_to_id.end()) {
                qid = name_v.size();
                name_v.push_back(qname);
                name_to_id[qname] = qid;
            } else {
                qid = name_to_id[qname];
            }
            // tname
            n++;
            i = Token::get_next(line, i, "\t", tname);
            if (i==string::npos) {
                return n;
            }
            if (name_to_id.find(tname)==name_to_id.end()) {
                tid = name_v.size();
                name_v.push_back(tname);
                name_to_id[tname] = tid;
            } else {
                tid = name_to_id[tname];
            }

            string s;
            // iden perc
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            iden_perc = atof(s.c_str());
            // align length
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            align_len = atol(s.c_str());
            // mismatch
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            mismatch = atol(s.c_str());
            // gap_open 
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            gap_open = atol(s.c_str());
            // query begin
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            qbegin = atol(s.c_str());
            // query end
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            qend = atol(s.c_str());
            int qs = '?';
            if (qbegin<qend) qs = '+';
            else {
                int t = qbegin;
                qbegin = qend;
                qend = t;
                qs = '-';
            }
            // target begin
            char ts = '?';
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            tbegin = atol(s.c_str());
            // target end
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            tend = atol(s.c_str());
            if (tbegin<tend) ts = '+';
            else {
                Pos_T t = tbegin;
                tbegin = tend;
                tend = t;
                ts = '-';
            }
            qbegin--;
            tbegin--;

            // E value
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                return n;
            }
            evalue = atof(s.c_str());
            // bit
            n++;
            i = Token::get_next(line, i, "\t", s);
            if (i==string::npos) {
                j = line.size();
            }
            bit = atof(s.c_str());

            /// set align
            // {{{
            align.set_seq_n(2);
            align.set_seq_id(0, tid);
            align.set_seq_id(1, qid);

            size_t qlen = qend - qbegin;
            size_t tlen = tend - tbegin;
            match = iden_perc/100 * align_len + 0.5;
            size_t qgap_len = align_len - qlen;
            size_t tgap_len = align_len - tlen;

            Range<Pos_T> & trange = align.get_range_ref(0);
            trange.set(tbegin, tend, ts);
            AlignStat & tstat = align.get_stat_ref(0);
            tstat.set_match(tlen);
            tstat.set_gap_len(tgap_len);
            tstat.set_score(evalue);

            Range<Pos_T> & qrange = align.get_range_ref(1);
            qrange.set(qbegin, qend, qs);
            AlignStat & qstat = align.get_stat_ref(1);
            qstat.set_len(align_len);
            qstat.set_match(match);
            qstat.set_mismatch(mismatch);
            qstat.set_gap_n(gap_open);
            qstat.set_gap_len(qgap_len);
            qstat.set_score(bit);
            // }}}

            return 0;
        }
        // }}}

        static int put_one(ostream & out, Align<Pos_T, ID_T> & align, vector<string> & name_v)
        // {{{
        {
            string line;
            if (out.good()) {
                put_one(line, align, name_v);
                out << line;
            }
            if (out.good()) return 0;
            else return -1;
        }
        // }}}

        static int put_one(string & line, Align<Pos_T, ID_T> & align, vector<string> & name_v)
        // {{{
        {
#ifndef NDEBUG
            assert(align.get_seq_n()==2);
#endif
            stringstream ss;
            ID_T qid = align.get_seq_id(1);
            ID_T tid = align.get_seq_id(0);
#ifndef NDEBUG
            assert(qid<name_v.size());
            assert(tid<name_v.size());
#endif
            ss << name_v[qid] << "\t" << name_v[tid];
            Range<Pos_T> & trange = align.get_range_ref(0);
            Range<Pos_T> & qrange = align.get_range_ref(1);
            AlignStat & tstat = align.get_stat_ref(0);
            AlignStat & qstat = align.get_stat_ref(1);
            float iden_perc = qstat.get_identity()*100;
            ss << "\t" << iden_perc << "\t" << qstat.get_len() << "\t" << qstat.get_mismatch() << "\t" << qstat.get_gap_n();
            Pos_T qbegin = qrange.get_begin();
            Pos_T tbegin = trange.get_begin();
            Pos_T qend = qrange.get_end();
            Pos_T tend = trange.get_end();
            char qstrand = qrange.get_strand();
            char tstrand = trange.get_strand();
            if (qstrand!='-') {
                ss << "\t" << qbegin+1 << "\t" << qend;
            } else {
                ss << "\t" << qend << "\t" << qbegin+1;
            }
            if (tstrand!='-') {
                ss << "\t" << tbegin+1 << "\t" << tend;
            } else {
                ss << "\t" << tend << "\t" << tbegin+1;
            }
            double evalue = tstat.score_m;
            float bit= qstat.score_m;
            ss << "\t" << evalue << "\t" << bit;

            line = ss.str();
            return 0;
        }
        // }}}

        static int skip_one(istream & in)
        // {{{
        {
            string line;
            if (in.good()) {
                getline(in, line);
            } else {
                return -1;
            }
        }
        // }}}

        static int skip(istream & in, size_t n)
        // {{{
        {
            string line;
            while (in.good() && n>0) {
                getline(in, line);
                n--;
            }
            if (in.good()) return n;
            else return -1;
        }
        // }}}
    };
    // }}}

    /* class IOBLAT<Pos_T, ID_T>
     */
    template <typename Pos_T, typename ID_T>
    class IOBLAT {
        typedef Align<Pos_T, ID_T> Align_T;
        typedef AlignBlock<Pos_T, ID_T> AlignBlock_T;
        typedef typename Align_T::Range_T Range_T;
    public:
        static int get_one(istream & in, Align<Pos_T, ID_T> & align, map<string, ID_T> & name_to_id, vector<string> & name_v, vector<Pos_T> & seq_len_v)
        // {{{
        {
            string line;
            if (in.good()) {
                getline(in, line);
                return get_one(line, align, name_to_id, name_v, seq_len_v);
            } else {
                return -1;
            }
        }
        // }}}

        static int get_one(string & line, Align<Pos_T, ID_T> & align, map<string, ID_T> & name_to_id, vector<string> & name_v, vector<Pos_T> & seq_len_v) 
        // {{{
        {
            int i=0;
            int n=0;
            int j;
            string qname, tname;
            size_t match, rep_match, mismatch, N_n;
            size_t qins_n, tins_n;
            size_t block_n;
            Pos_T align_len, qins_len, tins_len;
            Pos_T tbegin, tend, qbegin, qend;
            Pos_T tlen, qlen;
            ID_T qid, tid;
            const size_t nc = 21;
            char qstrand;
            const string col_name_v[nc] = {
                "match", "mismatch", "rep_match", "N_n", "qins_n", 
                "qins_len", "tins_n", "tins_len", "qstrand", "qname",
                "qlen", "qbegin", "qend", "tname", "tlen",
                "tbegin", "tend", "block_n", "block_lens", "qbegins",
                "tbegins"
            };
            const string col_type_v[nc] = {
                "size_t", "size_t", "size_t", "size_t", "size_t", 
                "size_t", "size_t", "size_t", "size_t", "string", 
                "size_t", "size_t", "size_t", "string", "size_t", 
                "size_t", "size_t", "size_t", "string", "string", 
                "string"
            };
            string s;
            vector<Pos_T> align_len_v, qbegin_v, tbegin_v;
            for (n=0; n<nc; n++) {
                if (i==string::npos) {
                    return n;
                }
                i = Token::get_next(line, i, "\t", s);
                switch (n) {
                    // match
                    case 0: match = atol(s.c_str()); break;
                    // mismatch
                    case 1: mismatch = atol(s.c_str()); break;
                    // rep_match
                    case 2: rep_match = atol(s.c_str()); break;
                    // N_n
                    case 3: N_n = atol(s.c_str()); break;
                    // qins_n
                    case 4: qins_n = atol(s.c_str()); break;
                    // qins_len
                    case 5: qins_len = atol(s.c_str()); break;
                    // tins_n
                    case 6: tins_n = atol(s.c_str()); break;
                    // tins_len
                    case 7: tins_len = atol(s.c_str()); break;
                    // qstrand
                    case 8: qstrand = s[0]; break;
                    // qname
                    case 9: qname = s; break;
                    case 10: qlen = atol(s.c_str()); break;
                    case 11: qbegin = atol(s.c_str()); break;
                    case 12: qend = atol(s.c_str()); break;
                    // tname
                    case 13: tname = s; break;
                    case 14: tlen = atol(s.c_str()); break;
                    case 15: tbegin = atol(s.c_str()); break;
                    case 16: tend = atol(s.c_str()); break;
                    // block count
                    case 17: block_n = atol(s.c_str()); break;
                    // block size
                    case 18: {
                        j=0;
                        for (int m=0; m<block_n; m++) {
							string s1;
                            j = Token::get_next(s, j, ",", s1);
							if (s1.empty()) break;
                            align_len_v.push_back(atol(s1.c_str()));
                        }
                        break;
                    }
                    // block qbegin
                    case 19: {
                        j=0;
                        for (int m=0; m<block_n; m++) {
							string s1;
                            j = Token::get_next(s, j, ",", s1);
							if (s1.empty()) break;
                            qbegin_v.push_back(atol(s1.c_str()));
                        }
                        break;
                    }
                    // block qbegin
                    case 20: {
                        j=0;
                        for (int m=0; m<block_n; m++) {
							string s1;
                            j = Token::get_next(s, j, ",", s1);
							if (s1.empty()) break;
                            tbegin_v.push_back(atol(s1.c_str()));
                        }
                        break;
                    }
                }
				s.clear();
            }
            if (name_to_id.find(qname)==name_to_id.end()) {
                qid = name_v.size();
                name_v.push_back(qname);
                seq_len_v.push_back(qlen);
                name_to_id[qname] = qid;
            } else {
                qid = name_to_id[qname];
#ifndef NDEBUG
                assert(qlen==seq_len_v[qid]);
#endif
            }
            if (name_to_id.find(tname)==name_to_id.end()) {
                tid = name_v.size();
                name_v.push_back(tname);
                seq_len_v.push_back(tlen);
                name_to_id[tname] = tid;
            } else {
                tid = name_to_id[tname];
#ifndef NDEBUG
                assert(tlen==seq_len_v[tid]);
#endif
            }

            size_t qgap_n=0, tgap_n=0;
            Pos_T qgap_len=0, tgap_len=0;

            /// set align
            // {{{
            align.set_seq_n(2);
            align.set_seq_id(0, tid);
            align.set_seq_id(1, qid);

            Range_T & trange = align.get_range_ref(0);
            Range_T & qrange = align.get_range_ref(1);
            AlignStat & tstat = align.get_stat_ref(0);
            AlignStat & qstat = align.get_stat_ref(1);

            for (i=0; i<block_n; i++) {
                Pos_T qbegin = qbegin_v[i];
                Pos_T tbegin = tbegin_v[i];
                Pos_T align_len = align_len_v[i];
                Pos_T qend = qbegin + align_len;
                Pos_T tend = tbegin + align_len;
                AlignBlock_T *p = new AlignBlock_T(2);
                Range_T & trange = p->get_range_ref(0);
                AlignStat & tstat = p->get_stat_ref(0);
                Range_T & qrange = p->get_range_ref(1);
                AlignStat & qstat = p->get_stat_ref(1);
                qrange.set(qbegin, qend, qstrand);
                trange.set(tbegin, tend, '+');
                qstat.set_len(align_len);
                tstat.set_len(align_len);
                align.push_back_child(p);
            }

            trange.set(tbegin, tend, '+');
            tstat.set_len(align_len);
            tstat.set_match(match);
            tstat.set_rep_match(rep_match);
            tstat.set_mismatch(mismatch);
            tstat.set_gap_n(qins_n);
            tstat.set_gap_len(qins_len);

            qrange.set(qbegin, qend, qstrand);
            qstat.set_len(align_len);
            qstat.set_match(match);
            tstat.set_rep_match(rep_match);
            qstat.set_mismatch(mismatch);
            qstat.set_gap_n(tins_n);
            qstat.set_gap_len(tins_len);
            qstat.set_N_n(N_n);
            // }}}

            return 0;
        }
        // }}}

        static int put_one(ostream & out, Align<Pos_T, ID_T> & align, vector<string> & name_v, vector<Pos_T> & seq_len_v)
        // {{{
        {
            string line;
            if (out.good()) {
                return put_one(line, align, name_v, seq_len_v);
                out << line;
            } else {
                return -1;
            }
        }
        // }}}

        static int put_one(string & line, Align<Pos_T, ID_T> & align, vector<string> & name_v, vector<Pos_T> seq_len_v)
        // {{{
        {
            int i, j, n;
            size_t qins_n, tins_n;
            size_t match, rep_match, mismatch, N_n;
            size_t block_n;
            Pos_T align_len, qins_len, tins_len;
            Pos_T tbegin, tend, qbegin, qend;
            ID_T qid, tid;
            char qstrand;
            size_t nc = 21;

            stringstream ss;

            Range_T & trange = align.get_range_ref(0);
            Range_T & qrange = align.get_range_ref(1);
            AlignStat & tstat = align.get_stat_ref(0);
            AlignStat & qstat = align.get_stat_ref(1);

            tbegin = trange.get_begin();
            tend = trange.get_end();
            qbegin = qrange.get_begin();
            qend = qrange.get_end();
            qstrand = qrange.get_strand();

            match = qstat.get_match();
            rep_match = qstat.get_rep_match();
            mismatch = qstat.get_mismatch();
            N_n = qstat.get_N_n();
            align_len = qstat.get_len();
            qins_n = tstat.get_gap_n();
            qins_len = tstat.get_gap_len();
            tins_n = qstat.get_gap_n();
            tins_len = qstat.get_gap_len();

            ss << match << "\t" << mismatch << "\t" << rep_match << "\t" << N_n;
            ss << "\t" << qins_n << "\t" << qins_len << "\t" << tins_n << "\t" << tins_len;
            ss << "\t" << qstrand;
            ss << "\t" << name_v[qid] << "\t" << seq_len_v[qid];
            ss << "\t" << qbegin << "\t" << qend;
            ss << "\t" << name_v[tid] << "\t" << seq_len_v[tid];
            ss << "\t" << tbegin << "\t" << tend;
            
            block_n = align.get_child_n();
            ss << "\t" << block_n << "\t";
            for (i=0; i<block_n; i++) {
                if (i>0) ss << ",";
                ss << align.get_child(i)->get_stat_ref(1).get_len();
            }
            ss << "\t";
            for (i=0; i<block_n; i++) {
                if (i>0) ss << ",";
                ss << align.get_child(i)->get_range_ref(1).get_begin();
            }
            ss << "\t";
            for (i=0; i<block_n; i++) {
                if (i>0) ss << ",";
                ss << align.get_child(i)->get_range_ref(0).get_begin();
            }
            return 0;
        }
        // }}}
    };
};

#endif
