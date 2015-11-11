#include "czl_bam.h"

/// class Bam
// {{{
/*
 *
 */
int Bam::split_by_chr_pos(string & bam_file, int is_sorted_by_coord, string & out_prefix, int length, int step, string * out_bed_file)
{
    /// open a bam file
    bamFile bam_fp;
    bam_fp = bam_open(bam_file.c_str(), "rb");
    if (bam_fp == NULL) {
        cerr << "Fail to open BAM file " << bam_file << endl;
        return 1;
    }
//  bam_index_t *bam_index = bam_index_load(bam_file.c_str());
//  if (bam_index == NULL) {
//      cerr << "Fail to open BAM index " << bam_file << endl;
//      return 1;
//  }
    bam_header_t *bam_header_p = bam_header_read(bam_fp);
    if (bam_header_p==NULL) {
        cerr << "Fail get header from BAM FILE " << bam_file << endl;
        return 2;
    }

    string file = out_prefix+"unmap.bam";
    string dir;
    bamFile unmap_fp = bam_open(file.c_str(), "wb");
    if (unmap_fp==NULL) {
        cerr << "Fail to open BAM FILE " << file << " for write." << endl;
        return 1;
    }
    bam_header_write(unmap_fp, bam_header_p);

    ofstream fout_bed;
    if (out_bed_file!=NULL) {
        fout_bed.open( out_bed_file->c_str() );
        if (!fout_bed.is_open()) {
            cerr << "Fail to open FILE " << *out_bed_file << endl;
            return 3;
        }
    }

    const int32_t nt = bam_header_p->n_targets;
    bam1_t *bam;
    uint32_t tid, begin, end;
    list<bam1_t*> bams;
    if ( is_sorted_by_coord ) {
        bamFile fp;
        bam = bam_init1();
        begin = 0;
        end = begin+length;
        tid = -1;
        int r = 0;
        while ( bam_read1(bam_fp, bam)>0 ) {
            if ( (bam->core.flag & BAM_FUNMAP)!=0 ) { // unmaped
                bam_write1(unmap_fp, bam);
            } else {
                uint32_t b = bam->core.pos;
            //  uint32_t e = bam_calend(&bam->core, bam1_cigar(bam));
                if ( bam->core.tid != tid || b > end) {
                    if (r!=0) { // if not first bam, close previous bam file
                        bam_close(fp);
                        if (tid!=bam->core.tid) {
                            tid = bam->core.tid;
                            begin = 0;
                            while ( !bams.empty() ) {
                                bam_destroy1(bams.front());
                                bams.pop_front();
                            }
                        } else {
                            begin += step;
                        }
                        end = begin + length;
                    } else {
                        tid = bam->core.tid;
                        r=1;
                    }
                    char *id = bam_header_p->target_name[tid];
                    stringstream ss;
                    ss << out_prefix << id << "_" << begin << "_" << end;
                    fp = bam_open(ss.str().c_str(), "wb");
                    if (fp==NULL) {
                        cerr << "Fail to open file " << ss.str() << endl;
                    }
                    bam_header_write(fp, bam_header_p);
                    if ( out_bed_file!=NULL ) {
                        fout_bed << id << "\t" << begin << "\t" << end << "\n";
                    }

					BOOST_FOREACH( bam1_t *bam1, bams) {
						if (bam1->core.pos<end) {
							bam_write1(fp, bams.front());
						}
					}
                    while ( !bams.empty() && bams.front()->core.pos < begin+step ) {
                        bam_destroy1(bams.front());
                        bams.pop_front();
                    }
                }
                bam_write1(fp, bam);
                if ( b >= begin+step) {
                    bams.push_back(bam);
                    bam = bam_init1();
                }
            }
        }
        if (!bams.empty()){
            while ( !bams.empty() && bams.front()->core.pos < end ) {
                bam_write1(fp, bams.front());
                bam_destroy1(bams.front());
                bams.pop_front();
            }
			bam_close(fp);
            while ( !bams.empty() ) {
                bam_destroy1(bams.front());
                bams.pop_front();
            }
        }
        bam_destroy1(bam);
    } else {
        cerr << "TO Develop" << endl;
    }
    if ( out_bed_file!=NULL ) {
        fout_bed.close();
    }
    bam_close(unmap_fp);

    return 0;
}

int Bam::split_by_chr(string & bam_file, string & out_prefix, vector<string> * out_chroms) 
{
    /// open a bam file
    bamFile bam_sn_fh;
    bam_sn_fh = bam_open(bam_file.c_str(), "rb");
    if (bam_sn_fh == NULL) {
        cerr << "Fail to open BAM file " << bam_file << endl;
        return 1;
    }
//  bam_index_t *bam_index = bam_index_load(bam_file.c_str());
//  if (bam_index == NULL) {
//      cerr << "Fail to open BAM index " << bam_file << endl;
//      return 1;
//  }
    bam_header_t *bam_header_p = bam_header_read(bam_sn_fh);

    string file = out_prefix+"unmap.bam";
    string dir;
    bamFile unmap_fp = bam_open(file.c_str(), "wb");
    bam_header_write(unmap_fp, bam_header_p);

    const int32_t nt = bam_header_p->n_targets;
    const int32_t n=100000;
    int32_t m=0, m1=0;
    vector<bam1_t*> bam1s[nt];
    vector<int32_t> tids;
    int64_t bam1_nums[nt];
    for (int32_t i=0; i<nt; i++) { bam1_nums[i] = 0; }
    string bam_files[nt];
    for (int32_t i=0; i<nt; i++) {
        bam_files[i] = out_prefix+bam_header_p->target_name[i];
        if (out_chroms!=NULL) {
            out_chroms->push_back(bam_header_p->target_name[i]);
        }
        dir = out_prefix+"tmp_"+boost::lexical_cast<string>(i);
        boost::filesystem::path p(dir);
        if (!boost::filesystem::exists(p)) {
            boost::filesystem::create_directory(p);
            if (!boost::filesystem::exists(p)) {
                cerr << "Can't create DIR " << dir;
                return 1;
            }
        }
    }
    if (out_chroms!=NULL) {
        out_chroms->push_back("unmap");
    }

    /// read a bam alignment
    map<int32_t, bamFile> tid_map_bamfp;
    bam1_t *bam1 = bam_init1();
    while ( bam_read1(bam_sn_fh, bam1) > 0 ) {
        if ( (bam1->core.flag & BAM_FUNMAP)!=0 ) {
            bam_write1(unmap_fp, bam1);
        } else {
            int32_t tid = bam1->core.tid;
            if (bam1s[tid].size()==0) tids.push_back(tid);
            bam1s[tid].push_back(bam1);
            m++;
            if ( m == n) {
                for (int32_t j=0; j<tids.size(); j++) {
                    tid = tids[j];
                    file = out_prefix+"tmp_"+boost::lexical_cast<string>(tid)+"/"+boost::lexical_cast<string>(bam1_nums[tid]++);
                    bamFile fp = bam_open(file.c_str(), "wb");
                //  bam_header_write(fp, bam_header_p);
                    BOOST_FOREACH(bam1, bam1s[tid]) {
                        bam_write1(fp, bam1);
                        bam_destroy1(bam1);
                    }
                    bam1s[tid].clear();
                    bam_close(fp);
                }
                tids.clear();
                m=0;
            }
            /*
            int32_t i = bam1_nums[tid]++;
            bam1s[tid][i] = bam1;
            if (bam1_nums[tid]==n) {
                bamFile fp = bam_open(bam_files[tid].c_str(), "ab");
                if (fp==NULL) {
                    cerr << "fail to open BAM file " << bam_files[tid] << endl;
                    return 1;
                }
                for (int32_t j=0; j<n; j++) {
                    bam_write1(fp, bam1s[tid][j]);
                    bam_destroy1(bam1s[tid][j]);
                }
                bam1_nums[tid] = 0;
                bam_close(fp);
            }
            */
            bam1 = bam_init1();
        }
    }
    bam_destroy1(bam1);

    for (int32_t j=0; j<tids.size(); j++) {
        int32_t tid = tids[j];
        file = out_prefix+"tmp_"+boost::lexical_cast<string>(tid)+"/"+boost::lexical_cast<string>(bam1_nums[tid]++);
        bamFile fp = bam_open(file.c_str(), "wb");
    //  bam_header_write(fp, bam_header_p);
        BOOST_FOREACH(bam1, bam1s[tid]) {
            bam_write1(fp, bam1);
            bam_destroy1(bam1);
        }
        bam1s[tid].clear();
        bam_close(fp);
    }
    tids.clear();
    /*
    for (int32_t tid=0; tid<nt; tid++) {
        if (bam1_nums[tid]>0) {
            bamFile fp = bam_open(bam_files[tid].c_str(), "wb");
            bam_seek(fp, 0, SEEK_END);
            for (int32_t j=0; j<bam1_nums[tid]; j++) {
                bam_write1(fp, bam1s[tid][j]);
                bam_destroy1(bam1s[tid][j]);
            }
            bam1_nums[tid] = 0;
            bam_close(fp);
        }
    }
    */
    bam_close(bam_sn_fh);
    bam_close(unmap_fp);

    bam1 = bam_init1();
    for (int32_t tid=0; tid<nt; tid++) {
        if (bam1_nums[tid]>0) {
            bamFile fpw = bam_open(bam_files[tid].c_str(), "wb");
            bam_header_write(fpw, bam_header_p);
            for (int64_t i=0; i<bam1_nums[tid]; i++) {
                file = out_prefix+"tmp_"+boost::lexical_cast<string>(tid)+"/"+boost::lexical_cast<string>(i);
                bamFile fpr = bam_open(file.c_str(), "rb");
                m=0;
                while ( bam_read1(fpr, bam1)>0 ) {
                    bam_write1(fpw, bam1);
                    bam1 = bam_init1();
                    m++;
                }
                bam_close(fpr);
                boost::filesystem::path p(file);
                boost::filesystem::remove(p);
            }
            bam_close(fpw);
            bam1_nums[tid]=0;
        }
        dir = out_prefix+"tmp_"+boost::lexical_cast<string>(tid)+"/";
        boost::filesystem::path p(dir);
        boost::filesystem::remove(p);
    }
    bam_destroy1(bam1);

    return 0;
}

struct BamTmpA {
    int k;  // current pos in name
    int32_t i:24, flag:8;  // i: split file id
    int parent;
    int read_num;  // read number in file
    string name;
};
/**
 *  @brief sort bam by names
 */
int Bam::sort_by_name(string & bam_file, string & out_bam_file, string & tmp_dir, int buf_n)
{
    return sort_by_name(bam_file.c_str(), out_bam_file.c_str(), tmp_dir.c_str(), buf_n);
}

int Bam::sort_by_name(const char *bam_file, const char *out_bam_file, const char *tmp_dir, int buf_n)
{
    int k = 0;
    string tmp_dir_s(tmp_dir);
    tmp_dir_s+="/t";
    string null_s(tmp_dir);
    null_s += "/n";

    vector<BamTmpA*> stack(1);
    stack[0] = new BamTmpA;
    stack[0]->k = 0;
    stack[0]->i = 0;
    stack[0]->flag = 0;
    stack[0]->parent = -1;
    stack[0]->name = tmp_dir_s;
    vector<bam1_t*> bams(buf_n);
    vector<bam1_t*> out_bams;
    bam1_t *bam;
    string file(bam_file);
    int n=0;
    const int K=256;

    bamFile fp;
    fp = bam_open(bam_file, "rb");
    bam_header_t *bam_header_p = bam_header_read(fp);
    bam_close(fp);

    int m=0;
    int k0=0;
    char nz_chr[K];
    int nnz=0;
    int count[K];
    bamFile bam_fps[K];
    int read_num=0;
    string name;
    while (1) {
        for (int i=0; i<K; i++) bam_fps[i]=NULL;
        for (int i=0; i<K; i++) count[i]=0;
        fp = bam_open(bam_file, "rb");
        bam_header_destroy(bam_header_read(fp));
        read_num=0;
        bam = bam_init1();
        char c;
        bam_fps[0] = bam_open(null_s.c_str(), "wb");
        bam_header_write(bam_fps[0], bam_header_p);
        m=0;
        nnz=0;
        while ( bam_read1(fp, bam)>0 ) {
            if (m>=buf_n) {
                for (int j=0; j<m; j++) {
                    bam1_t* bam1 = bams[j];
                    char *qname = bam1_qname(bam1);
                    c = qname[k0];
                    if ( bam_fps[c]==NULL ) {
                        name.push_back(c);
                        file = tmp_dir_s + name;
                        bam_fps[c] = bam_open(file.c_str(), "wb");
                        bam_header_write(bam_fps[c], bam_header_p);
                        name.erase(name.end()-1);
                        nnz++;
                    }
                    bam_write1(bam_fps[c], bam1);
                    bam_destroy1(bam1);
                    count[c]++;
                }
                m=0;
            }
            bams[m++] = bam;
            bam = bam_init1();
            read_num++;
        }
        bam_destroy1(bam);
        bam_close(fp);
        bam_close(bam_fps[0]);

        if ( read_num <= buf_n ) {
            boost::filesystem::remove(null_s);

            sort_by_name(bams.begin(), bams.begin()+read_num, k0);
            bamFile fpw = bam_open(out_bam_file, "wb");
            bam_header_write(fpw, bam_header_p);
            for (int j=0; j<read_num; j++) {
                bam1_t* bam1 = bams[j];
                bam_write1(fpw, bam1);
                bam_destroy1(bam1);
            }
            bam_close(fpw);
            return 0;
        } else {
            for (int j=0; j<m; j++) {
                bam1_t* bam1 = bams[j];
                char *qname = bam1_qname(bam1);
                c = qname[k0];
                if ( bam_fps[c]==NULL ) {
                    name.push_back(c);
                    file = tmp_dir_s + name;
                    bam_fps[c] = bam_open(file.c_str(), "wb");
                    bam_header_write(bam_fps[c], bam_header_p);
                    name.erase(name.end()-1);
                    nnz++;
                }
                bam_write1(bam_fps[c], bam1);
                bam_destroy1(bam1);
                count[c]++;
            }
            if (nnz==0) {
                boost::filesystem::remove(null_s);

                fp = bam_open(bam_file, "rb");
                bam_header_destroy(bam_header_read(fp));
                bamFile fpw = bam_open(out_bam_file, "wb");
                bam_header_write(fpw, bam_header_p);
                bam = bam_init1();
                while ( bam_read1(fp, bam)>0 ) {
                    bam_write1(fpw, bam);
                }
                bam_destroy1(bam);
                bam_close(fpw);
                bam_close(fp);
                return 0;
            } else if (nnz==1 && count[0]==0) {
                boost::filesystem::remove(null_s);
                name.push_back(c);
                file = tmp_dir_s+name;
                boost::filesystem::remove(file);
                k0++;
            } else {
                stack[0]->name = name;
                stack[0]->k = k0;
                file = tmp_dir_s+name;
                boost::filesystem::rename(null_s, file);
                int j=0;
                for (int16_t c=1; c<K; c++) {
                    if (count[c]>0) {
                        BamTmpA *a = new BamTmpA;
                        a->k = k0+1;
                        a->i = j++;
                        a->flag = 0;
                        a->read_num = count[c];
                        a->parent = 0;
                        a->name = name;
                        a->name.push_back(c);
                        stack.push_back(a);
                        bam_close(bam_fps[c]);
                        n++;
                    }
                }
                break;
            }
        }
    }
    while (!stack.empty()) {
        BamTmpA *a0 = stack[n];
        if ( (a0->flag & 0x01) != 0 ) { // if finished
            if ( n==0 ) break;
            if ( a0->i==0 ) { //   if first one in current level
            // merge to parent, set parent finished, pop all child
                int pn = a0->parent;
                file = tmp_dir_s+stack[pn]->name;
                boost::filesystem::rename(file, null_s);
                bamFile fpw = bam_open(file.c_str(), "wb");
                bam_header_write(fpw, bam_header_p);
                /// write 'null_s'
                {
                    fp = bam_open(null_s.c_str(), "rb");
                    bam_header_destroy(bam_header_read(fp));
                    bam = bam_init1();
                    while ( bam_read1(fp, bam)>0 ) {
                        bam_write1(fpw, bam);
                    }
                    bam_destroy1(bam);
                    bam_close(fp);
                    boost::filesystem::remove(null_s);
                }
                for (int j=n; j<stack.size(); j++) {
                    file = tmp_dir_s + stack[j]->name;
                    fp = bam_open(file.c_str(), "rb");
                    bam_header_destroy(bam_header_read(fp));
                    bam = bam_init1();
                    while ( bam_read1(fp, bam)>0 ) {
                        bam_write1(fpw, bam);
                    }
                    bam_destroy1(bam);
                    bam_close(fp);
                    boost::filesystem::remove(file);
                }
                bam_close(fpw);
                stack[pn]->flag |= 0x01; // set finished
                for (int j = stack.size()-1; j>=n; j--) {
                    delete stack[j];
                    stack.pop_back();
                }
                n = pn;
            } else {
                n--;
            }
        } else {
            if (a0->read_num <= buf_n) {
            // sort and write back, set finished
                file = tmp_dir_s+a0->name;
                fp = bam_open(file.c_str(), "rb");
                bam_header_destroy(bam_header_read(fp));
                bam = bam_init1();
                m=0;
                while ( bam_read1(fp, bam)>0 ) {
                    bams[m++] = bam;
                    bam = bam_init1();
                }
                bam_close(fp);

                sort_by_name(bams.begin(), bams.begin()+m, a0->k);

                fp = bam_open(file.c_str(), "wb");
                bam_header_write(fp, bam_header_p);
                for (int j=0; j<m; j++) {
                    bam_write1(fp, bams[j]);
                    bam_destroy1(bams[j]);
                }
                bam_destroy1(bam);
                bam_close(fp);
                a0->flag |= 0x01;
            } else {
                // split and add to stack
                string name = a0->name;
                while (1) {
                    for (int i=0; i<K; i++) bam_fps[i]=NULL;
                    for (int i=0; i<K; i++) count[i]=0;
                    bam_fps[0] = bam_open(null_s.c_str(), "wb");
                    bam_header_write(bam_fps[0], bam_header_p);
                    read_num=0;
                    nnz=0;
                    file = tmp_dir_s+a0->name;
                    fp = bam_open(file.c_str(), "rb");
                    bam_header_destroy(bam_header_read(fp));
                    bam = bam_init1();
                    char c;
                    while ( bam_read1(fp, bam)>0 ) {
                        char *qname = bam1_qname(bam);
                        c = qname[a0->k];
                        if ( bam_fps[c]==NULL ) {
                            name.push_back(c);
                            file = tmp_dir_s+name;
                            bam_fps[c] = bam_open(file.c_str(), "wb");
                            bam_header_write(bam_fps[c], bam_header_p);
                            name.erase(name.end()-1);
                            nnz++;
                        }
                        bam_write1(bam_fps[c], bam);
                        count[c]++;
                    }
                    bam_destroy1(bam);
                    bam_close(fp);
                    bam_close(bam_fps[0]);
                    if (nnz==0) {
                        boost::filesystem::remove(null_s);

                        a0->flag |= 0x01;
                        break;
                    } else if (nnz==1 && count[0]==0) {
                        bam_close(bam_fps[c]);
                        boost::filesystem::remove(null_s);
                        name.push_back(c);
                        file = tmp_dir_s+name;
                        boost::filesystem::remove(file);
                        a0->k++;
                    } else {
                        file = tmp_dir_s+a0->name;
                        boost::filesystem::remove(file);

                        a0->name = name;
                        file = tmp_dir_s + name;
                        boost::filesystem::rename(null_s, file);
                        int j=0;
                        int np = n;
                        for (int16_t c=1; c<K; c++) {
                            if (count[c]>0) {
                                BamTmpA *a = new BamTmpA;
                                a->k = a0->k+1;
                                a->i = j++;
                                a->flag = 0;
                                a->read_num = count[c];
                                a->parent = np;
                                a->name = a0->name;
                                a->name.push_back(c);
                                stack.push_back(a);
                                bam_close(bam_fps[c]);
                                n++;
                            }
                        }
                        break;
                    }
                }
            }
        }
    }
    string out_bam_s(out_bam_file);
    file = tmp_dir_s + stack[0]->name;
    boost::filesystem::rename(file, out_bam_s);
    delete stack[0];

    return 0;
}

int Bam::sort_by_name(bam1_t** bams, int n__bam, bam1_t** out_bams, int first_pos)
{
    for (int i=0; i<n__bam; i++) {
        out_bams[i] = bams[i];
    }
    return sort_by_name(out_bams, n__bam, first_pos);
}

struct BamTmpB{
    int k;
    int begin;
    int end;
};
/**
 *  @brief sort bam by names in memory, using dictionary sort
 *  @author zelin chen
 *  @date 2014-08-27
 *  @param  out_bams  input and output bams
 *  @param  n__bam    total number of bam records
 *  @param  first_pos first position in name for comparison
 *  @return  0 if success, currently always 0
 */
int Bam::sort_by_name(bam1_t** out_bams, int n__bam, int first_pos)
{
    if (n__bam<2) return 0;

    vector< BamTmpB > stack(1);
    stack[0].k    =first_pos;
    stack[0].begin=0;
    stack[0].end  =n__bam; 
    int n = 0;
    int j;
    const int K = 256;
    int begin[K+1], end[K+1];
    int nz_idx[K+1];
    while (!stack.empty()) {
        int k = stack[n].k;
        int b = stack[n].begin;
        int e = stack[n].end;
        k++;
        for (int i=0; i<K; i++) end[i]=0;
        for (int i=b; i<e; i++) {
            int32_t l_qname = out_bams[i]->core.l_qname;
            char *qname = bam1_qname(out_bams[i]);
            int16_t c = qname[k];
            end[c]++;
        }
        int m=0;
        for (int c=0; c<K; c++) {
            if (end[c]>0) nz_idx[m++] = c;
        }
        if (m==1) {
            if (nz_idx[0]==0) {
                stack.pop_back();
                n--;
            } else {
                stack[n].k++;
            }
        } else {
            stack.pop_back();
            n--;
            /// calculate the begin index in the vector for each alphabet
            // {{{
            begin[nz_idx[0]] =b;
            end[nz_idx[0]]+=b;
            for (j=1; j<m; j++) {
                begin[nz_idx[j]] = end[nz_idx[j-1]];
                end[nz_idx[j]] += end[nz_idx[j-1]];
            }
            for (j=0; j<m; j++) {
                BamTmpB tmpb;
                if (nz_idx[j]>0 && end[nz_idx[j]] > begin[nz_idx[j]]+1) {
                    tmpb.k    =k;
                    tmpb.begin=begin[nz_idx[j]];
                    tmpb.end  =end[nz_idx[j]];
                    stack.push_back(tmpb);
                    n++;
                }
            }
            // }}}
            j=0;
            for (int i=b; i<e; ) {
                while (i==end[nz_idx[j]]) {
                    j++;
                    if (j>=m) break;
                    i=begin[nz_idx[j]];
                }
                if (j>=m) break;
                int32_t l_qname = out_bams[i]->core.l_qname;
                char *qname = bam1_qname(out_bams[i]);
                int16_t c = qname[k];
                int i2=begin[c];
                if (i2==i) {
                    i++;
                } else {
                    bam1_t* bam = out_bams[i];
                    out_bams[i] = out_bams[i2];
                    out_bams[i2] = bam;
                }
                begin[c]++;
            }
        }
    }
    return 0;
}

int Bam::sort_by_name(vector<bam1_t*> & bams, vector<bam1_t*> & out_bams, int first_pos)
{
    out_bams = bams;
    return sort_by_name(bams, out_bams, first_pos);
}

int Bam::sort_by_name(vector<bam1_t*> & out_bams, int first_pos)
{
    return sort_by_name(out_bams.begin(), out_bams.end(), first_pos);
}

int Bam::sort_by_name(vector<bam1_t*>::iterator out_bams_begin, vector<bam1_t*>::iterator out_bams_end, int first_pos)
{
    int n__bam = out_bams_end - out_bams_begin;
    if (n__bam<2) return 0;

    vector< BamTmpB > stack(1);
    stack[0].k    =first_pos;
    stack[0].begin=0;
    stack[0].end  =n__bam; 
    int n = 0;
    int j;
    const int K = 256;
    int begin[K+1], end[K+1];
    int nz_idx[K+1];
    while (!stack.empty()) {
        int k = stack[n].k;
        int b = stack[n].begin;
        int e = stack[n].end;
        k++;
        for (int i=0; i<K; i++) end[i]=0;
        for (int i=b; i<e; i++) {
            char *qname = bam1_qname(*(out_bams_begin+i));
            int16_t c = qname[k];
            end[c]++;
        }
        int m=0;
        for (int c=0; c<K; c++) {
            if (end[c]>0) nz_idx[m++] = c;
        }
        if (m==1) {
            if (nz_idx[0]==0) {
                stack.pop_back();
                n--;
            } else {
                stack[n].k++;
            }
        } else {
            stack.pop_back();
            n--;
            /// calculate the begin index in the vector for each alphabet
            // {{{
            begin[nz_idx[0]] =b;
            end[nz_idx[0]]+=b;
            for (j=1; j<m; j++) {
                begin[nz_idx[j]] = end[nz_idx[j-1]];
                end[nz_idx[j]] += end[nz_idx[j-1]];
            }
            for (j=0; j<m; j++) {
                BamTmpB tmpb;
                if (nz_idx[j]>0 && end[nz_idx[j]] > begin[nz_idx[j]]+1) {
                    tmpb.k    =k;
                    tmpb.begin=begin[nz_idx[j]];
                    tmpb.end  =end[nz_idx[j]];
                    stack.push_back(tmpb);
                    n++;
                }
            }
            // }}}
            j=0;
            for (int i=b; i<e; ) {
                while (i==end[nz_idx[j]]) {
                    j++;
                    if (j>=m) break;
                    i=begin[nz_idx[j]];
                }
                if (j>=m) break;
                char *qname = bam1_qname((*(out_bams_begin+i)));
                int16_t c = qname[k];
                int i2=begin[c];
                if (i2==i) {
                    i++;
                } else {
                    bam1_t* bam = (*(out_bams_begin+i));
                    (*(out_bams_begin+i)) = (*(out_bams_begin+i2));
                    (*(out_bams_begin+i2)) = bam;
                }
                begin[c]++;
            }
        }
    }
    return 0;
}
// }}}

/**
  * converte bam to fastq (p1, p2, s)
  */
int bam_to_fastq(const string & bam_file, const string fastq_files[3])
// {{{
{
    bamFile bam_fp = bam_open(bam_file.c_str(), "rb");
    if (bam_fp == NULL) {
        cerr << "Fail to open BAM file " << bam_file << endl;
        return 1;
    }
    ofstream fastq_fs[3];
    string file;
    for (int i=0; i<3; i++) {
        fastq_fs[i].open(fastq_files[i].c_str());
        if (!fastq_fs[i].is_open()) {
            cerr << "Fail to open fastq file " << fastq_files[i] << endl;
            return 1;
        }
    }
    
    bam_header_t *bam_header_p = bam_header_read(bam_fp);
    vector<bam1_t*> bam1s[2], se_bam1s;
    bam1_t *bam1, *prev_bam1=NULL;
    bam1 = bam_init1();
	int r1=0, r=0;
    while ( r == 0 ) {
		r1 = bam_read1(bam_fp, bam1);
        if (prev_bam1==NULL) {
        } else if ( r1>0 && strcmp(bam1_qname(prev_bam1), bam1_qname(bam1))==0 ) {
        } else { // not the same qname
            int n0, n1;
            int n = se_bam1s.size();
            if (n>0) {
                int i1=0;
                int32_t l_qseq=0;
                for (int j=0; j<se_bam1s.size(); j++) {
					if (se_bam1s[j]->core.l_qseq > l_qseq) { i1=j; l_qseq = se_bam1s[i1]->core.l_qseq; }
                }
                char *qname = bam1_qname(se_bam1s[i1]);
                uint8_t *qual0 = bam1_qual(se_bam1s[i1]);
                char *qseq = new char[l_qseq+1];
                char* qual33 = new char[l_qseq+1];
                if ( bam1_strand(se_bam1s[i1]) ) { // reverse
                    int32_t j1=l_qseq;
                    qseq[j1]=0;
                    qual33[j1]=0;
                    for (int32_t j=0; j<l_qseq; j++) {
                        int8_t a = bam1_seqi(bam1_seq(se_bam1s[i1]), j);
                        switch(a) {
                            case 0x1: a=0x08; break;
                            case 0x2: a=0x04; break;
                            case 0x4: a=0x02; break;
                            case 0x8: a=0x01; break;
							default: break;
                        }
                        j1--;
                        qseq[j1] = bam_nt16_rev_table[a];
                        qual33[j1] = qual0[j]+33;
                    }
                } else {
                    int32_t j;
                    for (j=0; j<l_qseq; j++) {
                        int8_t a = bam1_seqi(bam1_seq(se_bam1s[i1]), j);
                        qseq[j] = bam_nt16_rev_table[a];
                        qual33[j] = qual0[j]+33;
                    }
                    qseq[j]=0;
                    qual33[j]=0;
                }
                fastq_fs[0] << "@" << qname << "\n";
                fastq_fs[0] << qseq << "\n" << "+" << "\n" << qual33 << "\n";
                delete qseq;
                delete qual33;
                BOOST_FOREACH (bam1_t* b, se_bam1s) {
                    bam_destroy1(b);
                }
                se_bam1s.clear();
            }

            int k[2];
            n0 = bam1s[0].size();
            n1 = bam1s[1].size();
			if ( n0>0 || n1>0 ) {
				if (n0>0 && n1==0) {
					k[0]=0;
				} else if (n1>0 && n0==0) {
					k[1]=0;
				} else if (n0>0 && n1>0) {
					k[0]=1; k[1]=2;
				}

				for (int i=0; i<2; i++) {
					if (bam1s[i].size()==0) continue;
					int i1=0;
					int32_t l_qseq=0;
					for (int j=0; j<bam1s[i].size(); j++) {
						if (bam1s[i][j]->core.l_qseq > l_qseq) { i1=j; l_qseq = bam1s[i][i1]->core.l_qseq; }
					}
					char *qname = bam1_qname(bam1s[i][i1]);
					uint8_t *qual0 = bam1_qual(bam1s[i][i1]);
					char *qseq = new char[l_qseq+1];
					char* qual33 = new char[l_qseq+1];
					if ( bam1_strand(bam1s[i][i1]) ) { // reverse
						int32_t j1=l_qseq;
						qseq[j1]=0;
						qual33[j1]=0;
						for (int32_t j=0; j<l_qseq; j++) {
							int8_t a = bam1_seqi(bam1_seq(bam1s[i][i1]), j);
							switch(a) {
								case 0x1: a=0x08; break;
								case 0x2: a=0x04; break;
								case 0x4: a=0x02; break;
								case 0x8: a=0x01; break;
								default: break;
							}
							j1--;
							qseq[j1] = bam_nt16_rev_table[a];
							qual33[j1] = qual0[j]+33;
						}
					} else {
						int32_t j;
						for (j=0; j<l_qseq; j++) {
							int8_t a = bam1_seqi(bam1_seq(bam1s[i][i1]), j);
							qseq[j] = bam_nt16_rev_table[a];
							qual33[j] = qual0[j]+33;
						}
						qseq[j]=0;
						qual33[j]=0;
					}
					fastq_fs[k[i]] << "@" << qname << "\n";
					fastq_fs[k[i]] << qseq << "\n" << "+" << "\n" << qual33 << "\n";
					delete qseq;
					delete qual33;
					BOOST_FOREACH (bam1_t* b, bam1s[i]) {
						bam_destroy1(b);
					}
					bam1s[i].clear();
				}
			}
        }
		if ( r1>0 ) {
			if ( (bam1->core.flag & BAM_FPAIRED) == 0) {
				se_bam1s.push_back(bam1);
			} else {
				if ( (bam1->core.flag & BAM_FREAD1) !=0 ) {
					bam1s[0].push_back(bam1);
				} else if ( (bam1->core.flag & BAM_FREAD2) !=0 ) {
					bam1s[1].push_back(bam1);
				} else {
				}
			}
			prev_bam1 = bam1;
			bam1 = bam_init1();
		} else {
			r = 1;
		}
    }
    bam_destroy1(bam1);

    bam_close(bam_fp);
    for (int i=0; i<3; i++) {
        fastq_fs[i].close();
    }
    return 0;
}
// }}}

/**
  * converte bam to fasta and qual (p1, p2, s)
  */
int bam_to_fasta_qual(const string & bam_file, const string fasta_files[3], const string qual_files[3])
// {{{
{
    bamFile bam_fp = bam_open(bam_file.c_str(), "rb");
    if (bam_fp == NULL) {
        cerr << "Fail to open BAM file " << bam_file << endl;
        return 1;
    }

    const uint32_t buf_l=16384;
    uint32_t fasta_buf_ul[3] = {0,0,0};
    uint32_t qual_buf_ul[3] = {0,0,0};
    char fasta_buf[3][buf_l+1];
    char qual_buf[3][buf_l+1];

    ofstream fasta_fs[3];
    ofstream qual_fs[3];
    /*
    for (int i=0; i<3; i++) {
        fasta_fs[i].rdbuf()->pubsetbuf(fasta_buf[i], buf_l);
        qual_fs[i].rdbuf()->pubsetbuf(qual_buf[i], buf_l);
    }
    */
    string file;
    for (int i=0; i<3; i++) {
        fasta_fs[i].open(fasta_files[i].c_str());
        if (!fasta_fs[i].is_open()) {
            cerr << "Fail to open fasta file " << fasta_files[i] << endl;
            return 1;
        }
        qual_fs[i].open(qual_files[i].c_str());
        if (!qual_fs[i].is_open()) {
            cerr << "Fail to open qual file " << qual_files[i] << endl;
            return 1;
        }
    }
    
    bam_header_t *bam_header_p = bam_header_read(bam_fp);
    vector<bam1_t*> bam1s[2], se_bam1s;
    bam1_t *bam1, *prev_bam1=NULL;
    bam1 = bam_init1();
	int r1=0, r=0;
    while ( r==0 ) {
		r1 = bam_read1(bam_fp, bam1);
        if (prev_bam1==NULL) {
        } else if ( r1>0 && strcmp(bam1_qname(prev_bam1), bam1_qname(bam1))==0 ) {
        } else { // not the same qname
            int n0, n1;
            int n = se_bam1s.size();
            if (n>0) {
                int i1=0;
                int l0=0, l;
                int32_t l_qseq=0;
                int32_t l_qname = se_bam1s[0]->core.l_qname;
                for (int j=0; j<se_bam1s.size(); j++) {
					if (se_bam1s[j]->core.l_qseq > l_qseq) { i1=j; l_qseq = se_bam1s[i1]->core.l_qseq; }
                }
                char *qname = bam1_qname(se_bam1s[i1]);
                uint8_t *qual0 = bam1_qual(se_bam1s[i1]);
                char *qseq = new char[l_qseq+1];
                if ( bam1_strand(se_bam1s[i1]) ) { // reverse
                    int32_t j1=l_qseq;
                    qseq[j1]=0;
                    for (int32_t j=0; j<l_qseq; j++) {
                        int8_t a = bam1_seqi(bam1_seq(se_bam1s[i1]), j);
                        switch(a) {
                            case 0x1: a=0x08; break;
                            case 0x2: a=0x04; break;
                            case 0x4: a=0x02; break;
                            case 0x8: a=0x01; break;
							default: break;
                        }
                        j1--;
                        qseq[j1] = bam_nt16_rev_table[a];
                    }
                } else {
                    int32_t j;
                    for (j=0; j<l_qseq; j++) {
                        int8_t a = bam1_seqi(bam1_seq(se_bam1s[i1]), j);
                        qseq[j] = bam_nt16_rev_table[a];
                    }
                    qseq[j]=0;
                }

                l0=l_qname+2 + l_qseq+1;
                l = fasta_buf_ul[0];
                if (l + l0 > buf_l) {
                    fasta_fs[0].write(fasta_buf[0], l*sizeof(char));
                    fasta_buf_ul[0] = 0;
                    l = 0;
                }
                fasta_buf[0][l++] = '>';
                memcpy(fasta_buf[0]+l, qname, l_qname);
                l+=l_qname;
                fasta_buf[0][l++] = '\n';
                memcpy(fasta_buf[0]+l, qseq, l_qseq);
                l+=l_qseq;
                fasta_buf[0][l++] = '\n';
                fasta_buf_ul[0] = l;

                l0 = l_qname+2 + l_qseq*4;
                l = qual_buf_ul[0];
                if (l + l0 > buf_l) {
                    qual_fs[0].write(qual_buf[0], l*sizeof(char));
                    qual_buf_ul[0] = 0;
                    l = 0;
                }
                qual_buf[0][l++] = '>';
                memcpy(qual_buf[0]+l, qname, l_qname);
                l+=l_qname;
                qual_buf[0][l++] = '\n';
                for (int32_t i=0; i<l_qseq; i++) {
                    if (i>0) qual_buf[0][l++]=' ';
                    int8_t a = qual0[i], ai=0;
                    char b[3];
                    while (a>=10) {
                        b[ai++] = (a%10)|0x30;
                        a/=10;
                    }
                    b[ai] = a|0x30;
                    while (ai>=0) { qual_buf[0][l++]=b[ai--]; }
                }
                qual_buf[0][l++] = '\n';
                qual_buf_ul[0] = l;

           //   fasta_fs[0] << ">" << qname << "\n" << qseq << "\n";
           //   qual_fs[0] << ">" << qname << "\n";
           //   for (int32_t i=0; i<l_qseq; i++) {
           //       if (i>0) qual_fs[0] << " ";
           //       qual_fs[0] << (int16_t)qual0[i];
           //   }
           //   qual_fs[0] << "\n";
                delete qseq;
                BOOST_FOREACH (bam1_t* b, se_bam1s) {
                    bam_destroy1(b);
                }
                se_bam1s.clear();
            }

            int k[2];
            n0 = bam1s[0].size();
            n1 = bam1s[1].size();
			if (n0>0 || n1>0) {
				if (n0>0 && n1==0) {
					k[0]=0;
				} else if (n1>0 && n0==0) {
					k[1]=0;
				} else if (n0>0 && n1>0) {
					k[0]=1; k[1]=2;
				}

				for (int i=0; i<2; i++) {
					if (bam1s[i].size()==0) continue;
					int l0=0, l=0;
					int i1=0;
					int32_t l_qseq=0;
					int32_t l_qname = bam1s[i][0]->core.l_qname;
					for (int j=0; j<bam1s[i].size(); j++) {
						if (bam1s[i][j]->core.l_qseq > l_qseq) { i1=j; l_qseq = bam1s[i][i1]->core.l_qseq; }
					}
					char *qname = bam1_qname(bam1s[i][i1]);
					uint8_t *qual0 = bam1_qual(bam1s[i][i1]);
					char *qseq = new char[l_qseq+1];
					if ( bam1_strand(bam1s[i][i1]) ) { // reverse
						int32_t j1=l_qseq;
						qseq[j1]=0;
						for (int32_t j=0; j<l_qseq; j++) {
							int8_t a = bam1_seqi(bam1_seq(bam1s[i][i1]), j);
							switch(a) {
								case 0x1: a=0x08; break;
								case 0x2: a=0x04; break;
								case 0x4: a=0x02; break;
								case 0x8: a=0x01; break;
								default: break;
							}
							j1--;
							qseq[j1] = bam_nt16_rev_table[a];
						}
					} else {
						int32_t j;
						for (j=0; j<l_qseq; j++) {
							int8_t a = bam1_seqi(bam1_seq(bam1s[i][i1]), j);
							qseq[j] = bam_nt16_rev_table[a];
						}
						qseq[j]=0;
					}

					int ki = k[i];
					l0= l_qname+2 + l_qseq+1;
					l = fasta_buf_ul[ki];
					if (l + l0 > buf_l) {
						fasta_fs[ki].write(fasta_buf[ki], l*sizeof(char));
						fasta_buf_ul[ki] = 0;
						l = 0;
					}
					fasta_buf[ki][l++] = '>';
					memcpy(fasta_buf[ki]+l, qname, l_qname);
					l+=l_qname;
					fasta_buf[ki][l++] = '\n';
					memcpy(fasta_buf[ki]+l, qseq, l_qseq);
					l+=l_qseq;
					fasta_buf[ki][l++] = '\n';
					fasta_buf_ul[ki] = l;

					l0 = l_qname+2 + l_qseq*4;
					l = qual_buf_ul[ki];
					if (l + l0 > buf_l) {
						qual_fs[ki].write(qual_buf[ki], l*sizeof(char));
						qual_buf_ul[ki] = 0;
						l = 0;
					}
					qual_buf[ki][l++] = '>';
					memcpy(qual_buf[ki]+l, qname, l_qname);
					l+=l_qname;
					qual_buf[ki][l++] = '\n';
					for (int32_t j=0; j<l_qseq; j++) {
						if (j>0) qual_buf[ki][l++]=' ';
						int8_t a = qual0[i], ai=0;
						char b[3];
						while (a>=10) {
							b[ai++] = (a%10)|0x30;
							a/=10;
						}
						b[ai] = a|0x30;
						while (ai>=0) { qual_buf[ki][l++]=b[ai--]; }
					//  string s = boost::lexical_cast<string>((int16_t)qual0[j]);
					//  memcpy(qual_buf[ki]+l, s.c_str(), s.size());
					//  l+=s.size();
					}
					qual_buf[ki][l++] = '\n';
					qual_buf_ul[ki] = l;

				//  fasta_fs[k[i]] << ">" << qname << "\n" << qseq << "\n";
				//  qual_fs[k[i]] << ">" << qname << "\n";
				//  for (int32_t j=0; j<l_qseq; j++) {
				//      if (j>0) qual_fs[k[i]] << " ";
				//      qual_fs[k[i]] << (int16_t)qual0[j];
				//  }
				//  qual_fs[k[i]] << "\n";
					delete qseq;
					BOOST_FOREACH (bam1_t* b, bam1s[i]) {
						bam_destroy1(b);
					}
					bam1s[i].clear();
				}
			}
        }
		if (r1>0) {
			if ( (bam1->core.flag & BAM_FPAIRED) == 0) {
				se_bam1s.push_back(bam1);
			} else {
				if ( (bam1->core.flag & BAM_FREAD1) !=0 ) {
					bam1s[0].push_back(bam1);
				} else if ( (bam1->core.flag & BAM_FREAD2) !=0 ) {
					bam1s[1].push_back(bam1);
				} else {
				}
			}
			prev_bam1 = bam1;
			bam1 = bam_init1();
		} else {
			r=1;
		}
    }
    bam_destroy1(bam1);

    bam_close(bam_fp);

    for (int i=0; i<3; i++) {
		fasta_fs[i].write(fasta_buf[i], fasta_buf_ul[i]*sizeof(char));
		qual_fs[i].write(qual_buf[i], qual_buf_ul[i]*sizeof(char));
        fasta_fs[i].close();
        qual_fs[i].close();
    }
    return 0;
}
// }}}

