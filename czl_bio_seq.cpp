#include "czl_bio_seq.hpp"
#include <numeric>

namespace czl_bio {
    
int detect_base_qual_type(vector<string> const & qual_v)
{
    int min=100000, max=0;
    int r;
    int i, j;
    for (i=0; i<qual_v.size(); i++) {
        string const & qual = qual_v[i];
        for (j=0; j<qual.size(); j++) {
            if ( qual[j]<min ) min = qual[j];
            if ( qual[j]>max ) max = qual[j];
        }
    }
    if ( min < 33 ) {
        r = 0;
    } else if ( min<64 ) {
        r = 33;
    } else {
        r = 64;
    }
    return r;
}

int trim_seq_by_qual(string const & qual,
           int win_size, int step_size, 
           int qual_thres, float win_frac_thres, int win_min,
           int cut_pos[2])
{
    return trim_seq_by_qual(qual, 0, qual.size(), win_size, step_size,
               qual_thres, win_frac_thres, win_min, cut_pos);
}

int trim_seq_by_qual(string const & qual, int begin, int length,
           int win_size, int step_size, 
           int qual_thres, float win_frac_thres, int win_min,
           int cut_pos[2])
// {{{
{
    int i, j, k;
    if ( win_min > win_size ) win_min = win_size;
    if ( step_size > win_size ) step_size = win_size;
    int wn=0;
    int end = begin+length;
    if ( begin >= end ) {
        cut_pos[0] = cut_pos[1] = 0;
        return 0;
    }
    if ( end>qual.size() ) end = qual.size();

    int b = begin, e = begin+win_size;
    if ( e > end ) e = end;
    for ( i=b; i < e; i++ ) {
        if ( qual[i] >= qual_thres ) {
            wn++;
        }
    }
    // if sequense size < window size
    if ( length <= win_size ) {
        if ( wn >= win_frac_thres*qual.size() && wn >= win_min) {
            int i0 = begin;
            while ( qual[i0]<qual_thres ) i0++;
            int i1 = end-1;
            while ( qual[i1]<qual_thres ) i1--;
            cut_pos[0] = i0;
            cut_pos[1] = i1;
        } else {
            cut_pos[0] = cut_pos[1] = 0;
        }
        return 0;
    }
    //

    // trim left
    while ( e<end ) {
        if ( wn >= win_frac_thres*(e-b) && wn >= win_min) {
            int i0 = b;
            while ( qual[i0]<qual_thres ) i0++;
            if ( i0==b ) {
                while ( i0>=0 && qual[i0]>=qual_thres ) i0--;
                i0++;
            }
            cut_pos[0] = i0;
            break;
        } else {
            for ( j=0; j<step_size && e<end; j++ ) {
                if ( qual[b] >= qual_thres ) {
                    wn--;
                }
                b++;
                if ( qual[e] >= qual_thres ) {
                    wn++;
                }
                e++;
            }
        }
    }

    if ( cut_pos[0] >=begin ) {
        // trim right
        b = end-win_size;
        e = end;
        wn = 0;
        for ( i=b; i<e; i++ ) {
            if ( qual[i] >= qual_thres ) {
                wn++;
            }
        }
        while ( b>=cut_pos[0] ) {
            if ( wn >= win_frac_thres*(e-b) && wn >= win_min) {
                int i1 = e;
                if ( qual[i1-1] < qual_thres ) {
                    while ( qual[i1-1]<qual_thres ) i1--;
                } else {
                    while ( i1<end && qual[i1]>=qual_thres ) i1++;
                }
                cut_pos[1] = i1;
                break;
            } else {
                for ( j=0; j<step_size && b>=begin; j++ ) {
                    if ( qual[e] >= qual_thres ) {
                        wn--;
                    }
                    e--;
                    if ( qual[b] >= qual_thres ) {
                        wn++;
                    }
                    b--;
                }
            }
        }
        //
    }
    if ( cut_pos[0]>=cut_pos[1] ) {
        cut_pos[0] = cut_pos[1] = 0;
    }

    return 0;
}
// }}}

typedef struct TmpA {
    int b;
    int e;
    int gn;
//    int gap_bn; // length of bad gap next to this region
} TmpA;

/// this method is not for pair-end
/// merge nearby good block to get trimmed sequence
int trim_seq_by_qual_m2(string const & qual, int begin, int length,
            int qual_thres, float frac_thres, int size_thres,
            int cut_pos[2])
// {{{
{
    vector<int> pos_v;
    int r;
    int i, j, k;
    int b, e;
    int n=0;

    i = begin;
    int end = begin+length;
    while ( i<end && qual[i] < qual_thres ) i++;
    if ( i==end ) {
        cut_pos[0] = 0;
        cut_pos[1] = 0;
        return 0;
    }

    int m = 1;
    for ( ; i<end; i++ ) {
        if ( qual[i]<qual_thres ) {
            if ( qual[i-1]<qual_thres ) {
            } else {
                pos_v.push_back(i);
                m = 0;
            }
        } else {
            if ( qual[i-1]<qual_thres ) {
                pos_v.push_back(i);
                m = 0;
            } else {
            }
        }
        m++;
    }
    pos_v.push_back(i);

    list<TmpA*> v;
    list<int> bv;
    for ( j=0; j<pos_v.size()-1; j+=2 ) {
        TmpA *p = new TmpA;
        p->b = pos_v[j];
        p->e = pos_v[j+1];
        p->gn = p->e - p->b;
//        if ( j+2>=pos_v.size() ) {
//            p->gap_bn = 0;
//        } else {
//            p->gap_bn = pos_v[j+2]-pos_v[j+1];
//        }
        v.push_back(p);
    }
    list<TmpA*>::iterator v_it = v.begin(), v_it_end, v_it1, max_v_it;
    float xf;
    while ( v.size()>2 ) {
        xf=0;
        v_it = v.begin();
        v_it_end = --v.end();
        while ( v_it!=v_it_end ) {
            v_it1 = v_it;
            v_it1++;
            int gn = (*v_it)->gn + (*v_it1)->gn;
            int l = (*v_it1)->e - (*v_it)->b;
            if ( (float)gn/l > xf ) {
                max_v_it = v_it;
                xf = (float)gn/l;
            }
            v_it++;
        }
        if ( xf >= frac_thres ) {
            v_it1 = max_v_it;
            v_it1++;
            (*max_v_it)->e = (*v_it1)->e;
            (*max_v_it)->gn += (*v_it1)->gn;
//            (*max_v_it)->gap_bn = (*v_it1)->gap_bn;
            v.erase(v_it1);
        } else {
            break;
        }
    }
    xf = 0;
    v_it = v.begin();
    while ( v_it!=v.end()) {
        if ( (*v_it)->gn>0 ) {
            int l = (*v_it)->e - (*v_it)->b;
            float f = (float)(*v_it)->gn / l;
            if ( l>=size_thres && f>=frac_thres && f>xf ) {
                xf = f;
                max_v_it = v_it;
            }
        }
        v_it++;
    }
    if ( xf>=frac_thres ) {
        cut_pos[0] = (*max_v_it)->b;
        cut_pos[1] = (*max_v_it)->e;
    } else {
        cut_pos[0] = 0;
        cut_pos[1] = 0;
    }
    for ( v_it = v.begin(); v_it!=v.end(); v_it++ ) {
        delete (*v_it);
    }

    return 0;
}
// }}}

int trim_seq_by_qual_m3(string const & qual,
           int win_size, int step_size, 
           int qual_thres, float win_frac_thres,
           int cut_pos[2])
{
    return trim_seq_by_qual_m3(qual, 0, qual.size(), win_size, step_size,
            qual_thres, win_frac_thres, cut_pos);
}

int trim_seq_by_qual_m3(string const & qual, int begin, int length,
           int win_size, int step_size, 
           int qual_thres, float win_frac_thres,
           int cut_pos[2])
// {{{
{
    int i, j, k;
    if ( step_size > win_size ) step_size = win_size;
    int wn=0;
    int end = begin+length;
    assert( end>begin );
    if ( end>qual.size() ) end = qual.size();
    for ( i=begin; i<end && i<win_size; i++ ) {
        if ( qual[i] >= qual_thres ) {
            wn++;
        }
    }
    // if sequense size < window size
    if ( length <= win_size ) {
        if ( wn >= win_frac_thres*qual.size() ) {
            int i0 = begin;
            while ( qual[i0]<qual_thres ) i0++;
            int i1 = end-1;
            while ( qual[i1]<qual_thres ) i1--;
            cut_pos[0] = i0;
            cut_pos[1] = i1;
        } else {
            cut_pos[0] = cut_pos[1] = 0;
        }
        return 0;
    }
    //

    // get good fraction for each windows
    int b = begin, e = begin+win_size;
    vector<int> pos_v;
    float prev_f = (float)wn/win_size;
    while ( e<end && prev_f < win_frac_thres ) {
        for ( j=0; j<step_size && e<end; j++ ) {
            if ( qual[b] >= qual_thres ) {
                wn--;
            }
            b++;
            if ( qual[e] >= qual_thres ) {
                wn++;
            }
            e++;
        }
        prev_f = (float)wn/(e-b);
    }
    pos_v.push_back(b);
    int max_i, max_l=0;
    while ( e<end ) {
        // calculate number of good bp in current window
        for ( j=0; j<step_size && e<end; j++ ) {
            if ( qual[b] >= qual_thres ) {
                wn--;
            }
            b++;
            if ( qual[e] >= qual_thres ) {
                wn++;
            }
            e++;
        }
        float f = (float)wn/(e-b);
        if ( (prev_f < win_frac_thres && f >= win_frac_thres) \
                || (prev_f >= win_frac_thres && f < win_frac_thres) ) {
            if ( prev_f >= win_frac_thres && b-pos_v.back()>max_l ) {
                max_l = b-pos_v.back();
                max_i = pos_v.size()-1;
            }
            pos_v.push_back(b);
            prev_f = f;
        }
    }
    pos_v.push_back(end);

    if ( max_l==0 ) {
        cut_pos[0] = 0;
        cut_pos[1] = 0;
    } else {
        b = pos_v[max_i];
        e = pos_v[max_i+1]+win_size;
        if ( e>end ) e = end;
        if ( qual[b] < qual_thres ) {
            while ( qual[b] < qual_thres ) b++;
        } else {
            while ( qual[b] >= qual_thres ) b--;
        }
        if ( qual[e-1] < qual_thres ) {
            while ( qual[e-1] < qual_thres ) e--;
        } else {
            while ( qual[e-1] >= qual_thres ) e++;
        }
        cut_pos[0] = b;
        cut_pos[1] = e;
    }

    return 0;
}
// }}}

int trim_seq_by_qual_m4(string const & qual, 
           int qual_thres, float frac_thres,
           int cut_pos[4])
{
    return trim_seq_by_qual_m4(qual, 0, qual.size(), qual_thres, frac_thres, cut_pos);
}

int trim_seq_by_qual_m4(string const & qual, int begin, int length,
           int qual_thres, float frac_thres,
           int cut_pos[4])
// {{{
{
    vector<int> pos_v;
    int r;
    int i, j, k;
    int b, e;
    int n=0;

    i = begin;
    int end = begin+length;
    while ( i<end && qual[i] < qual_thres ) i++;
    if ( i==end ) {
        cut_pos[0] = 0;
        cut_pos[1] = 0;
        cut_pos[2] = 0;
        cut_pos[3] = 0;
        return 0;
    }

    int m = 1;
    pos_v.push_back(i);
    i++;
    for ( ; i<end; i++ ) {
        if ( qual[i]<qual_thres ) {
            if ( qual[i-1]<qual_thres ) {
            } else {
                pos_v.push_back(i);
                m = 0;
            }
        } else {
            if ( qual[i-1]<qual_thres ) {
                pos_v.push_back(i);
                m = 0;
            } else {
            }
        }
        m++;
    }
    pos_v.push_back(i);

    int max_b=0, max_e=0, max_l = 0;
    int max_b1=0, max_e1=0, max_h1 = 0;
    float max_f = 0.0, max_f1 = 0.0;
    for ( i=0; i<pos_v.size()-1; i+=2 ) {
        int h=0;
        for ( j=i; j<pos_v.size()-1; j+=2 ) {
            int l = pos_v[j+1]-pos_v[i];
            h += pos_v[j+1] - pos_v[j];
            float f = (float)h/l;
            if ( f >= frac_thres ) {
                if (l > max_l || l==max_l && f>max_f) {
                    max_b = pos_v[i];
                    max_e = pos_v[j+1];
                    max_l = l;
                    max_f = f;
                }
                if (h > max_h1 || h==max_h1 && f>max_f) {
                    max_b1 = pos_v[i];
                    max_e1 = pos_v[j+1];
                    max_h1 = h;
                    max_f1 = f;
                }
            }
        }
    }
    cut_pos[0] = max_b;
    cut_pos[1] = max_e;
    cut_pos[2] = max_b1;
    cut_pos[3] = max_e1;

    return 0;
}
// }}}

typedef struct TmpB {
    int b;
    int e;
    int h;
    float f;
} TmpB;
bool tmpb_great(TmpB const *a, TmpB const *b)
{
    return ( a->h > b->h || (a->h==b->h && a->f > b->f ) );
}

int trim_seq_by_qual_m4(string const & qual1, string const & qual2,
           int qual_thres,
           float frac_thres, float frac_thres1, float frac_thres2,
           int cut_pos[8])
{
    return trim_seq_by_qual_m4(qual1, 0, qual1.size(), qual2, 0, qual2.size(),
            qual_thres, frac_thres, frac_thres1, frac_thres2, cut_pos);
}

int trim_seq_by_qual_m4(string const & qual1, int begin1, int length1,
           string const & qual2, int begin2, int length2,
           int qual_thres,
           float frac_thres, float frac_thres1, float frac_thres2,
           int cut_pos[8])
// {{{
{
    vector<TmpB*> v1, v2;
    vector<int> pos_v;
    int r;
    int i, j, k, h, l;
    int b, e;
    int n=0;
    int end;
    float max_f, f;
    int max_h, max_l;
    int qual_thres1 = qual_thres;
    int qual_thres2 = qual_thres;

    // qual1
    // {{{
    i = begin1;
    end = begin1 + length1;
    while ( i<end && qual1[i] < qual_thres1 ) i++;
    if ( i==end ) {
        cut_pos[0] = 0;
        cut_pos[1] = 0;
        cut_pos[2] = 0;
        cut_pos[3] = 0;
    } else {
        pos_v.push_back(i);
        i++;
        for ( ; i<end; i++ ) {
            if ( qual1[i]<qual_thres1 ) {
                if ( qual1[i-1]<qual_thres1 ) {
                } else {
                    pos_v.push_back(i);
                }
            } else {
                if ( qual1[i-1]<qual_thres1 ) {
                    pos_v.push_back(i);
                } else {
                }
            }
        }
        pos_v.push_back(i);

        for ( i=0; i<pos_v.size()-1; i+=2 ) {
            int h = 0;
            for ( j=i; j<pos_v.size()-1; j+=2 ) {
                int l = pos_v[j+1]-pos_v[i];
                h += pos_v[j+1]-pos_v[j];
                float f = (float)h/l;
                if ( f >= frac_thres1 ) {
                    TmpB *p = new TmpB;
                    p->b = pos_v[i];
                    p->e = pos_v[j+1];
                    p->h = h;
                    p->f = f;
                    v1.push_back(p);
                }
            }
        }
    }
	if ( v1.size()>1 ) {
		sort(v1.begin(), v1.end(), tmpb_great);
		j = 0;
		max_f = v1[0]->f;
		for (i=1; i<v1.size(); i++) {
			if ( v1[i]->f > max_f ) {
				j++;
				max_f = v1[i]->f;
				if ( j!=i ) {
					delete v1[j];
					v1[j] = v1[i];
					v1[i] = NULL;
				}
			}
		}
		v1.resize(j+1);
	}
    // }}}

    // qual2
    // {{{
    pos_v.clear();
    i = begin2;
    end = begin2 + length2;
    while ( i<end && qual2[i] < qual_thres2 ) i++;
    if ( i==end ) {
        cut_pos[4] = 0;
        cut_pos[5] = 0;
        cut_pos[6] = 0;
        cut_pos[7] = 0;
    } else {
        pos_v.push_back(i);
        i++;
        for ( ; i<end; i++ ) {
            if ( qual2[i]<qual_thres2 ) {
                if ( qual2[i-1]<qual_thres2 ) {
                } else {
                    pos_v.push_back(i);
                }
            } else {
                if ( qual2[i-1]<qual_thres2 ) {
                    pos_v.push_back(i);
                } else {
                }
            }
        }
        pos_v.push_back(i);

        for ( i=0; i<pos_v.size()-1; i+=2 ) {
            int h = 0;
            for ( j=i; j<pos_v.size()-1; j+=2 ) {
                int l = pos_v[j+1]-pos_v[i];
                h += pos_v[j+1]-pos_v[j];
                float f = (float)h/l;
                if ( f >= frac_thres2 ) {
                    TmpB *p = new TmpB;
                    p->b = pos_v[i];
                    p->e = pos_v[j+1];
                    p->h = h;
                    p->f = f;
                    v2.push_back(p);
                }
            }
        }
    }
	if ( v2.size()>1 ) {
		sort(v2.begin(), v2.end(), tmpb_great);
		j = 0;
		max_f = v2[0]->f;
		for (i=1; i<v2.size(); i++) {
			if ( v2[i]->f > max_f ) {
				j++;
				max_f = v2[i]->f;
				if ( j!=i ) {
					delete v2[j];
					v2[j] = v2[i];
					v2[i] = NULL;
				}
			}
		}
		v2.resize(j+1);
	}
    // }}}

    // qual1 + qual2
    // {{{
    max_h = 0;
    max_l = 0;
    for ( i=0; i<8; i++ ) cut_pos[i] = 0;
	if ( v1.size()==0 ) {
		if ( v2.size() > 0 ) {
			for ( j=0; j<v2.size(); j++ ) {
				int l = v2[j]->e - v2[j]->b;
				h = v2[j]->h;
				f = (float)h / l;
				if ( f>frac_thres ) {
					if ( h > max_h ) {
						cut_pos[6] = v2[j]->b;
						cut_pos[7] = v2[j]->e;
						max_h = h;
					} 
					if ( l > max_l ) {
						cut_pos[4] = v2[j]->b;
						cut_pos[5] = v2[j]->e;
						max_l = l;
					}
				}
			}
		}
	} else {
		if ( v2.size() == 0 ) {
			for ( j=0; j<v1.size(); j++ ) {
				int l = v1[j]->e - v1[j]->b;
				h = v1[j]->h;
				f = (float)h / l;
				if ( f>frac_thres ) {
					if ( h > max_h ) {
						cut_pos[2] = v1[j]->b;
						cut_pos[3] = v1[j]->e;
						max_h = h;
					} 
					if ( l > max_l ) {
						cut_pos[0] = v1[j]->b;
						cut_pos[1] = v1[j]->e;
						max_l = l;
					}
				}
			}
		} else {
			for ( i=0; i<v1.size(); i++ ) {
				int l1 = v1[i]->e - v1[i]->b;
				for ( j=0; j<v2.size(); j++ ) {
					int l2 = v2[j]->e - v2[j]->b;
					h = v1[i]->h + v2[j]->h;
					f = (float)h / (l1+l2);
					if ( f>frac_thres ) {
						if ( h > max_h ) {
							cut_pos[2] = v1[i]->b;
							cut_pos[3] = v1[i]->e;
							cut_pos[6] = v2[j]->b;
							cut_pos[7] = v2[j]->e;
							max_h = h;
						} 
						if ( l1+l2 > max_l ) {
							cut_pos[0] = v1[i]->b;
							cut_pos[1] = v1[i]->e;
							cut_pos[4] = v2[j]->b;
							cut_pos[5] = v2[j]->e;
							max_l = l;
						}
					}
				}
			}
		}
	}
    // }}}

    for ( i=0; i<v1.size(); i++ ) delete v1[i];
    for ( i=0; i<v2.size(); i++ ) delete v2[i];

    return 0;
}
// }}}

short is_seq_low_qual(string const & qual,
           int qual_thres, float frac_thres, int min_len)
{
    int i;
    int wn;
    if ( qual.size() != qual.size() ) {
        string s="seq and qual size is not the same.";
        s += CZL_DBG_INFO;
        msg.error(s, 1);
        return -1;
    }

    for ( i=0; i<qual.size(); i++ ) {
        if ( qual[i] >= qual_thres ) {
            wn++;
        }
    }

    if ( wn>=frac_thres && wn>=min_len ) {
        return 0;
    } else {
        return 1;
    }

    return 0;
}

int get_high_qual_len(string const & qual,
           int begin, int len, int qual_thres)
{
    int i;
    int wn = 0;
    if ( qual.size() != qual.size() ) {
        string s="seq and qual size is not the same.";
        s += CZL_DBG_INFO;
        msg.error(s, 1);
        return -1;
    }
    for ( i=begin; i<begin+len; i++ ) {
        if ( qual[i] >= qual_thres ) {
            wn++;
        }
    }

    return wn;
}

int get_high_qual_len(string const & qual,
           int qual_thres)
{
    return get_high_qual_len( qual, 0, qual.size(), qual_thres);
}

float get_qual_mean(string const & qual, int begin, int len)
{
    assert(begin<qual.size());
    if ( len==0 || begin>qual.size() ) return 0; 
    float sum;
    int end = begin+len;
    assert(begin < end);
    if ( end>qual.size() ) {
        end = qual.size();
        len = end - begin;
    }
    for ( int i=begin; i<end; i++ ) {
        sum += qual[i];
    }
    sum /= len;
    return sum;
}

};
