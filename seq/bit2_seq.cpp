#include "bit2_seq.hpp"

namespace czl_bio {
// class Bit2Seq
// {{{

uint8_t Bit2Seq::nt1_to_bit2(char c)
{
    uint8_t b;
    switch(c) {
    case 'A':
    case 'a':
        b = 0x0;
        break;
    case 'G':
    case 'g':
        b = 0x1;
        break;
    case 'T':
    case 't':
        b = 0x2;
        break;
    case 'C':
    case 'c':
        b = 0x3;
        break;
    default:
        b = 0x4;
        break;
    }
    return b;
}

char Bit2Seq::bit2_to_nt1(uint8_t b)
{
    switch(b) {
    case 0x0:
        return 'A';
        break;
    case 0x1:
        return 'G';
        break;
    case 0x2:
        return 'T';
        break;
    case 0x3:
        return 'C';
        break;
    default:
        break;
    }
    return 'N';
}

uint8_t Bit2Seq::complement(uint8_t b)
{
    return b ^ 0x2;
}

int Bit2Seq::nt1_to_bit2(string const & nt, size_t *n, vector<uint8_t> & b)
{
    *n = nt.size();
    b.assign( (*n+3)>>2, 0 );
    for (size_t i=0; i < *n; i++) {
        b[i>>2] |= (nt1_to_bit2(nt[i])&0x3) << ((i&0x3)<<1);
    }
    return 0;
}

Bit2Seq::Bit2Seq(string const & nt)
{
    n_m = nt.size();
    seq_m.assign( (n_m+3)>>2, 0 );
    for (size_t i=0; i<n_m; i++) {
        seq_m[i>>2] |= (nt1_to_bit2(nt[i])&0x3) << ((i&0x3)<<1);
    }
}

Bit2Seq::Bit2Seq(Bit2Seq const & src)
{
    n_m   = src.n_m;
    seq_m = src.seq_m;
}

int Bit2Seq::copy_from_nt1(string const & nt)
{
    n_m = nt.size();
    seq_m.assign( (n_m+3)>>2, 0 );
    for (size_t i=0; i<n_m; i++) {
        seq_m[i>>2] |= (nt1_to_bit2(nt[i])&0x3) << ((i&0x3)<<1);
    }
	return 0;
}

int Bit2Seq::copy_to_nt1(string & nt)
{
    nt.resize(n_m);
    for (size_t i=0; i<n_m; i++) {
        nt[i] = bit2_to_nt1(get_nothrow(i));
    }
	return 0;
}

uint8_t Bit2Seq::get(size_t i) const
{
    if (i>=n_m) {
        out_of_range e(CZL_DBG_INFO);
        throw e;
    }
    return seq_m[i>>2] >> ((i&0x3)<<1) & 0x3; // seq_m[i/4]>>(i%4*2)
}

uint8_t Bit2Seq::get4(size_t i) const
{
    if (i>=n_m) {
        out_of_range e(CZL_DBG_INFO);
        throw e;
    }
    return seq_m[i>>2]; // seq_m[i/4]
}


uint8_t Bit2Seq::get_nothrow(size_t i) const
{
    return seq_m[i>>2] >> ((i&0x3)<<1) & 0x3; // seq_m[i/4]>>(i%4*2)
}

void Bit2Seq::set(size_t i, uint8_t b)
{
    if (i>=n_m) {
        out_of_range e(CZL_DBG_INFO);
        throw e;
    }
    seq_m[i>>2] &= ~(0x3 << ((i&0x3)<<1));
    seq_m[i>>2] |= (b&0x3) << ((i&0x3)<<1);
}

void Bit2Seq::set4(size_t i, uint8_t b)
{
    if (i>=n_m) {
        out_of_range e(CZL_DBG_INFO);
        throw e;
    }
    seq_m[i>>2] = b;
}

void Bit2Seq::set_nothrow(size_t i, uint8_t b)
{
    seq_m[i>>2] &= ~(0x3 << ((i&0x3)<<1));
    seq_m[i>>2] |= (b&0x3) << ((i&0x3)<<1);
}

int Bit2Seq::complement()
{
    uint8_t mask = 0xaa; // bit: 10101010
    size_t j=0;
    for (j=0; j < (n_m>>2); j++) {
        seq_m[j] ^= mask;
    }
    if ( (n_m&0x3) != 0 ) {
        seq_m[j] ^= mask;
    }
    return 0;
}

int Bit2Seq::reverse()
{
    size_t i1 = n_m--;
    for (size_t i=0; i < i1; i++) {
        uint8_t a = get_nothrow(i);
        set_nothrow(i, get_nothrow(i1));
        set_nothrow(i1, a);
        i1--;
    }
    return 0;
}

int Bit2Seq::rev_compl()
{
    size_t i = 0;
    size_t i1 = n_m--;
    for (i=0; i < i1; i++) {
        uint8_t a = get_nothrow(i) ^ 0x2;
        set_nothrow(i, get_nothrow(i1) ^ 0x2);
        set_nothrow(i1, a);
        i1--;
    }
    if (i==i1) {
        set_nothrow(i, get_nothrow(i) ^ 0x2);
    }
    return 0;
}

int Bit2Seq::subseq(Bit2Seq & b, size_t pos, size_t len)
{
	if ( len > n_m-pos ) len = n_m-pos;
    b.resize(len);
    size_t j=0;
    for (size_t i=pos, j=0; j<len && i<n_m ; i++, j++) {
        b.set_nothrow(j, get_nothrow(i));
    }
    return 0;
}

void Bit2Seq::resize(size_t n)
{
    if ( (n+3)>>2 != (n_m+3)>>2 ) seq_m.resize((n+3)>>2);
    n_m = n;
}

int Bit2Seq::apply_fun1(int (*pf)(void*))
{
    for ( int64_t i=0; i<n_m; i++ ) {
        uint8_t b = seq_m[i>>2]>>((i&0x3)<<1);
        pf(&b);
        seq_m[i>>2] &= ~(0x3<<((i&0x3)<<1));
        seq_m[i>>2] |= b<<((i&0x3)<<1);
    }
	return 0;
}

Bit2Seq::Iter Bit2Seq::begin()
{
    return Iter(this);
}

Bit2Seq::Iter Bit2Seq::end()
{
    return Iter(this)+n_m;
}


Bit2Seq & Bit2Seq::operator =(Bit2Seq const & src)
{
    n_m   = src.n_m;
    seq_m = src.seq_m;
    return *this;
}

Bit2Seq Bit2Seq::operator +(Bit2Seq const & b)
{
    Bit2Seq mb;
    mb.seq_m = seq_m;
    mb.n_m = n_m + b.n_m;
    mb.resize(mb.n_m);
    for (size_t i=0, j=n_m; i<b.n_m; i++, j++) {
        mb.set_nothrow(j, b.get_nothrow(i));
    }
    return mb;
}

Bit2Seq & Bit2Seq::operator +=(Bit2Seq const & b)
{
    size_t j = this->n_m;
    this->resize(this->n_m+b.n_m);
    for (size_t i=0; i<b.n_m; i++, j++) {
        this->set_nothrow(j, b.get_nothrow(i));
    }
    return *this;
}


bool Bit2Seq::operator ==(Bit2Seq const & src)
{
	if ( n_m == src.n_m ) {
		size_t j1 = n_m>>2;
		size_t j = 0;
		for (j=0; j<j1; j++) {
			if ( seq_m[j] != src.seq_m[j] ) break;
		}
		if ( j<j1 ) {
			return false;
		} else {
			uint8_t k = (n_m&0x3)<<1;
			if ( k==0 )  {
				return true;
			} else {
				uint8_t mask = 0xff>>(8-k);
				if ( (seq_m[j1]&mask) == (src.seq_m[j1]&mask) ) return true;
				else return false;
			}
		}
	} else {
		return false;
	}
}

bool Bit2Seq::operator !=(Bit2Seq const & src)
{
	return !( *this==src );
}

bool Bit2Seq::operator <(Bit2Seq const & src)
{
	for (size_t i=0; i<n_m && i<src.n_m; i++) {
		if ( get_nothrow(i) < src.get_nothrow(i) ) return true;
		else if ( get_nothrow(i) > src.get_nothrow(i) ) return false;
	}
	if ( n_m < src.n_m ) return true;
	else return false;
}

bool Bit2Seq::operator >(Bit2Seq const & src)
{
	for (size_t i=0; i<n_m && i<src.n_m; i++) {
		if ( get_nothrow(i) > src.get_nothrow(i) ) return true;
		else if ( get_nothrow(i) < src.get_nothrow(i) ) return false;
	}
	if ( n_m > src.n_m ) return true;
	else return false;
}

// }}}

// Bit2Seq::Iter
// {{{
Bit2Seq::Iter::Iter(Bit2Seq * b)
{
    bseq_m = b;
    i_m    = 0;
}

Bit2Seq::Iter::Iter(Iter const & src)
{
    bseq_m = src.bseq_m;
    i_m    = src.i_m;
}

Bit2Seq::Iter & Bit2Seq::Iter::operator =(Bit2Seq::Iter const & src)
{
    bseq_m = src.bseq_m;
    i_m    = src.i_m;
    return *this;
}

uint8_t Bit2Seq::Iter::operator *()
{
    if ( i_m>=bseq_m->n_m ) {
        string s("Current range ");
        s += itos(i_m) + ", (>=" + itos(bseq_m->size()) + ")";
        throw out_of_range(s);
    }
    return bseq_m->get_nothrow(i_m);
}

Bit2Seq::Iter Bit2Seq::Iter::operator +(int64_t n)
{
    Iter iter(*this);
    iter.i_m += n;
    if ( iter.i_m<0 ) iter.i_m = 0;
    else if ( iter.i_m>=iter.bseq_m->size() ) iter.i_m = iter.bseq_m->size();
    return iter;
}

Bit2Seq::Iter Bit2Seq::Iter::operator -(int64_t n)
{
    Iter iter(*this);
    iter.i_m -= n;
    if ( iter.i_m<0 ) iter.i_m = 0;
    else if ( iter.i_m>=iter.bseq_m->size() ) iter.i_m = iter.bseq_m->size();
    return iter;
}

int64_t Bit2Seq::Iter::operator -(Iter & iter)
{
    if ( bseq_m == iter.bseq_m ) return i_m - iter.i_m;
    else {
        runtime_error e("Iter.seq_m is different.");
        throw e;
    }
}

Bit2Seq::Iter & Bit2Seq::Iter::operator ++()
{
    if (i_m < bseq_m->size()) i_m++;
    return (*this);
}

Bit2Seq::Iter Bit2Seq::Iter::operator ++(int)
{
    Bit2Seq::Iter iter(*this);
    if (i_m < iter.bseq_m->size()) i_m++;
    return iter;
}

Bit2Seq::Iter & Bit2Seq::Iter::operator --()
{
    if (i_m>0) i_m--;
    return *this;
}

Bit2Seq::Iter Bit2Seq::Iter::operator --(int n)
{
    Bit2Seq::Iter iter(*this);
    if (i_m>0) i_m--;
    return iter;
}

Bit2Seq::Iter & Bit2Seq::Iter::operator +=(int64_t n)
{
    i_m += n;
    if ( i_m<0 ) i_m = 0;
    else if ( i_m>=bseq_m->n_m ) i_m = bseq_m->n_m;
    return (*this);
}

Bit2Seq::Iter & Bit2Seq::Iter::operator -=(int64_t n)
{
    i_m -= n;
    if ( i_m<0 ) i_m = 0;
    else if ( i_m>=bseq_m->n_m ) i_m = bseq_m->n_m;
    return (*this);
}

bool Bit2Seq::Iter::operator ==(Iter const & b)
{
    return ( &bseq_m == &b.bseq_m );
}

bool Bit2Seq::Iter::operator !=(Iter const & b)
{
    return !( &bseq_m == &b.bseq_m );
}

bool Bit2Seq::Iter::operator >(Iter const & b)
{
    if (bseq_m == b.bseq_m) {
        return i_m > b.i_m;
    } else {
        runtime_error e("Iter.seq_m or Iter.n_m is different.");
        throw e;
    }
    return false;
}

bool Bit2Seq::Iter::operator <(Bit2Seq::Iter const & b)
{
    if (bseq_m == b.bseq_m) {
        return i_m < b.i_m;
    } else {
        runtime_error e("Iter.seq_m is different.");
        throw e;
    }
    return false;
}

bool Bit2Seq::Iter::operator >=(Bit2Seq::Iter const & b)
{
    if (bseq_m == b.bseq_m) {
        return i_m >= b.i_m;
    } else {
        runtime_error e("Iter.seq_m is different.");
        throw e;
    }
    return false;
}

bool Bit2Seq::Iter::operator <=(Bit2Seq::Iter const & b)
{
    if (bseq_m == b.bseq_m) {
        return i_m <= b.i_m;
    } else {
        runtime_error e("Iter.seq_m is different.");
        throw e;
    }
    return false;
}
// }}}

};
