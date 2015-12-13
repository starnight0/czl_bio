#include "bit2_seq.hpp"

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

int Bit2Seq::complement(uint8_t *b)
{
    (*b) ^= 0x2;
}

void nt_to_bit2(string const & s, short & len, vector<czl_bio::uint8_t> & b)
{
    len = s.size();
    short B=2;
    short sz = 8;
    short Bn = 4;
    short n = s.size();
    short bn = (n+Bn-1)/Bn;
    short i, j, k;
    b.resize(bn);
    j = 0;
    k = 0;
    for (i=0; i<n; i++) {
        czl_bio::uint8_t a = nt_to_bit2(s[i]);
        b[j] |= (a << k);
        k += B;
        if (k>=sz) {
            j++;
            k=0;
        }
    }
}

void bit2_rev_compl(short n, vector<czl_bio::uint8_t> & b)
{
    int i, j=0;
    short j0=0, j1=(n-1)/4;
    short k0=0, k1=((n-1)%4)*2;
    for (i=0; i<n/2; i++) {
        uint8_t t0 = ( (b[j0]>>k0)&0x3 )^0x2; // ^0x2 is complement
        uint8_t t1 = ( (b[j1]>>k1)&0x3 )^0x2;
        b[j0] &= ~(0x3<<k0);
        b[j1] &= ~(0x3<<k1);
        b[j0] |= t1 << k0;
        b[j1] |= t0 << k1;
        k0+=2;
        if (k0==8) {j0++; k0=0;}
        if (k1==0) {j1--; k1=6;}
        else k1-=2;
    }
}

czl_bio::uint8_t bit2_get_seqi(vector<czl_bio::uint8_t> const & seq, short i)
{
    int N=4;
    if (i/N >= seq.size()) return 0x4;
    return (seq[i/N] >> ((i%N)<<2)) & 0x3;
}

int bit2_set_seqi(vector<czl_bio::uint8_t> & seq, short i, czl_bio::uint8_t a)
{
    int N=4;
    if (i/N >= seq.size()) return -1;
    a &= 0x3;
    seq[i/N] |= a << ((i%N)<<2);
    return 0;
}

int apply_fun1(size_t n, vector<uint8_t> & seq, int (*pf)(void*))
{
}

// Bit2Seq::Iter
// {{{
Bit2Seq::Iter::Iter(size_t n, vector<uint8_t> & seq)
{
    n_m   = n;
    seq_m = seq;
    i_m   = 0;
}

Bit2Seq::Iter::Iter(Iter const & src)
{
    n_m  = src.n_m;
    seq_m = src.seq_m;
    i_m = src.i_m;
}

Iter & Bit2Seq::Iter::operator =(Iter const & src)
{
    n_m  = src.n_m;
    seq_m = src.seq_m;
    i_m = src.i_m;
    return *this;
}

uint8_t Bit2Seq::Iter::operator *()
{
    if ( i_m>=n_m ) {
        string s("Current range ");
        s += itos(i_m) + ", (>=" + itos(n_m) + ")";
        throw out_of_range(s);
    }
    return seq_m[i_m>>2]>>((i_m&0x3)<<1) &0x3;
}

Iter Bit2Seq::Iter::operator +(int64_t n)
{
    Iter iter(n_m);
    iter.i_m += n;
    if ( iter.i_m<0 ) iter.i_m = 0;
    else if ( iter.i_m>=iter.n_m ) iter.i_m = iter.n_m;
    return iter;
}

Iter Bit2Seq::Iter::operator -(int64_t n)
{
    Iter iter(n_m);
    iter.i_m -= n;
    if ( iter.i_m<0 ) iter.i_m = 0;
    else if ( iter.i_m>=iter.n_m ) iter.i_m = iter.n_m;
    return iter;
}

int64_t Bit2Seq::Iter::operator -(Iter & iter)
{
    if ( &seq_m == &iter.seq_m ) return i_m - iter.i_m;
    else {
        runtime_error e("Iter.seq_m is different.");
        throw e;
    }
}

Iter Bit2Seq::Iter::operator ++()
{
    Iter iter(*this);
    if (i_m<n_m) i_m++;
    return iter;
}

Iter & Bit2Seq::Iter::operator ++(Iter & iter)
{
    if (i_m<n_m) i_m++;
    return (*this);
}

Iter Bit2Seq::Iter::operator --()
{
    Iter iter(*this);
    if (i_m>0) i_m--;
    return iter;
}

Iter & Bit2Seq::Iter::operator --(Itet & iter)
{
    if (i_m>0) i_m--;
    return *this;
}

Iter & Bit2Seq::Iter::operator +=(int n)
{
    i_m += n;
    if ( i_m<0 ) i_m = 0;
    else if ( i_m>=n_m ) i_m = n_m;
    return (*this);
}

Iter & Bit2Seq::Iter::operator -=(int n)
{
    i_m -= n;
    if ( i_m<0 ) i_m = 0;
    else if ( i_m>=n_m ) i_m = n_m;
    return (*this);
}

bool Bit2Seq::Iter::operator ==(Iter const & b)
{
    return ( &seq_m == &b.seq_m && n_m==b.n_m && i_m==b.i_m );
}

bool Bit2Seq::Iter::operator !=(Iter const & b)
{
    return !( &seq_m == &b.seq_m && n_m==b.n_m && i_m==b.i_m );
}

bool Bit2Seq::Iter::operator >(Iter const & b)
{
    if (&seq_m == &b.seq_m && n_m==b.n_m) {
        return i_m > b.i_m;
    } else {
        runtime_error e("Iter.seq_m or Iter.n_m is different.");
        throw e;
    }
    return false;
}

bool Bit2Seq::Iter::operator <(Iter const & b)
{
    if (&seq_m == &b.seq_m && n_m==b.n_m) {
        return i_m < b.i_m;
    } else {
        runtime_error e("Iter.seq_m or Iter.n_m is different.");
        throw e;
    }
    return false;
}

bool Bit2Seq::Iter::operator >=(Iter const & b)
{
    if (&seq_m == &b.seq_m && n_m==b.n_m) {
        return i_m >= b.i_m;
    } else {
        runtime_error e("Iter.seq_m or Iter.n_m is different.");
        throw e;
    }
    return false;
}

bool Bit2Seq::Iter::operator <=(Iter const & b)
{
    if (&seq_m == &b.seq_m && n_m==b.n_m) {
        return i_m <= b.i_m;
    } else {
        runtime_error e("Iter.seq_m or Iter.n_m is different.");
        throw e;
    }
    return false;
}
// }}}
