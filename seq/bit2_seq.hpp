#ifndef CZL_BIO_BIT2_SEQ_HPP
#define CZL_BIO_BIT2_SEQ_HPP

namespace czl_bio {

class Bit2Seq {
public:
    /**
     * @brief change A,C,G,T to two bit
     * @param  c  nucleic acid symbol (A,C,G,T)
     * @return  0x0 if c is 'A'
     *          0x1 if c is 'G'
     *          0x2 if c is 'T'
     *          0x3 if c is 'C'
     *          0x4 if c is other
     */
    static uint8_t nt1_to_bit2(char c);
    static char bit2_to_nt1(uint8_t b);
    static uint8_t nt1_to_bit2(char c);
    static uint8_t complement(uint8_t b);
    static Iter begin(size_t n, vector<uint8_t> & seq);
    static Iter end(size_t n, vector<uint8_t> & seq);

    class Iter {
    public:
        Iter(size_t n, vector<uint8_t> & seq);
        Iter(Iter const & src);
        Iter & operator =(Iter const & src);
        uint8_t operator *();
        Iter operator +(int64_t n);
        Iter operator -(int64_t n);
        int64_t operator -(Iter & iter);
        Iter & operator ++();
        Iter operator ++(Iter & iter);
        Iter & operator --();
        Iter operator --(Iter & iter);
        Iter & operator +=(int64_t n);
        Iter & operator -=(int64_t n);
        bool operator ==(Iter const & b);
        bool operator !=(Iter const & b);
        bool operator >(Iter const & b);
        bool operator <(Iter const & b);
        bool operator >=(Iter const & b);
        bool operator <=(Iter const & b);
    private:
        int64_t n_m;
        vector<uint8_t> & seq_m;
        int64_t i_m;
    };
};

};

#endif
