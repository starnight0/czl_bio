#ifndef CZL_BIO_BIT2_SEQ_HPP
#define CZL_BIO_BIT2_SEQ_HPP

#include "../czl_common.hpp"

namespace czl_bio {

class Bit2Seq {
public:
	class Iter;
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
    static uint8_t complement(uint8_t b);
	int nt1_to_bit2(string const & nt, size_t *n, vector<uint8_t> & b);

    Bit2Seq(): n_m(0) {}
    Bit2Seq(string const & nt);
    Bit2Seq(Bit2Seq const & src);
    ~Bit2Seq() {}

    int copy_from_nt1(string const & nt);
    int copy_to_nt1(string & nt);
    uint8_t get(size_t i) const;
    uint8_t get4(size_t i) const;
	inline uint8_t get_nothrow(size_t i) const;
    void set(size_t i, uint8_t b);
    void set4(size_t i, uint8_t b);
    inline void set_nothrow(size_t i, uint8_t b);
    int complement();
    int reverse();
    int rev_compl();
	int apply_fun1(int (*pf)(void*));
	int subseq(Bit2Seq & b, size_t pos, size_t len=(size_t)(-1));
	bool empty() { return n_m==0; }
	size_t size() { return n_m; }
	void resize(size_t n);
    Iter begin();
    Iter end();

	Bit2Seq & operator =(Bit2Seq const & src);
	Bit2Seq operator +(Bit2Seq const & src);
	Bit2Seq & operator +=(Bit2Seq const & src);

	bool operator ==(Bit2Seq const & src);
	bool operator !=(Bit2Seq const & src);
	bool operator <(Bit2Seq const & src);
	bool operator >(Bit2Seq const & src);

    class Iter {
    public:
        Iter(Bit2Seq * seq);
        Iter(Iter const & src);
        Iter & operator =(Iter const & src);
        uint8_t operator *();
        Iter operator +(int64_t n);
        Iter operator -(int64_t n);
        int64_t operator -(Iter & iter);
        Iter & operator ++();
        Iter operator ++(int);
        Iter & operator --();
        Iter operator --(int);
        Iter & operator +=(int64_t n);
        Iter & operator -=(int64_t n);
        bool operator ==(Iter const & b);
        bool operator !=(Iter const & b);
        bool operator >(Iter const & b);
        bool operator <(Iter const & b);
        bool operator >=(Iter const & b);
        bool operator <=(Iter const & b);
    private:
		Bit2Seq * bseq_m;
        int64_t i_m;
    };

private:
    size_t n_m;
    vector<uint8_t> seq_m;
};

};

#endif
