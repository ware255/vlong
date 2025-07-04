/*
 *
 * Array orientation is little endian
 * a[] = { 0x01230000, 0x45670000, 0x89AB0000, 0xCDFF0000 };
 * print: CDEF000089AB00004567000001230000
 *
 * Relationship between arrays and the number of bits
 * +--------+------+
 * | Arrays | Bits |
 * +--------+------+
 * |      1 |   32 |
 * |      2 |   64 |
 * |      4 |  128 |
 * |      8 |  256 |
 * |     16 |  512 |
 * |     32 | 1024 |
 * |     64 | 2048 |
 * |    128 | 4096 |
 * +--------+------+
 *
 * The class is divided into two for memory optimization.
 * 
 */

#ifndef _vlong_hpp_
#define _vlong_hpp_

// Macros for doing double precision multiply
#define BPU   ( 8*sizeof(unsigned int) ) // Number of bits in an unsigned
#define lo(x) ( (x) & ((1<<(BPU/2))-1) ) // lower half of unsigned
#define hi(x) ( (x) >> (BPU/2) )         // upper half
#define lh(x) ( (x) << (BPU/2) )         // make upper half

typedef unsigned int uint32_t;

class vlong_value {
    uint32_t * a;                        // array of units
    uint32_t z;                          // units allocated
    uint32_t n;                          // used units (read-only)
    void clear();                        // set n to zero
    void reserve(uint32_t x);            // storage hint
    void copy_range(const vlong_value& x, uint32_t start, uint32_t end);
    void add_shifted(const vlong_value& x, uint32_t shift);
public:
    uint32_t share; // share count, used by vlong to delay physical copying
    uint32_t get(uint32_t i) const;      // get ith unsigned
    void set(uint32_t i, uint32_t x);    // set ith unsigned
    int is_zero() const;
    int test(uint32_t i) const;
    uint32_t bits() const;
    void init(uint32_t x);
    void copy(vlong_value& x);
    int cf(vlong_value& x) const;
    void shl();
    void shr();
    void shr(uint32_t x);
    void add(vlong_value& x);
    void sub(vlong_value& x);
    operator unsigned();                 // Unsafe conversion to unsigned
    vlong_value();
    ~vlong_value();
    // Time critical routine
    void fast_mul(vlong_value& x, vlong_value& y, uint32_t keep);
    void karatsuba_mul(vlong_value& x, vlong_value& y);
    void mul(vlong_value& x, vlong_value& y);
    void div(vlong_value& x, vlong_value& y, vlong_value& rem);
};

class vlong { // very long integer - can be used like long
    int cf(const vlong x) const;
    void docopy();
public:
    int negative;
    class vlong_value *value;

    // Standard arithmetic operators
    friend vlong operator +( const vlong& x, const vlong& y );
    friend vlong operator -( const vlong& x, const vlong& y );
    friend vlong operator *( const vlong& x, const vlong& y );
    friend vlong operator /( const vlong& x, const vlong& y );
    friend vlong operator %( const vlong& x, const vlong& y );
    vlong& operator +=( const vlong& x );
    vlong& operator -=( const vlong& x );

    // Standard comparison operators
    friend inline int operator !=( const vlong& x, const vlong& y ){ return x.cf( y ) != 0; }
    friend inline int operator ==( const vlong& x, const vlong& y ){ return x.cf( y ) == 0; }
    friend inline int operator >=( const vlong& x, const vlong& y ){ return x.cf( y ) >= 0; }
    friend inline int operator <=( const vlong& x, const vlong& y ){ return x.cf( y ) <= 0; }
    friend inline int operator > ( const vlong& x, const vlong& y ){ return x.cf( y ) > 0; }
    friend inline int operator < ( const vlong& x, const vlong& y ){ return x.cf( y ) < 0; }

    // Construction and conversion operations
    vlong(uint32_t x = 0);
    vlong(const vlong& x);
    ~vlong();
    operator unsigned ();
    vlong& operator =(const vlong& x);
};

vlong modexp(const vlong& x, const vlong& e, const vlong& m); // m must be odd
vlong gcd(const vlong& X, const vlong& Y); // greatest common denominator
vlong modinv(const vlong& a, const vlong& m); // modular inverse

#include "vlong.cpp"

#endif // _vlong_hpp_
