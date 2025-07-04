#ifdef DEBUG
#include <cstdio>
#endif
#include "vlong.hpp"

uint32_t vlong_value::get(uint32_t i) const {
    if (i >= n) return 0;
    return a[i];
}

void vlong_value::clear() {
    n = 0;
}

vlong_value::vlong_value() {
    z = 0;
    a = 0;
    n = 0;
    share = 0;
}

vlong_value::~vlong_value() {
    uint32_t i = z;
    while (i) { i -= 1; a[i] = 0; } // burn
    delete [] a;
}

void vlong_value::reserve(uint32_t x) {
    if (x > z) {
        uint32_t *na = new uint32_t[x];
        for (uint32_t i = 0; i < n; ++i) na[i] = a[i];
        delete [] a;
        a = na;
        z = x;
    }
}

void vlong_value::set(uint32_t i, uint32_t x) {
    if (i < n) {
        a[i] = x;
        if (x == 0) while (n && a[n-1] == 0) n -= 1; // normalise
    }
    else if (x) {
        reserve(i+1);
        for (uint32_t j = n; j < i; ++j) a[j] = 0;
        a[i] = x;
        n = i+1;
    }
}

vlong_value::operator unsigned() {
    return get(0);
}

int vlong_value::is_zero() const {
    return n == 0;
}

int vlong_value::test(uint32_t i) const {
    return (get(i / BPU) & (1 << (i % BPU))) != 0;
}

uint32_t vlong_value::bits() const {
    uint32_t x = n * BPU;
    while (x && test(x-1) == 0) x -= 1;
    return x;
}

void vlong_value::init(uint32_t x) {
    clear();
    set(0, x);
}

void vlong_value::copy(vlong_value& x) {
    clear();
    uint32_t i = x.n;
    while (i) {
        i -= 1;
        set(i, x.get(i));
    }
}

int vlong_value::cf(vlong_value& x) const {
    if (n > x.n) return +1;
    if (n < x.n) return -1;
    uint32_t i = n;
    while (i) {
        i -= 1;
        if (a[i] > x.a[i]) return +1;
        if (a[i] < x.a[i]) return -1;
    }
    return 0;
}

void vlong_value::add(vlong_value& x) {
    uint32_t carry = 0;
    uint32_t max, u, ux;
    max = (n > x.n) ? n : x.n;
    reserve(max);
    for (uint32_t i = 0; i <= max; ++i) {
        u = get(i) + carry;
        carry = (u < carry);
        ux = x.get(i);
        u = u + ux;
        carry += (u < ux);
        set(i, u);
    }
}

void vlong_value::sub(vlong_value& x) {
    uint32_t carry = 0;
    uint32_t u, ux, nu;
    for (uint32_t i = 0; i < n; ++i) {
        ux = x.get(i);
        ux += carry;
        if (ux >= carry) {
            u = get(i);
            nu = u - ux;
            carry = (nu > u);
            set(i, nu);
        }
    }
}

void vlong_value::shr() {
    uint32_t carry = 0;
    uint32_t i = n, u;
    while (i) {
        i -= 1;
        u = get(i);
        set(i, (u >> 1) + carry);
        carry = u << (BPU-1);
    }
}

void vlong_value::shr(uint32_t x) {
    uint32_t delta = x/BPU, u; x %= BPU;
    for (uint32_t i = 0; i < n; ++i) {
        u = get(i+delta);
        if (x) {
            u >>= x;
            u += get(i+delta+1) << (BPU-x);
        }
        set(i,u);
    }
}

void vlong_value::shl() {
    uint32_t carry = 0;
    uint32_t N = n, u; // necessary, since n can change
    for (uint32_t i = 0; i <= N; ++i) {
        u = get(i);
        set(i, (u << 1) + carry);
        carry = u >> (BPU-1);
    }
}

void vlong_value::fast_mul(vlong_value& x, vlong_value& y, uint32_t keep) {
    // *this = (x*y) % (2**keep)
    uint32_t i, j, limit = (keep + BPU - 1) / BPU; // size of result in words
    uint32_t m, c, v, p, w;
    reserve(limit);
    for (i = 0; i < limit; ++i) a[i] = 0;
    uint32_t minx = (x.n > limit) ? limit : x.n, miny;
    for (i = 0; i < minx; ++i) {
        m = x.a[i];
        c = 0; // carry
        miny = i + y.n;
        if (miny > limit) miny = limit;
        for (j = i; j < miny; ++j) {
            // This is the critical loop
            v = a[j], p = y.a[j-i];
            v += c; c = ( v < c );
            w = lo(p)*lo(m); v += w; c += ( v < w );
            w = lo(p)*hi(m); c += hi(w); w = lh(w); v += w; c += ( v < w );
            w = hi(p)*lo(m); c += hi(w); w = lh(w); v += w; c += ( v < w );
            c += hi(p) * hi(m);
            a[j] = v;
        }
        while (c && j < limit) {
            a[j] += c; 
            c = (a[j] < c);
            j += 1;
        }
    }

    // eliminate unwanted bits
    keep %= BPU; if (keep) a[limit-1] &= (1 << keep) - 1;

    // calculate n
    while (limit && a[limit-1] == 0) limit -= 1;
    n = limit;
}

void vlong_value::copy_range(const vlong_value& x, uint32_t start, uint32_t end) {
    clear();
    if (start >= end || start >= x.n) return;

    reserve(end - start);
    for (uint32_t i = start; i < end && i < x.n; ++i)
        set(i - start, x.get(i));
}

void vlong_value::add_shifted(const vlong_value& x, uint32_t shift) {
    uint32_t carry = 0, u, ux, sum;
    for (uint32_t i = 0; i < x.n || carry; ++i) {
        ux = x.get(i);
        u = get(i + shift);
        sum = u + ux + carry;
        carry = (sum < u) || ((sum == u) && (ux || carry));
        set(i + shift, sum);
    }
}

void vlong_value::karatsuba_mul(vlong_value& x, vlong_value& y) {
    uint32_t n = (x.n >= y.n) ? x.n : y.n;

    if (n <= 32) {
        fast_mul(x, y, x.bits() + y.bits());
        return;
    }

    uint32_t m = n >> 1;

    // x = x1·B^m + x0, y = y1·B^m + y0
    vlong_value x0, x1, y0, y1;
    x0.copy_range(x, 0, m);
    x1.copy_range(x, m, x.n);
    y0.copy_range(y, 0, m);
    y1.copy_range(y, m, y.n);

    // z0 = x0 * y0
    vlong_value z0; z0.karatsuba_mul(x0, y0);

    // z2 = x1 * y1
    vlong_value z2; z2.karatsuba_mul(x1, y1);

    // (x0 + x1), (y0 + y1)
    vlong_value xsum, ysum;
    xsum.copy(x0); xsum.add(x1);
    ysum.copy(y0); ysum.add(y1);

    // z1 = (x0 + x1)*(y0 + y1) - z0 - z2
    vlong_value z1; z1.karatsuba_mul(xsum, ysum);
    z1.sub(z0); z1.sub(z2);

    // z2·B^(2m) + z1·B^m + z0
    clear();
    add_shifted(z0, 0);
    add_shifted(z1, m);
    add_shifted(z2, m << 1);
}

void vlong_value::mul(vlong_value& x, vlong_value& y) {
    //fast_mul(x, y, x.bits()+y.bits());
    karatsuba_mul(x, y);
}

void vlong_value::div(vlong_value& x, vlong_value& y, vlong_value& rem) {
    init(0);
    rem.copy(x);
    vlong_value m, s;
    m.copy(y);
    s.init(1);
    while (rem.cf(m) > 0) {
        m.shl();
        s.shl();
    }
    while (rem.cf(y) >= 0) {
        while (rem.cf(m) < 0) {
            m.shr();
            s.shr();
        }
        rem.sub(m);
        add(s);
    }
}

void vlong::docopy() {
    if (value->share) {
        value->share -= 1;
        vlong_value *nv = new vlong_value;
        nv->copy(*value);
        value = nv;
    }
}

int vlong::cf(const vlong x) const {
    int neg = negative && !value->is_zero();
    if (neg == (x.negative && !x.value->is_zero()))
        return value->cf(*x.value);
    else if (neg) return -1;
    else return +1;
}

vlong::vlong(uint32_t x) {
    value = new vlong_value;
    negative = 0;
    value->init(x);
}

vlong::vlong(const vlong& x) { // copy constructor
    negative = x.negative;
    value = x.value;
    value->share += 1;
}

vlong& vlong::operator =(const vlong& x) {
    if (value->share) value->share -= 1; else delete value;
    value = x.value;
    value->share += 1;
    negative = x.negative;
    return *this;
}

vlong::~vlong() {
    if (value->share) value->share -= 1; else delete value;
}

vlong::operator unsigned () { // conversion to unsigned
    return value->get(0);
}

vlong& vlong::operator +=(const vlong& x) {
    if (negative == x.negative) {
        docopy();
        value->add(*x.value);
    }
    else if (value->cf(*x.value) >= 0) {
        docopy();
        value->sub(*x.value);
    }
    else {
        vlong tmp = *this;
        *this = x;
        *this += tmp;
    }
    return *this;
}

vlong& vlong::operator -=(const vlong& x) {
    if (negative != x.negative) {
        docopy();
        value->add(*x.value);
    }
    else if (value->cf(*x.value) >= 0) {
        docopy();
        value->sub(*x.value);
    }
    else {
        vlong tmp = *this;
        *this = x;
        *this -= tmp;
        negative = 1 - negative;
    }
    return *this;
}

vlong operator +(const vlong& x, const vlong& y) {
    vlong result = x;
    result += y;
    return result;
}

vlong operator -(const vlong& x, const vlong& y) {
    vlong result = x;
    result -= y;
    return result;
}

vlong operator *(const vlong& x, const vlong& y) {
    vlong result;
    result.value->mul(*x.value, *y.value);
    result.negative = x.negative ^ y.negative;
    return result;
}

vlong operator /(const vlong& x, const vlong& y) {
    vlong result;
    vlong_value rem;
    result.value->div(*x.value, *y.value, rem);
    result.negative = x.negative ^ y.negative;
    return result;
}

vlong operator %(const vlong& x, const vlong& y) {
    vlong result;
    vlong_value divide;
    divide.div(*x.value, *y.value, *result.value);
    result.negative = x.negative; // not sure about this?
    return result;
}

vlong gcd(const vlong& X, const vlong& Y) {
    vlong x = X, y = Y;
    while (1) {
        if ( y == (vlong)0 ) return x;  // ambiguity fixed.
        x = x % y;
        if ( x == (vlong)0 ) return y;
        y = y % x;
    }
}

vlong modinv(const vlong& a, const vlong& m) // modular inverse
// returns i in range 1..m-1 such that i*a = 1 mod m
// a must be in range 1..m-1
{
    vlong j = 1, i = 0, b = m, c = a, x, y;
    while (c != (vlong)0) {
        x = b / c;
        y = b - x*c;
        b = c;
        c = y;
        y = j;
        j = i - j*x;
        i = y;
    }
    if (i < (vlong)0)
        i += m;
    return i;
}

#ifdef DEBUG
void vlong_value::print() {
    int point = n - 1;
    while (*(a + point) == 0) point--;
    printf("%X", a[point--]);
    while (point >= 0) {
        printf("%08X", a[point]);
        point--;
    }
    puts("");
}
#endif

class monty { // class for montgomery modular exponentiation
    vlong R, R1, m, n1;
    vlong T, k;   // work registers
    uint32_t N;   // bits for R
    void mul(vlong& x, const vlong& y);
public:
    vlong exp(const vlong& x, const vlong& e);
    monty(const vlong& M);
};

monty::monty(const vlong& M) {
    m = M;
    N = 0; R = 1; while ( R < M ) { R += R; N += 1; }
    R1 = modinv( R-m, m );
    n1 = R - modinv( m, R );
}

void monty::mul(vlong& x, const vlong& y) {
    // T = x*y;
    T.value->fast_mul(*x.value, *y.value, N*2);

    // k = ( T * n1 ) % R;
    k.value->fast_mul(*T.value, *n1.value, N);

    // x = ( T + k*m ) / R;
    x.value->fast_mul(*k.value, *m.value, N*2);
    x += T;
    x.value->shr(N);

    if (x >= m) x -= m;
}

vlong monty::exp(const vlong& x, const vlong& e) {
    vlong result = R-m, t = ( x * R ) % m;
    uint32_t bits = e.value->bits();
    uint32_t i = 0;
    while (1) {
        if (e.value->test(i))
            mul(result, t);
        i += 1;
        if (i == bits) break;
        mul(t, t);
    }
    return (result * R1) % m;
}

vlong modexp(const vlong& x, const vlong& e, const vlong& m) {
    monty me(m);
    return me.exp(x, e);
}


