#include "rsa.h"
#include <string.h>
#include <ctype.h>

/* ---- 내부 상수/헬퍼 ---- */
#define LIMB_BITS 32u
#define MAX_LIMBS BIGNUM_ARRAY_SIZE

static inline void bn_normalize(Bignum* a) {
    int s = (a->size < 0) ? 0 : a->size;
    if (s > MAX_LIMBS) s = MAX_LIMBS;
    while (s > 0 && a->limbs[s - 1] == 0) s--;
    a->size = s;
}
static inline void bn_zero(Bignum* a) {
    memset(a->limbs, 0, sizeof(a->limbs));
    a->size = 0;
}
static inline int bn_is_zero(const Bignum* a) { return a->size == 0; }
static inline int bn_is_odd(const Bignum* a) { return (a->size > 0) && (a->limbs[0] & 1u); }
static inline int bn_ucmp(const Bignum* a, const Bignum* b) { // unsigned compare
    if (a->size != b->size) return (a->size > b->size) ? 1 : -1;
    for (int i = a->size - 1; i >= 0; --i) {
        if (a->limbs[i] != b->limbs[i]) return (a->limbs[i] > b->limbs[i]) ? 1 : -1;
    }
    return 0;
}
static inline uint32_t u32_min(uint32_t a, uint32_t b){ return a < b ? a : b; }

// 비트 길이 
static int bn_bit_length(const Bignum* a) {
    if (a->size == 0) return 0;
    uint32_t ms = a->limbs[a->size - 1];
    int bits = (a->size - 1) * (int)LIMB_BITS;
    // 상위 limb의 실제 비트수
    int leading = 32;
    while (leading > 0 && ((ms >> (leading - 1)) & 1u) == 0u) leading--;
    return bits + leading;
}
// i번째 비트(0=LSB) 조회 
static int bn_get_bit(const Bignum* a, int bit_index) {
    if (bit_index < 0) return 0;
    int limb = bit_index / (int)LIMB_BITS;
    int off  = bit_index % (int)LIMB_BITS;
    if (limb >= a->size) return 0;
    return (int)((a->limbs[limb] >> off) & 1u);
}
// q의 i번째 비트를 1로 set (0→1) 
static void bn_set_bit(Bignum* a, int bit_index) {
    int limb = bit_index / (int)LIMB_BITS;
    int off  = bit_index % (int)LIMB_BITS;
    if (limb >= MAX_LIMBS) return; // 초과분 무시
    a->limbs[limb] |= (1u << off);
    if (a->size <= limb) a->size = limb + 1;
}
// 1비트 왼쪽 시프트 (a <<= 1)
static void bn_shl1(Bignum* a) {
    uint32_t carry = 0;
    for (int i = 0; i < a->size; ++i) {
        uint64_t v = ((uint64_t)a->limbs[i] << 1) | carry;
        a->limbs[i] = (uint32_t)v;
        carry = (uint32_t)(v >> 32);
    }
    if (carry && a->size < MAX_LIMBS) {
        a->limbs[a->size++] = carry;
    }
}
// 1비트 오른쪽 시프트 (a >>= 1) 
static void bn_shr1(Bignum* a) {
    uint32_t carry = 0;
    for (int i = a->size - 1; i >= 0; --i) {
        uint32_t new_carry = a->limbs[i] & 1u;
        a->limbs[i] = (a->limbs[i] >> 1) | (carry << 31);
        carry = new_carry;
        if (i == 0) break;
    }
    bn_normalize(a);
}
// r = a + b, carry 반환 
static uint32_t bn_uadd(Bignum* r, const Bignum* a, const Bignum* b) {
    uint64_t carry = 0;
    int n = (a->size > b->size) ? a->size : b->size;
    if (n > MAX_LIMBS) n = MAX_LIMBS;
    for (int i = 0; i < n; ++i) {
        uint64_t av = (i < a->size) ? a->limbs[i] : 0;
        uint64_t bv = (i < b->size) ? b->limbs[i] : 0;
        uint64_t s = av + bv + carry;
        r->limbs[i] = (uint32_t)s;
        carry = s >> 32;
    }
    r->size = n;
    if (carry && r->size < MAX_LIMBS) {
        r->limbs[r->size++] = (uint32_t)carry;
        carry = 0;
    }
    return (uint32_t)carry;
}
// r = a - b (a>=b 가정), borrow 반환(0 정상, 1 underflow) 
static uint32_t bn_usub(Bignum* r, const Bignum* a, const Bignum* b) {
    if (bn_ucmp(a,b) < 0) { // underflow
        bn_zero(r);
        return 1;
    }
    uint64_t borrow = 0;
    int n = a->size;
    for (int i = 0; i < n; ++i) {
        uint64_t av = a->limbs[i];
        uint64_t bv = (i < b->size) ? b->limbs[i] : 0;
        uint64_t d = av - bv - borrow;
        r->limbs[i] = (uint32_t)d;
        borrow = (d >> 63) & 1u; // 음수면 borrow=1
    }
    r->size = n;
    bn_normalize(r);
    return (uint32_t)borrow;
}
// r += small(0..2^32-1) 
static void bn_add_small(Bignum* r, uint32_t x) {
    uint64_t carry = x;
    int i = 0;
    while (carry != 0 && i < MAX_LIMBS) {
        uint64_t v = (uint64_t)r->limbs[i] + carry;
        r->limbs[i] = (uint32_t)v;
        carry = v >> 32;
        i++;
    }
    if (i > r->size) r->size = i;
    if (carry && r->size < MAX_LIMBS) {
        r->limbs[r->size++] = (uint32_t)carry;
    }
}
// r *= small(0..2^32-1) 
static void bn_mul_small(Bignum* r, uint32_t k) {
    if (k == 0 || bn_is_zero(r)) { bn_zero(r); return; }
    uint64_t carry = 0;
    for (int i = 0; i < r->size; ++i) {
        uint64_t v = (uint64_t)r->limbs[i] * k + carry;
        r->limbs[i] = (uint32_t)v;
        carry = v >> 32;
    }
    if (carry && r->size < MAX_LIMBS) {
        r->limbs[r->size++] = (uint32_t)carry;
    }
}

// r = (a + b) mod m, 가정: a<m, b<m
static void bn_mod_add(Bignum* r, const Bignum* a, const Bignum* b, const Bignum* m) {
    Bignum t; bn_zero(&t);
    bn_uadd(&t, a, b); // t = a + b
    // t < 2m 이므로 한 번만 감산 검사
    if (bn_ucmp(&t, m) >= 0) {
        bn_usub(&t, &t, m);
    }
    *r = t;
}
// r = (2a) mod m, 가정: a<m 
static void bn_mod_double(Bignum* r, const Bignum* a, const Bignum* m) {
    Bignum t; bn_zero(&t);
    uint32_t carry = 0;
    int n = a->size;
    for (int i = 0; i < n; ++i) {
        uint64_t v = ((uint64_t)a->limbs[i] << 1) | carry;
        t.limbs[i] = (uint32_t)v;
        carry = (uint32_t)(v >> 32);
    }
    t.size = n;
    if (carry && t.size < MAX_LIMBS) t.limbs[t.size++] = carry;
    // t < 2m → 한 번 감산
    if (bn_ucmp(&t, m) >= 0) bn_usub(&t, &t, m);
    *r = t;
}

// a %= m (이진 나눗셈 기반, 몫은 버림)
static void bn_mod(Bignum* a, const Bignum* m) {
    if (bn_is_zero(m) || bn_is_zero(a)) return;
    if (bn_ucmp(a, m) < 0) return;

    Bignum r; bn_zero(&r);
    int nbits = bn_bit_length(a);
    for (int i = nbits - 1; i >= 0; --i) {
        // r = (r<<1) + bit_i(a)
        bn_shl1(&r);
        if (bn_get_bit(a, i)) bn_add_small(&r, 1u);
        // if r >= m: r -= m
        if (bn_ucmp(&r, m) >= 0) bn_usub(&r, &r, m);
    }
    *a = r;
}

/* ---- 공개 API 구현 ---- */

void bignum_init(Bignum* bn) { bn_zero(bn); }

void bignum_copy(Bignum* dest, const Bignum* src) {
    if (dest == src) return;
    memcpy(dest->limbs, src->limbs, sizeof(uint32_t) * MAX_LIMBS);
    dest->size = src->size;
}

// hex → bignum (대소문자/0x 허용). 성공 0, 실패 -1
int bignum_from_hex(Bignum* bn, const char* hex_str) {
    bn_zero(bn);
    if (!hex_str) return -1;

    // 앞쪽 공백/0x/0X 스킵
    const char* p = hex_str;
    while (*p && isspace((unsigned char)*p)) p++;
    if (p[0] == '0' && (p[1] == 'x' || p[1] == 'X')) p += 2;

    int seen = 0;
    for (; *p; ++p) {
        if (isspace((unsigned char)*p)) continue;
        int v;
        if (*p >= '0' && *p <= '9') v = *p - '0';
        else if (*p >= 'a' && *p <= 'f') v = 10 + (*p - 'a');
        else if (*p >= 'A' && *p <= 'F') v = 10 + (*p - 'A');
        else return -1; // invalid
        // bn = bn*16 + v
        bn_mul_small(bn, 16u);
        bn_add_small(bn, (uint32_t)v);
        seen = 1;
    }
    if (!seen) { bn_zero(bn); return 0; }
    bn_normalize(bn);
    return 0;
}

// bignum → hex (소문자, 선행 0 제거), 호출자가 free() 필요 
char* bignum_to_hex(const Bignum* bn) {
    if (bn->size == 0) {
        char* s = (char*)malloc(2);
        if (s) { s[0] = '0'; s[1] = '\0'; }
        return s;
    }
    // 최악 길이: size*8 + 1
    size_t maxlen = (size_t)bn->size * 8u + 1u;
    char* buf = (char*)malloc(maxlen);
    if (!buf) return NULL;

    // 최상위 limb는 선행 0 제거, 이후 limb는 8자리로 채움
    int i = bn->size - 1;
    int n = snprintf(buf, maxlen, "%x", bn->limbs[i]);
    for (i = bn->size - 2; i >= 0; --i) {
        n += snprintf(buf + n, maxlen - n, "%08x", bn->limbs[i]);
    }
    return buf;
}

// 비교: a>b → 1, a<b → -1, a==b → 0
int bignum_compare(const Bignum* a, const Bignum* b) {
    return bn_ucmp(a, b);
}

// result = a + b 
void bignum_add(Bignum* result, const Bignum* a, const Bignum* b) {
    bn_zero(result);
    bn_uadd(result, a, b);
    bn_normalize(result);
}

// result = a - b (a<b이면 0으로 클리어) 
void bignum_subtract(Bignum* result, const Bignum* a, const Bignum* b) {
    bn_zero(result);
    if (bn_ucmp(a, b) < 0) { bn_zero(result); return; }
    bn_usub(result, a, b);
}

// result = a * b 
void bignum_multiply(Bignum* result, const Bignum* a, const Bignum* b) {
    Bignum tmp;
    bignum_init(&tmp);

    if (a->size == 0 || b->size == 0) {
        bignum_init(result);
        return;
    }

    // 최대 유효 길이: a.size + b.size 
    int max_out = a->size + b->size;
    if (max_out > MAX_LIMBS) max_out = MAX_LIMBS;

    for (int i = 0; i < a->size; ++i) {
        uint64_t carry = 0;
        for (int j = 0; j < b->size; ++j) {
            int k = i + j;
            if (k >= MAX_LIMBS) {
                // 더 이상 하위 limb에 쓸 공간이 없으니, 남은 곱 결과는 잘림
                break;
            }
            uint64_t sum = (uint64_t)tmp.limbs[k]
                         + (uint64_t)a->limbs[i] * (uint64_t)b->limbs[j]
                         + carry;
            tmp.limbs[k] = (uint32_t)sum;
            carry = sum >> 32;
        }
        // 남은 carry를 다음 limb들로 전파
        int k = i + b->size;
        while (carry != 0 && k < MAX_LIMBS) {
            uint64_t sum = (uint64_t)tmp.limbs[k] + carry;
            tmp.limbs[k] = (uint32_t)sum;
            carry = sum >> 32;
            ++k;
        }
    }

    tmp.size = max_out;
    bn_normalize(&tmp);
    *result = tmp;
}

// 이진(long) 나눗셈: a = b*q + r (q,r 반환). b==0이면 q=0, r=a 
void bignum_divide(Bignum* quotient, Bignum* remainder, const Bignum* a, const Bignum* b) {
    if (quotient) bn_zero(quotient);
    if (remainder) bn_zero(remainder);
    if (bn_is_zero(b)) {
        if (remainder) bignum_copy(remainder, a);
        return;
    }
    if (bn_ucmp(a, b) < 0) {
        if (remainder) bignum_copy(remainder, a);
        return;
    }

    Bignum q; bn_zero(&q);
    Bignum r; bn_zero(&r);

    int nbits = bn_bit_length(a);
    for (int i = nbits - 1; i >= 0; --i) {
        // r = (r<<1) + bit_i(a)
        bn_shl1(&r);
        if (bn_get_bit(a, i)) bn_add_small(&r, 1u);
        // if r>=b: r-=b, q.setbit(i)
        if (bn_ucmp(&r, b) >= 0) {
            bn_usub(&r, &r, b);
            bn_set_bit(&q, i);
        }
    }
    if (quotient) *quotient = q;
    if (remainder) *remainder = r;
}

// (a*b) mod m — 러시아 농민법(russian peasant method): 2N 임시버퍼 없이 동작 
// 이 방법은 속도가 느려 만약 속도 문제시, Montgomery Reduction 등으로 교체 고려
static void bn_mod_mul(Bignum* r, const Bignum* a, const Bignum* b, const Bignum* m) {
    Bignum A, B, R;
    A = *a; B = *b; bn_zero(&R);

    // A,B를 m으로 축소(나머지)
    bn_mod(&A, m);
    bn_mod(&B, m);

    while (!bn_is_zero(&B)) {
        if (bn_is_odd(&B)) bn_mod_add(&R, &R, &A, m);
        bn_mod_double(&A, &A, m);
        bn_shr1(&B);
    }
    *r = R;
}

// 모듈러 거듭제곱: result = base^exp mod modulus 
void bignum_mod_exp(Bignum* result, const Bignum* base, const Bignum* exp, const Bignum* modulus) {
    // 특수 케이스
    if (bn_is_zero(modulus)) { bn_zero(result); return; }

    Bignum a, e, r;
    a = *base; e = *exp; bn_zero(&r);

    // r = 1
    r.size = 1; r.limbs[0] = 1u;

    // a %= modulus
    bn_mod(&a, modulus);

    while (!bn_is_zero(&e)) {
        if (bn_is_odd(&e)) {
            Bignum t; bn_zero(&t);
            bn_mod_mul(&t, &r, &a, modulus); // r = r*a mod m
            r = t;
        }
        // a = a*a mod m
        Bignum t2; bn_zero(&t2);
        bn_mod_mul(&t2, &a, &a, modulus);
        a = t2;

        // e >>= 1
        bn_shr1(&e);
    }
    *result = r;
}
