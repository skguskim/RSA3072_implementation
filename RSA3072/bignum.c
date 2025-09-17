#include "bignum.h"
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdlib.h> // for malloc, free, calloc
#include <stdio.h>  // for snprintf, printf

#define LIMB_BITS 32u
#define MAX_LIMBS BIGNUM_ARRAY_SIZE

typedef unsigned char  u8;
typedef uint32_t       u32;
typedef uint64_t       u64;

// =====================================================================
//  공통 헬퍼 함수
// =====================================================================

// 상위 0 워드들을 잘라 유효 길이로 줄여 반환
static inline int  limbs_trim(const u32* x, int n){ 
    while(n>0 && x[n-1]==0u) --n; 
    return n; 
}

// 비트 길이 확인
// a: 비트 길이를 확인할 큰 수 포인터
// 반환값: 비트 길이
static int bn_bit_length(const Bignum* a) {
    if (a->size == 0) return 0;
    uint32_t ms = a->limbs[a->size - 1];
    int bits = (a->size - 1) * (int)LIMB_BITS;
    int leading = 32;
    while (leading > 0 && ((ms >> (leading - 1)) & 1u) == 0u) leading--;
    return bits + leading;
}

// 큰 수 a를 0으로 변경
// a: 0으로 바꿀 큰 수 포인터
// 반환값: 없음
static void bn_zero(Bignum* a) {
    memset(a->limbs, 0, sizeof(uint32_t) * MAX_LIMBS);
    a->size = 0;
}

// a가 0이면 true, 아니면 false 반환
static inline int bn_is_zero(const Bignum* a) { return a->size == 0; }

// 큰 수 a의 크기를 정규화(예시: 0032 -> 32)
// a: 정규화할 큰 수 포인터
// 반환값: 없음
static void bn_normalize(Bignum* a) {
    while (a->size > 0 && a->limbs[a->size - 1] == 0) {
        a->size--;
    }
}

// i번째 비트(0 = LSB) 조회 
// 반환값: 조회한 비트 값(0 또는 1)
static int bn_get_bit(const Bignum* a, const int bit_index) {
    if (bit_index < 0) return 0;
    int limb = bit_index / (int)LIMB_BITS;
    int off = bit_index % (int)LIMB_BITS;
    if (limb >= a->size) return 0;
    return (int)((a->limbs[limb] >> off) & 1u);
}

//  r <- r + x (여기서 x는 0..2^32-1의 작은 정수)
// 내부에서 64비트 누산으로 캐리를 처리하며, r가 충분히 크면 상위 워드로 캐리를 전파
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

//  r <- r * k (k는 0..2^32-1)
// 각 limb에 대해 곱하고 캐리를 전파하며 k==0 또는 r==0이면 r=0
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

// n비트 왼쪽 시프트 (a <<= n)
static void bn_shift_left(Bignum* a, const int shift) {
    if (shift <= 0 || a->size == 0) return;
    int limb_shift = shift / (int)LIMB_BITS;
    int bit_shift = shift % (int)LIMB_BITS;

    if (a->size + limb_shift >= MAX_LIMBS) {
        // 포화 방지: 가능한 한도까지만 밀고 나머지는 버림
        limb_shift = (MAX_LIMBS - 1) - a->size;
        if (limb_shift < 0) limb_shift = 0;
    }

    uint32_t tmp[MAX_LIMBS] = { 0 };
    uint32_t carry = 0;
    for (int i = 0; i < a->size && i + limb_shift < MAX_LIMBS; ++i) {
        uint32_t lo = (bit_shift == 0) ? a->limbs[i] : (a->limbs[i] << bit_shift);
        uint32_t hi = (bit_shift == 0) ? 0 : (a->limbs[i] >> (32 - bit_shift));
        uint64_t v = (uint64_t)lo + ((uint64_t)carry);
        tmp[i + limb_shift] = (uint32_t)v;
        carry = (uint32_t)(v >> 32);
        if (i + limb_shift + 1 < MAX_LIMBS) {
            // hi는 다음 워드로 더해짐
            uint64_t w = (uint64_t)tmp[i + limb_shift + 1] + hi;
            tmp[i + limb_shift + 1] = (uint32_t)w;
            if (w >> 32) carry += 1u;
        }
    }

    int new_size = a->size + limb_shift + 1;
    if (new_size > MAX_LIMBS) new_size = MAX_LIMBS;
    memcpy(a->limbs, tmp, sizeof(uint32_t) * new_size);
    a->size = new_size;
    bn_normalize(a);
}

// n비트 오른쪽 시프트 (a >>= n)
static void bn_shift_right(Bignum* a, const int shift) {
    if (shift <= 0 || a->size == 0) return;
    int limb_shift = shift / (int)LIMB_BITS;
    int bit_shift = shift % (int)LIMB_BITS;

    if (limb_shift >= a->size) { bn_zero(a); return; }

    uint32_t tmp[MAX_LIMBS] = { 0 };
    int j = 0;
    for (int i = limb_shift; i < a->size; ++i, ++j) {
        uint32_t cur = a->limbs[i];
        uint32_t prev = (i + 1 < a->size) ? a->limbs[i + 1] : 0;
        if (bit_shift == 0) {
            tmp[j] = cur;
        }
        else {
            tmp[j] = (cur >> bit_shift) | (prev << (32 - bit_shift));
        }
    }
    memcpy(a->limbs, tmp, sizeof(uint32_t) * j);
    a->size = j;
    bn_normalize(a);
}

// 카라츠바 배열 연산 헬퍼
static inline u32 add_n(u32* r, const u32* x, const u32* y, int n){
    u64 c=0;
    for(int i=0;i<n;++i){ u64 s=(u64)x[i]+(u64)y[i]+c; r[i]=(u32)s; c=s>>32; }
    return (u32)c;
}

// r = x + y (가변 길이), 반환: 실제 길이
static int add_var(u32* r, const u32* x, int xn, const u32* y, int yn){
    int n = (xn>yn?xn:yn);
    u64 c=0;
    for(int i=0;i<n;++i){
        u64 xi = (i<xn)?x[i]:0, yi=(i<yn)?y[i]:0;
        u64 s = xi+yi+c;
        r[i]=(u32)s; c=s>>32;
    }
    if(c){ r[n++]=(u32)c; }
    return n;
}

// r = x - y (가변 길이, x>=y 가정), 반환: 실제 길이
static int sub_var(u32* r, const u32* x, int xn, const u32* y, int yn){
    u64 b=0;
    for(int i=0;i<xn;++i){
        u64 xi=x[i], yi=(i<yn)?y[i]:0;
        u64 d = xi - yi - b;
        r[i]=(u32)d;
        b = (d>>63)&1u;
    }
    int rn = xn; while(rn>0 && r[rn-1]==0u) --rn;
    return rn;
}

// 스쿨북 곱: r 길이 >= an+bn, r는 사전에 0으로 클리어되어 있어야 안정적
static void mul_schoolbook(u32* r, const u32* a, int an, const u32* b, int bn){
    for(int i=0;i<an;++i){
        u64 carry=0, ai=a[i];
        int k=i;
        for(int j=0;j<bn;++j,++k){
            u64 t = (u64)r[k] + ai*(u64)b[j] + carry;
            r[k]=(u32)t; carry=t>>32;
        }
        // 남은 carry 전파
        while(carry){
            u64 t=(u64)r[k] + carry;
            r[k]=(u32)t; carry=t>>32; ++k;
        }
    }
}

// r += x (오프셋 off에서 시작)
static void add_into_at(u32* r, int rlen, const u32* x, int xn, int off){
    u64 c=0;
    int i=0, k=off;
    for(; i<xn; ++i, ++k){
        u64 t=(u64)r[k] + (u64)x[i] + c;
        r[k]=(u32)t; c=t>>32;
    }
    while(c && k<rlen){
        u64 t=(u64)r[k] + c;
        r[k]=(u32)t; c=t>>32; ++k;
    }
}
// =====================================================================
//  큰 수 초기화, 복사, 해제, 변환 등
// =====================================================================

// ========= API 함수 =========
// 큰 수 bn을 0으로 초기화(선언 시에만 사용)
void bignum_init(Bignum* bn) { bn_zero(bn); }

// 큰 수 bn을 0으로 변경
void bignum_set_zero(Bignum* bn) { bn_zero(bn); }

// 큰 수를 복사
// dest: 복사 대상
// src: 복사할 원본
// 반환값: 없음
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

void print_bignum(const char* name, const Bignum* bn) {
    char* hex = bignum_to_hex(bn);
    printf(" - %s: %s\n", name, hex);
    free(hex);
}

// =====================================================================
//  큰 수 비교
// =====================================================================

// ========= 헬퍼 함수 =========
static inline int bn_ucmp(const Bignum* a, const Bignum* b) {
    if (a->size != b->size) return (a->size > b->size) ? 1 : -1;
    for (int i = a->size - 1; i >= 0; --i) {
        if (a->limbs[i] != b->limbs[i]) return (a->limbs[i] > b->limbs[i]) ? 1 : -1;
    }
    return 0;
}

// ========= API 함수 =========
// 비교: a>b → 1, a<b → -1, a==b → 0
int bignum_compare(const Bignum* a, const Bignum* b) {
    return bn_ucmp(a, b);
}

// =====================================================================
//  큰 수 덧셈
// =====================================================================

// ========= 헬퍼 함수 =========
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

// ========= API 함수 =========
// result = a + b 
void bignum_add(Bignum* result, const Bignum* a, const Bignum* b) {
    bn_zero(result);
    bn_uadd(result, a, b);
    bn_normalize(result);
}

// =====================================================================
//  큰 수 뺄셈
// =====================================================================

// ========= 헬퍼 함수 =========
// r = a - b (a>=b 가정), borrow 반환(0 정상, 1 underflow) 
static uint32_t bn_usub(Bignum* r, const Bignum* a, const Bignum* b) {
    if (bn_ucmp(a, b) < 0) { // underflow
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

// ========= API 함수 =========
// result = a - b (a<b이면 0으로 클리어) 
void bignum_subtract(Bignum* result, const Bignum* a, const Bignum* b) {
    bn_zero(result);
    if (bn_ucmp(a, b) < 0) { bn_zero(result); return; }
    bn_usub(result, a, b);
}

// =====================================================================
//  큰 수 곱셈 (카라츠바 알고리즘)
// =====================================================================

// ========= 헬퍼 함수 =========
// r = a * b
static int add_var(u32* r, const u32* x, int xn, const u32* y, int yn){
    int n = (xn>yn?xn:yn);
    u64 c=0;
    for(int i=0;i<n;++i){
        u64 xi = (i<xn)?x[i]:0, yi=(i<yn)?y[i]:0;
        u64 s = xi+yi+c;
        r[i]=(u32)s; c=s>>32;
    }
    if(c){ r[n++]=(u32)c; }
    return n;
}

static int sub_var(u32* r, const u32* x, int xn, const u32* y, int yn){
    u64 b=0;
    for(int i=0;i<xn;++i){
        u64 xi=x[i], yi=(i<yn)?y[i]:0;
        u64 d = xi - yi - b;
        r[i]=(u32)d;
        b = (d>>63)&1u;
    }
    return limbs_trim(r, xn);
}

static void mul_schoolbook(u32* r, const u32* a, int an, const u32* b, int bn){
    if (an == 0 || bn == 0) return;
    memset(r, 0, sizeof(u32) * (an + bn));
    for(int i=0;i<an;++i){
        u64 carry=0, ai=a[i];
        for(int j=0;j<bn;++j){
            u64 t = (u64)r[i+j] + ai*(u64)b[j] + carry;
            r[i+j]=(u32)t; carry=t>>32;
        }
        if(carry) {
            r[i+bn] += (u32)carry;
        }
    }
}

static void add_into_at(u32* r, int rlen, const u32* x, int xn, int off){
    if (xn == 0) return;
    u64 c=0;
    int i=0, k=off;
    for(; i<xn && k<rlen; ++i, ++k){
        u64 t=(u64)r[k] + (u64)x[i] + c;
        r[k]=(u32)t; c=t>>32;
    }
    while(c && k<rlen){
        u64 t=(u64)r[k] + c;
        r[k]=(u32)t; c=t>>32; ++k;
    }
}

static void karatsuba_recursive(u32* r, const u32* a, int an, const u32* b, int bn, u32* ws) {
    an = limbs_trim(a, an);
    bn = limbs_trim(b, bn);

    const int KARATSUBA_CUTOFF = 32;
    if (an <= KARATSUBA_CUTOFF || bn <= KARATSUBA_CUTOFF || an == 0 || bn == 0) {
        if (an > 0 && bn > 0) {
            mul_schoolbook(r, a, an, b, bn);
        }
        return;
    }

    int m = (an > bn) ? (an / 2) : (bn / 2);
    
    const u32* a0 = a; int an0 = (an > m) ? m : an;
    const u32* a1 = a + m; int an1 = (an > m) ? (an - m) : 0;
    const u32* b0 = b; int bn0 = (bn > m) ? m : bn;
    const u32* b1 = b + m; int bn1 = (bn > m) ? (bn - m) : 0;

    u32* sA = ws;
    u32* sB = sA + (m + 1);
    u32* z1_prod = sB + (m + 1);
    u32* next_ws = z1_prod + (2 * m + 2);

    u32* z0 = r;
    u32* z2 = r + 2 * m;
    memset(z0, 0, sizeof(u32) * 2 * m);
    memset(z2, 0, sizeof(u32) * (an - m + bn - m + 2)); 
    karatsuba_recursive(z0, a0, an0, b0, bn0, next_ws);
    karatsuba_recursive(z2, a1, an1, b1, bn1, next_ws);
    
    int z0n = limbs_trim(z0, 2 * m);
    int z2n = limbs_trim(z2, an - m + bn - m + 2); 
    
    int sAn = add_var(sA, a0, an0, a1, an1);
    int sBn = add_var(sB, b0, bn0, b1, bn1);
    memset(z1_prod, 0, sizeof(u32) * (sAn + sBn + 1));
    karatsuba_recursive(z1_prod, sA, sAn, sB, sBn, next_ws);
    
    int z1n = limbs_trim(z1_prod, sAn + sBn);
    z1n = sub_var(z1_prod, z1_prod, z1n, z0, z0n);
    z1n = sub_var(z1_prod, z1_prod, z1n, z2, z2n);

    add_into_at(r, an + bn + 2, z1_prod, z1n, m);
    add_into_at(r, an + bn + 2, z2, z2n, 2 * m);
}

// ========= API 함수 =========
// result = a * b
void bignum_multiply(Bignum* result, const Bignum* a, const Bignum* b) {
    if (!result || !a || !b) return;
    if (bn_is_zero(a) || bn_is_zero(b)) { 
        bn_zero(result); 
        return; 
    }

    int an = a->size;
    int bn = b->size;
    int rn = an + bn;

    u32* rbuf = (u32*)calloc(rn + 2, sizeof(u32));
    if (!rbuf) { bn_zero(result); return; }

    int n = (an > bn ? an : bn);
    size_t ws_size = 4 * n + 8;
    u32* workspace = (u32*)calloc(ws_size, sizeof(u32));
    if (!workspace) { free(rbuf); bn_zero(result); return; }
    
    karatsuba_recursive(rbuf, a->limbs, an, b->limbs, bn, workspace);

    int trim = limbs_trim(rbuf, rn + 2);
    if (trim > MAX_LIMBS) trim = MAX_LIMBS;
    
    bn_zero(result);
    memcpy(result->limbs, rbuf, sizeof(u32) * trim);
    result->size = trim;

    free(rbuf);
    free(workspace);
}

// =====================================================================
//  큰 수 나눗셈
// =====================================================================
void bignum_divide(Bignum* quotient, Bignum* remainder,
    const Bignum* a, const Bignum* b) {
    bn_zero(quotient);
    bn_zero(remainder);

    // divide by zero: 조용히 리턴 (필요시 에러 처리)
    if (bn_is_zero(b) || bn_is_zero(a)) return;

    int cmp = bn_ucmp(a, b);
    if (cmp < 0) { // a < b
        bignum_copy(remainder, a);
        return;
    }
    if (cmp == 0) { // a == b
        quotient->limbs[0] = 1;
        quotient->size = 1;
        return;
    }

    // 1워드 divisor 빠른 경로
    if (b->size == 1) {
        uint32_t divisor = b->limbs[0];
        uint64_t rem = 0;

        quotient->size = a->size;
        for (int i = a->size - 1; i >= 0; --i) {
            uint64_t cur = (rem << 32) | a->limbs[i];
            quotient->limbs[i] = (uint32_t)(cur / divisor);
            rem = cur % divisor;
        }
        bn_normalize(quotient);

        if (rem) { remainder->limbs[0] = (uint32_t)rem; remainder->size = 1; }
        else { remainder->size = 0; }
        return;
    }

    // Knuth Algorithm D

    // 정규화: 제수 최상위 워드에 leading 1이 오도록 왼쪽으로 s비트 시프트
    Bignum U, V;
    bignum_copy(&U, a);
    bignum_copy(&V, b);

    uint32_t v_high = V.limbs[V.size - 1];
    int s = 0;
    if (v_high < 0x80000000u) {
        while ((v_high & 0x80000000u) == 0u) { v_high <<= 1; s++; }
    }
    if (s) {
        bn_shift_left(&U, s);
        bn_shift_left(&V, s);
    }

    int n = V.size;                   // length of divisor
    int m = U.size - n;               // 몫 자리수-1 (j는 m..0)
    if (m < 0) {                      // 방어
        bignum_copy(remainder, a);
        return;
    }

    // 배열로 작업 (Up: 길이 n+m+1, Vp: 길이 n)
    uint32_t Up[MAX_LIMBS + 1] = { 0 };
    uint32_t Vp[MAX_LIMBS] = { 0 };

    for (int i = 0; i < U.size && i < MAX_LIMBS + 1; ++i) Up[i] = U.limbs[i];
    for (int i = 0; i < n && i < MAX_LIMBS; ++i)          Vp[i] = V.limbs[i];

    // 몫 워드 버퍼
    uint32_t Qp[MAX_LIMBS] = { 0 };
    int q_len = m + 1;
    if (q_len > MAX_LIMBS) q_len = MAX_LIMBS;

    // 메인 루프: j = m down to 0
    for (int j = m; j >= 0; --j) {
        // D3: qhat, rhat 추정 (상위 두 워드 사용)
        uint64_t u2 = Up[j + n];          // 상위 워드
        uint64_t u1 = Up[j + n - 1];
        uint64_t u0 = (n >= 2) ? Up[j + n - 2] : 0;

        uint64_t v1 = Vp[n - 1];
        uint64_t v0 = (n >= 2) ? Vp[n - 2] : 0;

        uint64_t uhat = (u2 << 32) | u1;
        uint64_t qhat = uhat / v1;
        uint64_t rhat = uhat % v1;
        if (qhat > 0xFFFFFFFFull) qhat = 0xFFFFFFFFull;

        // 과추정 보정
        if (n > 1) {
            while (qhat * v0 > ((rhat << 32) | u0)) {
                qhat--;
                rhat += v1;
                if (rhat >= (1ull << 32)) break; // rhat overflow 시 중지
            }
        }

        // D4: Up[j..j+n] -= qhat * Vp[0..n-1]
        uint64_t borrow = 0;
        uint64_t carrymul = 0;
        for (int i = 0; i < n; ++i) {
            uint64_t p = qhat * Vp[i] + carrymul;
            carrymul = p >> 32;
            uint64_t diff = (uint64_t)Up[j + i] - (uint32_t)p - borrow;
            Up[j + i] = (uint32_t)diff;
            borrow = (diff >> 63) & 1u;
        }
        // 상위 워드까지 반영
        uint64_t diff = (uint64_t)Up[j + n] - carrymul - borrow;
        Up[j + n] = (uint32_t)diff;
        borrow = (diff >> 63) & 1u;

        // D5: borrow면 add back, qhat--
        if (borrow) {
            uint64_t carry2 = 0;
            for (int i = 0; i < n; ++i) {
                uint64_t ssum = (uint64_t)Up[j + i] + Vp[i] + carry2;
                Up[j + i] = (uint32_t)ssum;
                carry2 = ssum >> 32;
            }
            Up[j + n] += (uint32_t)carry2;
            qhat--;
        }

        // 몫 자리 기록
        if (j < q_len) Qp[j] = (uint32_t)qhat;
    }

    // 몫 복사
    quotient->size = q_len;
    for (int i = 0; i < q_len; ++i) quotient->limbs[i] = Qp[i];
    bn_normalize(quotient);

    // 나머지: Up[0..n-1] 를 비정규화 복구 (s비트 오른쪽 시프트)
    Bignum Rtmp; bn_zero(&Rtmp);
    Rtmp.size = n;
    if (Rtmp.size > MAX_LIMBS) Rtmp.size = MAX_LIMBS;
    for (int i = 0; i < Rtmp.size; ++i) Rtmp.limbs[i] = Up[i];

    if (s) bn_shift_right(&Rtmp, s);
    bn_normalize(&Rtmp);
    bignum_copy(remainder, &Rtmp);
}


// =====================================================================
//  Montgomery Multiplication
// =====================================================================

// ========= 헬퍼 함수 =========
// 상수 시간 로드를 위한 마스크
static inline uint32_t ct_mask_if(int cond) {
    // cond != 0 -> 0xFFFFFFFF, cond == 0 -> 0x00000000
    return (uint32_t)0 - (uint32_t)(cond != 0);
}

// out = x - y  (n워드), 반환: borrow(0 또는 1). 분기 없음.
static uint32_t sub_n_ct(uint32_t* out, const uint32_t* x, const uint32_t* y, int n) {
    uint64_t borrow = 0;
    for (int i = 0; i < n; ++i) {
        uint64_t xi = x[i];
        uint64_t yi = y[i];
        uint64_t d = xi - yi - borrow;
        out[i] = (uint32_t)d;
        borrow = (d >> 63) & 1u;  // 음수면 borrow = 1
    }
    return (uint32_t)borrow;
}

// n0' 계산
static uint32_t mont_n0prime(uint32_t n0) {
    uint32_t x = 1;
    x *= 2 - n0 * x;
    x *= 2 - n0 * x;
    x *= 2 - n0 * x;
    x *= 2 - n0 * x;
    x *= 2 - n0 * x;
    return (uint32_t)(0u - x); // = -n0^{-1} mod 2^32
}

// CIOS 방법을 사용한 모듈러 거듭제곱 전용 몽고메리 곱셈
static void mont_mul(Bignum* r, const Bignum* a, const Bignum* b, const Bignum* N, uint32_t n0prime) {
    const int n = N->size;

    // N, a, b의 하위 n워드를 로컬에 준비 (상수시간 로드)
    uint32_t NN[MAX_LIMBS] = { 0 };
    uint32_t AA[MAX_LIMBS] = { 0 };
    uint32_t BB[MAX_LIMBS] = { 0 };

    for (int i = 0; i < n; ++i) {
        NN[i] = N->limbs[i]; // N은 모듈러스(공개)라 그냥 로드
        // a,b는 상수시간 마스크로 size 넘어가면 0 취급
        uint32_t am = ct_mask_if(i < a->size);
        uint32_t bm = ct_mask_if(i < b->size);
        AA[i] = a->limbs[i] & am;
        BB[i] = b->limbs[i] & bm;
    }

    // t: 길이 n+2 버퍼
    uint32_t t[MAX_LIMBS + 2] = { 0 };

    // CIOS 스캔 (고정 길이)
    for (int i = 0; i < n; ++i) {
        uint64_t bi = (uint64_t)BB[i];

        // t = t + a*bi
        uint64_t carry = 0;
        for (int j = 0; j < n; ++j) {
            uint64_t sum = (uint64_t)t[j] + (uint64_t)AA[j] * bi + carry;
            t[j] = (uint32_t)sum;
            carry = sum >> 32;
        }
        uint64_t acc = (uint64_t)t[n] + carry;
        t[n] = (uint32_t)acc;
        t[n + 1] = (uint32_t)(acc >> 32);

        // m = (t[0] * n0') mod 2^32  (하위워드만 필요)
        uint32_t m = (uint32_t)((uint64_t)t[0] * (uint64_t)n0prime);

        // t = t + m*N
        carry = 0;
        for (int j = 0; j < n; ++j) {
            uint64_t sum = (uint64_t)t[j] + (uint64_t)m * (uint64_t)NN[j] + carry;
            t[j] = (uint32_t)sum;
            carry = sum >> 32;
        }
        acc = (uint64_t)t[n] + carry;
        t[n] = (uint32_t)acc;
        t[n + 1] += (uint32_t)(acc >> 32);

        // t = t / b (워드 1칸 오른쪽 시프트) : 고정 길이
        for (int k = 0; k < n; ++k) t[k] = t[k + 1];
        t[n] = 0;
    }

    // 최종 보정
    uint32_t t_minus_N[MAX_LIMBS] = { 0 };
    uint32_t borrow = sub_n_ct(t_minus_N, t, NN, n); // t-N
    // borrow==0 → t>=N → t_minus_N 선택 / borrow==1 → t<N → t 선택
    uint32_t mask = (uint32_t)0 - (uint32_t)(1u - borrow); // 0xFFFFFFFF 또는 0x00000000

    // out = (t_minus_N & mask) | (t & ~mask)
    for (int i = 0; i < n; ++i) {
        uint32_t v = (t_minus_N[i] & mask) | (t[i] & ~mask);
        r->limbs[i] = v;
    }
    // 상위 워드는 0으로 정리
    for (int i = n; i < MAX_LIMBS; ++i) r->limbs[i] = 0;

    r->size = n;
    bn_normalize(r);
}

// 몽고메리 변환 헬퍼
static void mont_to(Bignum* r, const Bignum* x,
    const Bignum* N, uint32_t n0prime, const Bignum* RR) {
    mont_mul(r, x, RR, N, n0prime);               // r = x * R  (mod N)
}
static void mont_from(Bignum* r, const Bignum* x_bar,
    const Bignum* N, uint32_t n0prime) {
    Bignum one; bn_zero(&one); one.limbs[0] = 1u; one.size = 1;
    mont_mul(r, x_bar, &one, N, n0prime);         // r = x_bar * R^{-1} (mod N)
}

// RR = R^2 mod N (R = 2^(32 * N->size))
// 왼쪽 시프트 대신 "두 배 + 조건부 감산"을 2*n*32 번 반복
static void mont_compute_RR(Bignum* RR, const Bignum* N) {
    bn_zero(RR);
    RR->limbs[0] = 1u;
    RR->size = 1;

    const int nbits = 2 * N->size * (int)LIMB_BITS; // = 2*n*32
    for (int i = 0; i < nbits; ++i) {
        // RR <<= 1
        bn_shift_left(RR, 1);

        // RR = RR mod N (조건부 감산; 최악 1~몇 회)
        if (bn_ucmp(RR, N) >= 0) bn_usub(RR, RR, N);
    }
}

// a <- a mod N (단순 반복 감산; a가 N보다 훨씬 크면 비효율적)
static void bn_reduce_simple(Bignum* a, const Bignum* N) {
    bn_normalize(a);
    while (bn_ucmp(a, N) >= 0) {
        bn_usub(a, a, N);
    }
}

// ========= API 함수 =========
// 모듈러 곱셈 (몽고 메리 곱셈)
// result = a*b mod modulus
void bignum_mod_mul(Bignum* result, const Bignum* a, const Bignum* b, const Bignum* modulus) {
    if (!result || !a || !b || !modulus || modulus->size == 0) {
        bignum_set_zero(result);
        return;
    }
    if ((modulus->limbs[0] & 1u) == 0) { // 몽고메리 전제: N 홀수
        // 필요 시 비-몽고 경로를 따로 구현
        bignum_set_zero(result);
        return;
    }

    // 1) 파라미터
    const uint32_t n0prime = mont_n0prime(modulus->limbs[0]);
    Bignum RR;
    mont_compute_RR(&RR, modulus);

    // 2) 입력을 N으로 감축
    Bignum A;
    bignum_copy(&A, a);
    bn_reduce_simple(&A, modulus);
    Bignum B;
    bignum_copy(&B, b);
    bn_reduce_simple(&B, modulus);

    // 3) 몽고메리 영역으로 진입
    Bignum Abar, Bbar, Zbar;
    mont_to(&Abar, &A, modulus, n0prime, &RR);
    mont_to(&Bbar, &B, modulus, n0prime, &RR);

    // 4) 곱(몽고메리 상태)
    mont_mul(&Zbar, &Abar, &Bbar, modulus, n0prime);

    // 5) 복귀
    mont_from(result, &Zbar, modulus, n0prime);
}

// =====================================================================
//  모듈러 거듭 제곱
// =====================================================================

// ========= 헬퍼 함수 =========


// ========= API 함수 =========
// 모듈러 거듭제곱 result = base^exp mod modulus
// result: 결과를 저장할 Bignum 포인터
// base: 밑
// exp: 지수
// modulus: 모듈러스
void bignum_mod_exp(Bignum* result, const Bignum* base, const Bignum* exp, const Bignum* modulus) {
    // 예외/엣지 처리
    if (!result || !base || !exp || !modulus || modulus->size == 0) {
        bignum_set_zero(result);
        return;
    }
    // 짝수 모듈러스는 몽고메리 불가 (RSA 모듈러스는 홀수이므로 OK)
    if ((modulus->limbs[0] & 1u) == 0u) {
        // 필요하면 여기서 다른(비몽고) 루틴으로 fallback 하세요.
        bignum_set_zero(result);
        return;
    }

    // --- 몽고메리 사전 준비 ---
    const uint32_t n0prime = mont_n0prime(modulus->limbs[0]);

    Bignum RR;
    mont_compute_RR(&RR, modulus);  // RR = R^2 mod N

    // a = base mod N (mont_mul은 a,b < N 가정으로 사용하는 게 안전)
    Bignum a;
    bignum_copy(&a, base);
    bn_reduce_simple(&a, modulus);

    // a_bar = a * R mod N = MonPro(a, RR)
    Bignum a_bar;
    bignum_init(&a_bar);
    mont_mul(&a_bar, &a, &RR, modulus, n0prime);

    // res_bar = 1 * R mod N  (몽고메리의 '1')
    Bignum one;
    bignum_init(&one);
    one.limbs[0] = 1u;
    one.size = 1;
    Bignum res_bar;
    bignum_init(&res_bar);
    mont_mul(&res_bar, &one, &RR, modulus, n0prime);

    // --- 지수승 루프 (LSB-first) ---
    const int bits = bn_bit_length(exp);
    for (int i = 0; i < bits; ++i) {
        if (bn_get_bit(exp, i)) {
            // res_bar = res_bar * a_bar (mod N)
            mont_mul(&res_bar, &res_bar, &a_bar, modulus, n0prime);
        }
        // a_bar = a_bar^2 (mod N)
        mont_mul(&a_bar, &a_bar, &a_bar, modulus, n0prime);
    }

    // --- 몽고메리 영역 → 일반 영역 복귀: res = res_bar * 1 * R^{-1} mod N
    mont_mul(result, &res_bar, &one, modulus, n0prime);
}