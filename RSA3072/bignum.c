#include "bignum.h"
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h> // 추가
#include <complex.h> // 추가

// 카라츠바 알고리즘에서 연산 속도를 높이기 위해 FFT 사용
// FFT에서 회전인자 e^{±iθ} 계산 시 필요한 π
// 일부 구현에서는 <math.h>가 M_PI를 제공하지 않으므로, 없으면 직접 정의
// 접미사 L은 long double 정밀도(cosl/sinl)와 맞추기 위함
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif


#define LIMB_BITS 32u
#define MAX_LIMBS BIGNUM_ARRAY_SIZE

// FFT 전환 임계값(32비트 limb 개수 기준), 환경에 맞게 조정
#define MUL_FFT_THRESHOLD_LIMBS 64

// =====================================================================
//  공통 헬퍼 함수
// =====================================================================

// 비트 길이 
static int bn_bit_length(const Bignum* a) {
    if (a->size == 0) return 0;
    uint32_t ms = a->limbs[a->size - 1];
    int bits = (a->size - 1) * (int)LIMB_BITS;
    int leading = 32;
    while (leading > 0 && ((ms >> (leading - 1)) & 1u) == 0u) leading--;
    return bits + leading;
}

// i번째 비트(0 = LSB) 조회 
// 반환값이 뭔지 설명
static int bn_get_bit(const Bignum* a, int bit_index) {
    if (bit_index < 0) return 0;
    int limb = bit_index / (int)LIMB_BITS;
    int off = bit_index % (int)LIMB_BITS;
    if (limb >= a->size) return 0;
    return (int)((a->limbs[limb] >> off) & 1u);
}

// r += small(0..2^32-1) 설명 다시
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

// r *= small(0..2^32-1) 설명 다시
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
static void bn_shift_left(Bignum* a, int shift) {

}

// n비트 오른쪽 시프트 (a >>= n)
static void bn_shift_right(Bignum* a, int shift) {

}

// =====================================================================
//  큰 수 초기화, 복사, 해제, 변환 등
// =====================================================================

// ========= API 함수 =========
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

// =====================================================================
//  큰 수 비교
// =====================================================================

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
// r = a * b (FFT 사용)

// 상위 0(limb)들 제거 및 유효한 길이를 반환
static inline int bn_limbs_trim(const uint32_t* x, int n) {    // x: limb 배열, n: 배열 길이
    int s = n;                                                 // s에 현재 길이를 복사
    while (s > 0 && x[s - 1] == 0u) --s;                       // 최상위에서부터 0인 limb를 줄여나감
    return s;                                                  // 유효 limb 개수 반환
}

// 반복형 Cooley–Tukey FFT (long double complex 버전)
static void fft_ldc(long double complex* a, int n, int invert) { // a: 복소수 배열, n: 길이 (2의 거듭제곱), invert: 0=정변환, 1=역변환
    // 비트 반전(bit-reversal) 퍼뮤테이션: 나중의 나비연산이 연속 메모리에서 진행되도록 재배열
    for (int i = 1, j = 0; i < n; ++i) {                       // i: 진행 인덱스, j: 비트 반전된 인덱스
        int bit = n >> 1;                                      // 최상위 비트부터 검사
        for (; j & bit; bit >>= 1) j ^= bit;                   // j에서 연속된 1을 0으로 만들며 위쪽 비트 이동
        j ^= bit;                                              // 해당 비트를 1로 세움
        if (i < j) {                                           // i < j이면 두 위치를 교환
            long double complex t = a[i];                      // 임시 저장
            a[i] = a[j];                                       // 스왑
            a[j] = t;                                          // 스왑
        }
    }
    // size=2,4,8,... 단계별 나비 연산
    for (int len = 2; len <= n; len <<= 1) {                   // len: 현재 스테이지 블록 길이
        long double ang = 2.0L * M_PI / (long double)len * (invert ? -1.0L : 1.0L); // 각도 θ = ±2π/len
        long double complex wlen = cosl(ang) + sinl(ang)*I;    // wlen = e^{±iθ} (회전 인자)
        int half = len >> 1;                                   // 블록의 절반 길이
        for (int i = 0; i < n; i += len) {                     // 블록 단위로 순회
            long double complex w = 1.0L + 0.0L*I;             // w = 1 (회전 인자의 누적 곱)
            for (int j = 0; j < half; ++j) {                   // 블록 내에서 나비 연산 수행
                long double complex u = a[i + j];              // u: 상단 값
                long double complex v = a[i + j + half] * w;   // v: 하단 값 × 현재 회전 인자
                a[i + j]        = u + v;                       // 상단 갱신: u+v
                a[i + j + half] = u - v;                       // 하단 갱신: u-v
                w *= wlen;                                     // 다음 위치를 위한 회전 인자 갱신
            }
        }
    }
    if (invert) {                                              // 역변환일 때는 1/n로 나눠서 평균
        for (int i = 0; i < n; ++i) a[i] /= (long double)n;    // 크기 정규화
    }
}

// 32비트 limb 배열 a(길이 an), b(길이 bn)를 곱해 out_limbs/out_size로 반환 (FFT 사용)
static void bn_mul_fft_16bit_u32(                              // FFT 내부 기수는 2^16(=65536) 사용
    uint32_t** out_limbs, int* out_size,                       // 출력: 동적할당된 32비트 limb 배열, 그 길이
    const uint32_t* a, int an,                                 // 입력: 피승수 a, limb 개수 an
    const uint32_t* b, int bn)                                 // 입력: 승수   b, limb 개수 bn
{
    if (an == 0 || bn == 0) {                                  // 둘 중 하나가 0 길이면 결과는 0
        *out_limbs = (uint32_t*)calloc(1, sizeof(uint32_t));   // 길이 0에 맞춰 최소 배열 할당 (여기선 1개 0)
        *out_size  = 0;                                        // 유효 길이는 0
        return;                                                // 종료
    }

    // (1) 32비트 limb → 16비트 digit 분해: 각 limb을 [low16, high16] 두 자리로 쪼갬
    int la = an * 2;                                           // a의 16비트자리 개수
    int lb = bn * 2;                                           // b의 16비트자리 개수

    // (2) FFT 길이 n: la+lb 이상의 2의 거듭제곱으로 맞춤 (컨볼루션 결과 길이를 커버)
    int n = 1;                                                 // n 초기값
    while (n < la + lb) n <<= 1;                               // n을 2배씩 늘려 충분히 크게

    // (3) 복소수 버퍼 할당 및 16비트 digit 로드
    long double complex* fa = (long double complex*)calloc(n, sizeof(long double complex)); // a의 스펙트럼 버퍼
    long double complex* fb = (long double complex*)calloc(n, sizeof(long double complex)); // b의 스펙트럼 버퍼

    for (int i = 0; i < an; ++i) {                             // a의 각 32비트 limb에 대해
        uint32_t x = a[i];                                     // limb 값 x
        fa[2*i + 0] = (long double)(x & 0xFFFFu);              // 하위 16비트 → 16진수 자리 0
        fa[2*i + 1] = (long double)(x >> 16);                  // 상위 16비트 → 16진수 자리 1
    }
    for (int i = 0; i < bn; ++i) {                             // b도 동일하게 분해
        uint32_t x = b[i];                                     // limb 값 x
        fb[2*i + 0] = (long double)(x & 0xFFFFu);              // 하위 16비트
        fb[2*i + 1] = (long double)(x >> 16);                  // 상위 16비트
    }

    // (4) FFT → element-wise multiplication → IFFT 
    fft_ldc(fa, n, 0);                                         // a에 FFT 적용
    fft_ldc(fb, n, 0);                                         // b에 FFT 적용
    for (int i = 0; i < n; ++i) fa[i] *= fb[i];                // 주파수 영역에서 element-wise multiplication
    fft_ldc(fa, n, 1);                                         // IFFT로 시간영역으로 되돌림

    // (5) 반올림 및 캐리 전파 (기수 2^16)
    uint32_t* conv16 = (uint32_t*)calloc(n + 2, sizeof(uint32_t)); // 16비트 자리 결과 저장 버퍼
    long long carry = 0;                                        // 16비트 기수에서의 캐리
    int len16 = n;                                              // 초기 길이는 n (IFFT 결과 길이)
    for (int i = 0; i < n; ++i) {                               // 각 16비트 자리마다
        long long v = llroundl(creall(fa[i])) + carry;          // 실수부를 반올림하여 정수화 + 이전 캐리 합산
        if (v >= 0) {                                           // 비음수면 일반 처리
            conv16[i] = (uint32_t)(v & 0xFFFFLL);               // 하위 16비트 저장
            carry     = v >> 16;                                // 상위 비트는 다음 자리로 캐리
        } else {                                                // 음수가 나올 경우 보정
            long long need = ((-v) + 65535LL) >> 16;            // 올림을 통해 음수 제거에 필요한 단위 수
            v += need * 65536LL;                                // 그만큼 2^16을 더해서 양수화
            carry = -need;                                      // 캐리는 그만큼 음수로 설정
            conv16[i] = (uint32_t)(v & 0xFFFFLL);               // 하위 16비트 저장
        }
    }
    while (carry > 0) {                                         // 마지막까지 남은 캐리를 처리
        conv16[len16++] = (uint32_t)(carry & 0xFFFFLL);         // 16비트 단위로 잘라 저장
        carry >>= 16;                                           // 다음 캐리
    }
    while (len16 > 0 && conv16[len16 - 1] == 0u) --len16;       // 최상위의 불필요한 0 자리 제거

    // (6) 16비트 두 자리 → 32비트 limb로 재조립 (LSB-first)
    int rn = (len16 + 1) / 2;                                   // 16비트 두 자리가 32비트 한 limb
    uint32_t* r = (uint32_t*)calloc(rn, sizeof(uint32_t));      // 결과 32비트 limb 배열 할당
    for (int i = 0; i < rn; ++i) {                              // 각 32비트 limb에 대해
        uint32_t lo = conv16[2*i + 0];                          // 하위 16비트 자리
        uint32_t hi = (2*i + 1 < len16) ? conv16[2*i + 1] : 0u; // 상위 16비트 자리 (없으면 0)
        r[i] = lo | (hi << 16);                                 // hi:lo 를 합쳐 32비트 limb 구성
    }

    free(fa); free(fb); free(conv16);                           // 임시 버퍼 해제

    rn = bn_limbs_trim(r, rn);                                  // 상위 0 limb 제거 후 유효 길이 구함
    *out_limbs = r;                                             // 출력 포인터에 결과 배열 전달
    *out_size  = rn;                                            // 결과 길이 설정
}

// ========= API 함수 =========
// result = a * b
void bignum_multiply(Bignum* result, const Bignum* a, const Bignum* b) { // result: 곱, a/b: 피연산자
    if (a->size == 0 || b->size == 0) {                        // 한쪽이 0이면 곱은 0
        bn_zero(result);                                       // result를 0으로 초기화
        return;                                                // 종료
    }

    uint32_t* tmp = NULL;                                      // FFT 결과(동적 배열) 포인터
    int tmpsize = 0;                                           // FFT 결과의 limb 개수
    bn_mul_fft_16bit_u32(&tmp, &tmpsize, a->limbs, a->size, b->limbs, b->size); // FFT 컨볼루션 수행

    bn_zero(result);                                           // result를 0으로 초기화
    int copy = tmpsize;                                        // 복사할 limb 수 결정
    if (copy > MAX_LIMBS) copy = MAX_LIMBS;                    // result의 고정 용량을 넘지 않도록 제한
    if (copy > 0) {                                            // 복사할 게 있으면
        memcpy(result->limbs, tmp, sizeof(uint32_t) * copy);   // 하위 limb부터 copy개를 복사
        result->size = copy;                                   // 결과 길이 설정
    } else {                                                   // 복사할 게 없다면
        result->size = 0;                                      // 길이 0
    }
    bn_normalize(result);                                      // 혹시 남은 상위 0 limb 정리 (표준형 유지 위해)

    if (tmp) free(tmp);                                        // FFT 임시 결과 메모리 해제
}

// =====================================================================
//  큰 수 나눗셈
// =====================================================================



// =====================================================================
//  큰 수 모듈러 연산
// =====================================================================


