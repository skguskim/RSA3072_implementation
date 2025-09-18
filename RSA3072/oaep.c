/*
 * RSA-3072 OAEP 패딩 (SHA-256, 테스트용 dummy hash)
 *
 * 고정 크기 버퍼만 사용 (동적 할당 없음, VLA 사용 안 함)
 *
 * 계산 근거(바이트 단위 표기):
 *  - k = 3072 bits = 384 bytes  // 다이어그램의 nLen(비트) ↔ 여기서는 바이트 k
 *  - hLen = 32 bytes
 *  - m_max = k - 2*hLen - 2 = 318 bytes        // 다이어그램의 Len(M) ≤ nLen − 2hLen − 16(bits)
 *  - db_max = k - hLen - 1 = 351 bytes         // 다이어그램의 nLen − hLen − 8(bits)
 *
 * [DIAGRAM GAP] Generate seed / MSB_hLen:
 *  - 다이어그램에는 'Generate seed'와 'MSB_hLen'이 표시되지만,
 *    이 구현은 seed를 내부에서 생성하지 않습니다.
 *    → seed는 호출자가 CSPRNG로 hLen(32B) 난수를 만들어 전달해야 합니다.
 *    → 'MSB_hLen' 단계는 불필요(이미 정확히 hLen 바이트를 받기 때문).
 *
 * [DIAGRAM GAP] Label L:
 *  - 다이어그램의 L(레이블)→Hash(L)은 일반 레이블을 허용합니다.
 *    여기서는 테스트 단순화를 위해 L = ""(빈 레이블)로 고정하여 lHash = H("")만 계산합니다.
 *    필요 시 API를 (const uint8_t* label, size_t label_len)로 확장해 Hash(L)로 교체하세요.
 *
 * 사용법:
 *  - out 버퍼은 반드시 k(=384) 바이트 이상이어야 함.
 *  - seed 길이는 hLen(=32) 바이트여야 함.  // [DIAGRAM] seed는 항상 hLen
 *  - msg_len <= m_max 이어야 함.
 *
 * 주의:
 *  - 이 코드는 테스트용 dummy_sha256 사용. 실제 사용 시에는 SHA256() 등으로 교체.
 *  - seed는 반드시 **CSPRNG**로 매 호출마다 새로 생성해 전달.
 */

#include "rsa.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define RSA3072_K_BYTES 384
#define OAEP_HLEN 32
#define OAEP_DB_MAX (RSA3072_K_BYTES - OAEP_HLEN - 1) /* 351 */
#define OAEP_M_MAX (RSA3072_K_BYTES - 2*OAEP_HLEN - 2) /* 318 */
#define SHA256_DIGEST_LEN OAEP_HLEN

/* 버퍼 상한을 DB 쪽 최대치 기준으로 재정의
 * (다이어그램의 두 번째 MGF 입력: seed = MaskedDB 이므로 길이 = db_len ≤ 351) */
#undef  MGF1_COMBINED_MAX
#define MGF1_COMBINED_MAX (OAEP_DB_MAX + 4) /* 351 + 4 = 355 */

/* --- 테스트용 더미 SHA-256 (실사용 시 실제 SHA256으로 교체) --- */
static void dummy_sha256(const uint8_t *data, size_t len, uint8_t *out32) {
    (void)data; (void)len;
    static const uint8_t dummy_hash[SHA256_DIGEST_LEN] = {
        0x94,0x11,0x58,0x6e,0x48,0x76,0x6f,0x3d,
        0x56,0x7b,0x14,0x98,0x77,0x05,0x77,0x9a,
        0x32,0x18,0x18,0x8a,0x47,0x47,0xa1,0x01,
        0x0c,0xf0,0x6f,0x56,0x90,0x24,0x18,0x86
    };
    memcpy(out32, dummy_hash, SHA256_DIGEST_LEN);
}

/* --- MGF1 (SHA-256 기반) ---
 * [DIAGRAM] MGF 블록: 입력 seed 길이는 32B(hLen) 또는 db_len(351B) 등 가변.
 * - 1차 호출: seed = seed(hLen), mask_len = db_len
 * - 2차 호출: seed = maskedDB(db_len), mask_len = hLen
 */
static int rsa_mgf1_fixed(uint8_t *mask, size_t mask_len, const uint8_t *seed, size_t seed_len) {
    if (!mask || !seed) return 0;
    if (seed_len == 0 || seed_len > OAEP_DB_MAX) return 0; /* 안전 상한 */

    uint8_t digest[SHA256_DIGEST_LEN];
    uint8_t combined[MGF1_COMBINED_MAX]; /* seed(≤351) + counter(4) */
    size_t generated = 0;
    uint32_t counter = 0;

    while (generated < mask_len) {
        /* [DIAGRAM] combined = seed || COUNTER (big-endian) */
        memcpy(combined, seed, seed_len);
        combined[seed_len + 0] = (uint8_t)((counter >> 24) & 0xFF);
        combined[seed_len + 1] = (uint8_t)((counter >> 16) & 0xFF);
        combined[seed_len + 2] = (uint8_t)((counter >> 8) & 0xFF);
        combined[seed_len + 3] = (uint8_t)(counter & 0xFF);

        dummy_sha256(combined, seed_len + 4, digest);

        size_t to_copy = (mask_len - generated < SHA256_DIGEST_LEN)
                           ? (mask_len - generated)
                           : SHA256_DIGEST_LEN;
        memcpy(mask + generated, digest, to_copy);
        generated += to_copy;
        counter++;
    }
    return 1;
}

static void secure_zero(void *p, size_t n) {
    volatile uint8_t *v = (volatile uint8_t*)p;
    while (n--) *v++ = 0;
}

/* --- OAEP 인코딩 (고정 버퍼 사용, 동적할당 없음)
 * [DIAGRAM] 이 함수는 OAEP_Encode에서 EM까지를 생성합니다.
 *  - RSA 모듈러 지수 연산(EM^e mod n)은 포함되지 않습니다(다이어그램 밖).
 *
 * [DIAGRAM GAP] Label L:
 *  - 현재 구현은 L = "" 고정. lHash = H("").
 *  - 필요한 경우 (label, label_len) 인자를 추가해 Hash(L)로 변경하세요.
 */
int rsa_oaep_encode(uint8_t *EM, const uint8_t *L, const Bignum *M, const size_t k, const uint8_t *seed) {
    if (!EM || !L || !M || !seed) return -1;
    if (k != RSA3072_K_BYTES) return -1; /* RSA-3072 전용 */
    if (M->size > OAEP_M_MAX) return -1;
    /* [DIAGRAM] seed 길이는 hLen이어야 함(호출자 책임, Generate seed는 외부에서 수행) */
      /* NOTE: 이 함수는 seed를 내부에서 생성하지 않음.
     *       → 호출자가 CSPRNG로 hLen(32B) 난수를 만들어 전달해야 함.
     *       → 다이어그램의 MSB_hLen 절차는 불필요(이미 정확히 32B를 받음). */

    /* lHash = Hash(L), 여기서는 L = "" (빈 레이블) */
    uint8_t lHash[OAEP_HLEN];
    dummy_sha256((const uint8_t *)"", 0, lHash);

    /* [DIAGRAM] DB = lHash || 0^PS || 0x01 || M
     * ps_len = k - msg_len - 2*hLen - 2              // 다이어그램의 … − 16(bits) 대응
     * db_len = hLen + ps_len + 1 + msg_len = k - hLen - 1
     */
    size_t ps_len = k - M->size - 2*OAEP_HLEN - 2;
    size_t db_len = OAEP_HLEN + ps_len + 1 + M->size;
    if (db_len > OAEP_DB_MAX) return -1;

    /* 내부 버퍼들 (DB 길이 = k - hLen - 1) */
    uint8_t DB[OAEP_DB_MAX];
    uint8_t dbMask[OAEP_DB_MAX];
    uint8_t maskedDB[OAEP_DB_MAX];
    uint8_t seedMask[OAEP_HLEN];
    uint8_t maskedSeed[OAEP_HLEN];

    /* DB 채우기 */
    uint8_t *p = DB;
    memcpy(p, lHash, OAEP_HLEN); p += OAEP_HLEN;
    if (ps_len) memset(p, 0x00, ps_len); p += ps_len;
    *p++ = 0x01; /* 다이어그램의 0^7 | 1 한 바이트(0x01) */
    memcpy(p, M->limbs, M->size);

    /* [DIAGRAM] dbMask = MGF(seed, nLen − hLen − 8) ≡ MGF(seed, db_len bytes) */
    if (!rsa_mgf1_fixed(dbMask, db_len, seed, OAEP_HLEN)) return -1;

    /* maskedDB = DB ⊕ dbMask */
    for (size_t i = 0; i < db_len; ++i) maskedDB[i] = DB[i] ^ dbMask[i];

    /* [DIAGRAM] seedMask = MGF(maskedDB, hLen) */
    if (!rsa_mgf1_fixed(seedMask, OAEP_HLEN, maskedDB, db_len)) return -1;

    /* maskedSeed = seed ⊕ seedMask */
    for (size_t i = 0; i < OAEP_HLEN; ++i) maskedSeed[i] = seed[i] ^ seedMask[i];

    /* [DIAGRAM] EM = 0^8 || maskedSeed || maskedDB (총 길이 k 바이트) */
    EM[0] = 0x00;                      /* 0^8 */
    memcpy(EM + 1, maskedSeed, OAEP_HLEN);
    memcpy(EM + 1 + OAEP_HLEN, maskedDB, db_len);

    /* 민감 중간값 소거 */
    secure_zero(DB, db_len);
    secure_zero(dbMask, db_len);
    secure_zero(maskedDB, db_len);
    secure_zero(seedMask, OAEP_HLEN);
    secure_zero(maskedSeed, OAEP_HLEN);
    secure_zero(lHash, OAEP_HLEN);

    return 0;
}

<<<<<<< HEAD
// OAEP_DECODE
int rsa_oaep_decode(uint8_t* M, size_t* mLen, const size_t mMax, const uint8_t* EM, const size_t k, const uint8_t* L, const size_t Llen) {
    if (!M || !mLen || !EM) return -1;
    const size_t db_len = k - OAEP_HLEN - 1;                   /* nLen − hLen − 8 (bits) ↔ k − hLen − 1 (bytes) */

    /* [Y | maskedSeed | maskedDB] 파싱 */
    uint8_t Y = EM[0];
    const uint8_t* maskedSeed = EM + 1;
    const uint8_t* maskedDB = EM + 1 + OAEP_HLEN;

    /* --- 버퍼 준비 --- */
    uint8_t seed[OAEP_HLEN];
    uint8_t seedMask[OAEP_HLEN];
    uint8_t dbMask[OAEP_DB_MAX];        /* = db_len */
    uint8_t DB[OAEP_DB_MAX];            /* 복원된 DB */
    uint8_t lHash[OAEP_HLEN];
    int bad = 0;                         /* 형식 오류 플래그 */

    /* Y == 0x00 ? */
    bad |= (Y != 0);

    /* lHash = Hash(L) */
    if (L && Llen > 0) dummy_sha256(L, Llen, lHash);
    else               dummy_sha256((const uint8_t*)"", 0, lHash);

    /* seedMask = MGF(maskedDB, hLen) */
    if (!rsa_mgf1_fixed(seedMask, OAEP_HLEN, maskedDB, db_len)) { bad = 1; goto finish; }

    /* seed = maskedSeed XOR seedMask */
    for (size_t i = 0; i < OAEP_HLEN; ++i) seed[i] = maskedSeed[i] ^ seedMask[i];

    /* dbMask = MGF(seed, db_len) */
    if (!rsa_mgf1_fixed(dbMask, db_len, seed, OAEP_HLEN)) { bad = 1; goto finish; }

    /* DB = maskedDB XOR dbMask */
    for (size_t i = 0; i < db_len; ++i) DB[i] = maskedDB[i] ^ dbMask[i];

    /* DB = lHash || PS || 0x01 || M  */

    /* 1) lHash 비교  */
    {
        uint8_t diff = 0;
        for (size_t i = 0; i < OAEP_HLEN; ++i) diff |= (uint8_t)(DB[i] ^ lHash[i]);
        bad |= (diff != 0);
    }

    // 2) PS는 전부 0x00 이어야 하고, 그 다음 바이트가 0x01 이어야 함.
    size_t i = OAEP_HLEN;
    int seen_one = 0;
    size_t one_pos = 0;

    while (i < db_len) {
        uint8_t b = DB[i];
        if (!seen_one) {
            if (b == 0x00) {
                /* PS 영역 계속 */
            }
            else if (b == 0x01) {
                seen_one = 1;
                one_pos = i;
                break;          /* 0x01 발견 → 이후는 메시지 */
            }
            else {
                /* PS 구간에 0x00/0x01 이외의 값 */
                bad = 1;
                break;
            }
        }
        ++i;
    }
    if (!seen_one) bad = 1;            /* 0x01 구분자를 못 찾음 */

    /* 3) 메시지 추출 */
    size_t m_off = seen_one ? (one_pos + 1) : db_len;
    size_t m_len = (m_off <= db_len) ? (db_len - m_off) : 0;

    /* 출력 버퍼 체크 */
    if (m_len > mMax) { bad = 1; goto finish; }

    /* 메시지 복사 */
    if (m_len && m_off < db_len) memcpy(M, DB + m_off, m_len);
    *mLen = m_len;

finish:
    /* 민감 데이터 소거 */
    secure_zero(seed, sizeof(seed));
    secure_zero(seedMask, sizeof(seedMask));
    secure_zero(dbMask, db_len);
    secure_zero(DB, db_len);
    secure_zero(lHash, sizeof(lHash));

    return bad ? -1 : 0;
}
=======
/*
 * rsa_oaep_decode - OAEP 디코딩 (RSA-3072, SHA-256, fixed buffers)
 *
 * 파라미터:
 *  - M_out : 복원된 메시지를 담을 Bignum (호출자는 M_out->limbs 버퍼가 OAEP_M_MAX 이상임을 보장)
 *  - L     : 레이블 (현재 구현은 L=""를 가정해 lHash = H("") 처리)
 *  - EM    : 인코딩된 메시지(길이 k 바이트, 보통 RSA 연산 후 얻은 바이트 스트림)
 *  - k     : EM 길이(바이트) — 이 구현은 RSA3072_K_BYTES(384)만 허용
 *
 * 반환:
 *   0  : 성공 (M_out->limbs에 복원된 메시지, M_out->size에 길이 저장)
 *  -1  : 실패 (포맷/검증 오류)
 */
int rsa_oaep_decode(Bignum *M_out, const uint8_t *L, const uint8_t *EM, size_t k) {
    if (!M_out || !L || !EM) return -1;
    if (k != RSA3072_K_BYTES) return -1;

    /* EM 구조: 0x00 || maskedSeed(hLen) || maskedDB(db_len) */
    if (EM[0] != 0x00) return -1;

    size_t db_len = k - OAEP_HLEN - 1;
    if (db_len > OAEP_DB_MAX) return -1;

    /* 고정 내부 버퍼 */
    uint8_t DB[OAEP_DB_MAX];
    uint8_t dbMask[OAEP_DB_MAX];
    uint8_t maskedDB[OAEP_DB_MAX];
    uint8_t seedMask[OAEP_HLEN];
    uint8_t maskedSeed[OAEP_HLEN];
    uint8_t seed[OAEP_HLEN];
    uint8_t lHash[OAEP_HLEN];

    /* 분리 */
    memcpy(maskedSeed, EM + 1, OAEP_HLEN);
    memcpy(maskedDB,   EM + 1 + OAEP_HLEN, db_len);

    /* seedMask = MGF(maskedDB, hLen) */
    if (!rsa_mgf1_fixed(seedMask, OAEP_HLEN, maskedDB, db_len)) {
        secure_zero(maskedDB, db_len);
        secure_zero(maskedSeed, OAEP_HLEN);
        secure_zero(seedMask, OAEP_HLEN);
        return -1;
    }

    /* seed = maskedSeed ⊕ seedMask */
    for (size_t i = 0; i < OAEP_HLEN; ++i) seed[i] = maskedSeed[i] ^ seedMask[i];

    /* dbMask = MGF(seed, db_len) */
    if (!rsa_mgf1_fixed(dbMask, db_len, seed, OAEP_HLEN)) {
        secure_zero(DB, db_len);
        secure_zero(dbMask, db_len);
        secure_zero(maskedDB, db_len);
        secure_zero(seedMask, OAEP_HLEN);
        secure_zero(maskedSeed, OAEP_HLEN);
        secure_zero(seed, OAEP_HLEN);
        return -1;
    }

    /* DB = maskedDB ⊕ dbMask */
    for (size_t i = 0; i < db_len; ++i) DB[i] = maskedDB[i] ^ dbMask[i];

    /* lHash = Hash(L) -- 현재 구현은 L=""를 가정(빈 레이블) */
    dummy_sha256((const uint8_t *)"", 0, lHash);

    /* lHash 비교 */
    if (memcmp(DB, lHash, OAEP_HLEN) != 0) {
        secure_zero(DB, db_len);
        secure_zero(dbMask, db_len);
        secure_zero(maskedDB, db_len);
        secure_zero(seedMask, OAEP_HLEN);
        secure_zero(maskedSeed, OAEP_HLEN);
        secure_zero(seed, OAEP_HLEN);
        secure_zero(lHash, OAEP_HLEN);
        return -1;
    }

    /* 0x01 구분자 찾기 (lHash 뒤의 PS(0x00...0x00) 다음에 0x01이 와야 함) */
    size_t idx = OAEP_HLEN;
    while (idx < db_len && DB[idx] == 0x00) idx++;
    if (idx >= db_len || DB[idx] != 0x01) {
        secure_zero(DB, db_len);
        secure_zero(dbMask, db_len);
        secure_zero(maskedDB, db_len);
        secure_zero(seedMask, OAEP_HLEN);
        secure_zero(maskedSeed, OAEP_HLEN);
        secure_zero(seed, OAEP_HLEN);
        secure_zero(lHash, OAEP_HLEN);
        return -1;
    }

    /* 메시지 위치 및 길이 계산 */
    size_t msg_offset = idx + 1;
    size_t msg_len = db_len - msg_offset;
    if (msg_len > OAEP_M_MAX) {
        secure_zero(DB, db_len);
        secure_zero(dbMask, db_len);
        secure_zero(maskedDB, db_len);
        secure_zero(seedMask, OAEP_HLEN);
        secure_zero(maskedSeed, OAEP_HLEN);
        secure_zero(seed, OAEP_HLEN);
        secure_zero(lHash, OAEP_HLEN);
        return -1;
    }

    /* 호출자 Bignum 버퍼로 복사 (호출자는 M_out->limbs의 용량을 보장해야 함) */
    memcpy(M_out->limbs, DB + msg_offset, msg_len);
    M_out->size = msg_len;

    /* 민감값 소거 */
    secure_zero(DB, db_len);
    secure_zero(dbMask, db_len);
    secure_zero(maskedDB, db_len);
    secure_zero(seedMask, OAEP_HLEN);
    secure_zero(maskedSeed, OAEP_HLEN);
    secure_zero(seed, OAEP_HLEN);
    secure_zero(lHash, OAEP_HLEN);

    return 0;
}
>>>>>>> e815dc3fc9f68d3747a3224fff4947363c3fbc14
