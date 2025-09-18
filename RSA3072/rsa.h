#pragma once

#include "bignum.h"

// RSA 키 비트 길이 정의 (3072 비트)
#define RSA_KEY_BITS 3072
// 소수 p, q의 비트 길이 (키 길이의 절반)
#define RSA_PRIME_BITS 128 // (RSA_KEY_BITS / 2)
// 밀러-라빈 테스트 반복 횟수
#define MILLER_RABIN_ROUNDS 64

typedef struct {
    Bignum n; // Modulus
    Bignum e; // Public Exponent
} RSA_PublicKey; // 공개키

typedef struct {
    Bignum n; // Modulus
    Bignum d; // Private Exponent
    // CRT를 위해 필요한 값
    Bignum p;
    Bignum q;
    Bignum dP; // d mod (p-1)
    Bignum dQ; // d mod (q-1)
    Bignum qInv; // q^(-1) mod p
} RSA_PrivateKey; // 개인키


// =============================================================================
// ## 2. 안전한 난수 생성기 (담당: 정태진, 김주영) 파일: random.c
// =============================================================================

/**
 * 암호학적으로 안전한 난수를 생성하는 함수
 * buffer: 난수를 저장할 버퍼
 * size: 생성할 바이트 수
 */
int generate_secure_random(unsigned char* buffer, size_t size);

// =============================================================================
// ## 3. 소수 판별법 (담당: 김성우, 김민수) 파일: prime.c
// =============================================================================

/**
 * 밀러-라빈 소수 판별법
 * n: 판별할 큰 수 (소수 후보)
 * k: 테스트 반복 횟수
 * return 소수일 확률이 높으면 1, 아니면 0
 */
int is_probably_prime(const Bignum* n, int k);

// =============================================================================
// ## 4. 소수 생성 모듈 (담당: 이나원) 파일: prime.c
// =============================================================================

/**
 * 지정된 비트 길이의 소수를 생성하는 함수
 * prime: 생성된 소수를 저장할 Bignum 포인터
 * bits: 원하는 소수의 비트 길이 (지금 과제는 3072)
 */
void generate_prime(Bignum* prime, int bits);


// =============================================================================
// ## 5. 파라미터 계산 (담당: 김민수) 파일: rsa.c
// =============================================================================

/**
 * p와 q로부터 RSA 키 쌍(공개키, 개인키)을 생성
 * pub_key: 생성된 공개키를 저장할 구조체 포인터
 * priv_key: 생성된 개인키를 저장할 구조체 포인터
 */
void rsa_generate_keys(RSA_PublicKey* pub_key, RSA_PrivateKey* priv_key, const Bignum* p, const Bignum* q);

// =============================================================================
// ## 6. 암호화 및 복호화 (담당: 김영진) 파일: rsa.c
// =============================================================================

/**
 * RSA 암호화
 * ciphertext: 암호화된 결과를 저장할 Bignum 포인터
 * message: 암호화할 원문 Bignum
 * pub_key: 공개키
 */
int rsa_encrypt(uint8_t* ciphertext, const Bignum* message, const RSA_PublicKey* pub_key);

/**
 * RSA 복호화
 * message: 복호화된 결과를 저장할 Bignum 포인터
 * ciphertext: 복호화할 암호문 Bignum
 * priv_key: 개인키
 */
int rsa_decrypt(Bignum* message, const Bignum* ciphertext, const RSA_PrivateKey* priv_key);

// =============================================================================
// ## 7. 테스트 벡터 (담당: 김강민) 파일: main.c
// =============================================================================

/*
* 테스트 벡터 검증 -> main.c, 하단의 함수는 RSA-OAEP 패딩을 위한 구현체 
*/


// 원본 메시지를 OAEP 규칙에 맞게 패딩하여 암호화에 사용할 메시지 블록(EM)을 생성

// out: 패딩된 메시지(EM)를 저장할 버퍼
// msg: 원본 메시지
// msg_len: 원본 메시지 길이 (바이트)
// k: RSA 모듈러스 길이 (바이트)
// seed: OAEP에 사용될 시드(SHA256 해시 길이인 32바이트여야 함)

int rsa_oaep_encode(uint8_t *out, const uint8_t *L, const Bignum *M, const size_t k, const uint8_t *seed);
int rsa_oaep_decode(uint8_t* M, size_t* mLen, const size_t mMax, const uint8_t* EM, const size_t k, const uint8_t* L, const size_t Llen);

