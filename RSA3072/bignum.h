#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> // uint32_t와 같은 고정 너비 정수를 사용하기 위함

// 큰 수 저장을 위한 배열 크기 (32비트 정수 기준)
#define BIGNUM_ARRAY_SIZE ((5000 / 32) + 1)

typedef struct {
    uint32_t limbs[BIGNUM_ARRAY_SIZE];
    size_t size;
} Bignum; // 큰 수

// =============================================================================
// ## 1. 큰 수 연산 모듈 (담당: 김나현, 김영진, 이나원) 파일: bignum.c
// =============================================================================

// Bignum 초기화, 복사, 해제 함수
void bignum_init(Bignum* bn);
void bignum_set_zero(Bignum* bn);
void bignum_copy(Bignum* dest, const Bignum* src);

// 문자열과 Bignum 간 변환 함수 (16진수 권장)
int bignum_from_hex(Bignum* bn, const char* hex_str);
char* bignum_to_hex(const Bignum* bn);

// Bignum 출력 함수
// name: 같이 출력할 설명
void print_bignum(const char* name, const Bignum* bn);

// Bignum 비교 함수 (a > b -> 1, a < b -> -1, a == b -> 0)
int bignum_compare(const Bignum* a, const Bignum* b);

// 기본 사칙연산
void bignum_add(Bignum* result, const Bignum* a, const Bignum* b);
void bignum_subtract(Bignum* result, const Bignum* a, const Bignum* b);
void bignum_multiply(Bignum* result, const Bignum* a, const Bignum* b);
void bignum_divide(Bignum* quotient, Bignum* remainder, const Bignum* a, const Bignum* b);

// 모듈러 곱셈 (몽고 메리 곱셈)
void bignum_mod_mul(Bignum* result, const Bignum* a, const Bignum* b, const Bignum* modulus);

// 모듈러 거듭제곱 (RSA의 핵심 연산)
// result = base^exp mod modulus
void bignum_mod_exp(Bignum* result, const Bignum* base, const Bignum* exp, const Bignum* modulus);