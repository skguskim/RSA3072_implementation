#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> // uint32_t�� ���� ���� �ʺ� ������ ����ϱ� ����

// ū �� ������ ���� �迭 ũ�� (32��Ʈ ���� ����)
#define BIGNUM_ARRAY_SIZE ((5000 / 32) + 1)

typedef struct {
    uint32_t limbs[BIGNUM_ARRAY_SIZE];
    size_t size;
} Bignum; // ū ��

// =============================================================================
// ## 1. ū �� ���� ��� (���: �質��, �迵��, �̳���) ����: bignum.c
// =============================================================================

// Bignum �ʱ�ȭ, ����, ���� �Լ�
void bignum_init(Bignum* bn);
void bignum_set_zero(Bignum* bn);
void bignum_copy(Bignum* dest, const Bignum* src);

// ���ڿ��� Bignum �� ��ȯ �Լ� (16���� ����)
int bignum_from_hex(Bignum* bn, const char* hex_str);
char* bignum_to_hex(const Bignum* bn);

// Bignum ��� �Լ�
// name: ���� ����� ����
void print_bignum(const char* name, const Bignum* bn);

// Bignum �� �Լ� (a > b -> 1, a < b -> -1, a == b -> 0)
int bignum_compare(const Bignum* a, const Bignum* b);

// �⺻ ��Ģ����
void bignum_add(Bignum* result, const Bignum* a, const Bignum* b);
void bignum_subtract(Bignum* result, const Bignum* a, const Bignum* b);
void bignum_multiply(Bignum* result, const Bignum* a, const Bignum* b);
void bignum_divide(Bignum* quotient, Bignum* remainder, const Bignum* a, const Bignum* b);

// ��ⷯ ���� (���� �޸� ����)
void bignum_mod_mul(Bignum* result, const Bignum* a, const Bignum* b, const Bignum* modulus);

// ��ⷯ �ŵ����� (RSA�� �ٽ� ����)
// result = base^exp mod modulus
void bignum_mod_exp(Bignum* result, const Bignum* base, const Bignum* exp, const Bignum* modulus);