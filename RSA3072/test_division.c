#include "bignum.h"
#include <stdio.h>

int test_bignum_divide() {
    Bignum a, b, quotient, remainder;

    // Initialize
    bignum_init(&a);
    bignum_init(&b);
    bignum_init(&quotient);
    bignum_init(&remainder);

    // Test case 1: 100 / 7 = 14 remainder 2
    bignum_from_hex(&a, "64"); // 100 in hex
    bignum_from_hex(&b, "7");  // 7 in hex

    bignum_divide(&quotient, &remainder, &a, &b);

    printf("Test 1: 100 / 7\n");
    print_bignum("Dividend", &a);
    print_bignum("Divisor", &b);
    print_bignum("Quotient", &quotient);
    print_bignum("Remainder", &remainder);
    printf("\n");

    // Test case 2: Large number division
    bignum_from_hex(&a, "123456789ABCDEF0");
    bignum_from_hex(&b, "FEDCBA987654321");

    bignum_divide(&quotient, &remainder, &a, &b);

    printf("Test 2: Large number division\n");
    print_bignum("Dividend", &a);
    print_bignum("Divisor", &b);
    print_bignum("Quotient", &quotient);
    print_bignum("Remainder", &remainder);
    printf("\n");

    // Test case 3: Division by 1
    bignum_from_hex(&a, "DEADBEEF");
    bignum_from_hex(&b, "1");

    bignum_divide(&quotient, &remainder, &a, &b);

    printf("Test 3: Division by 1\n");
    print_bignum("Dividend", &a);
    print_bignum("Divisor", &b);
    print_bignum("Quotient", &quotient);
    print_bignum("Remainder", &remainder);
    printf("-------------------------\n");

    return 0;
}

int test_bignum_modexp() {
    Bignum a, b, result, modulus;

    // Initialize
    bignum_init(&a);
    bignum_init(&b);
    bignum_init(&result);
    bignum_init(&modulus);

    // Test case 1
    bignum_from_hex(&a, "2"); // 100 in hex
    bignum_from_hex(&b, "6");  // 7 in hex
    bignum_from_hex(&modulus, "101");

    bignum_mod_exp(&result, &a, &b, &modulus);

    printf("Test 1: 2^6 mod 100\n");
    print_bignum("base", &a);
    print_bignum("exp", &b);
    print_bignum("mod", &modulus);
    print_bignum("res(40)", &result);
    printf("\n");

    // Test case 2: Large number
    bignum_from_hex(&a, "123456789ABCD");
    bignum_from_hex(&b, "FEDCBA987654321");
    bignum_from_hex(&modulus, "234321");

    bignum_mod_exp(&result, &a, &b, &modulus);

    printf("Test 2: Large number division\n");
    print_bignum("base", &a);
    print_bignum("exp", &b);
    print_bignum("mod", &modulus);
    print_bignum("res(f18ca)", &result);
    printf("\n");

    // Test case 3: exp by 1
    bignum_from_hex(&a, "DEADBEEF");
    bignum_from_hex(&b, "1");
    bignum_from_hex(&modulus, "61");

    bignum_mod_exp(&result, &a, &b, &modulus);

    printf("Test 3: Division by 1\n");
    print_bignum("base", &a);
    print_bignum("exp", &b);
    print_bignum("mod", &modulus);
    print_bignum("res(28)", &result);
    printf("-------------------------\n");

    return 0;
}