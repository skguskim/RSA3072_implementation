#include "RSA3072/bignum.h"
#include <stdio.h>

int main() {
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

    return 0;
}