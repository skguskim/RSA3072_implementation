#include "rsa.h"
#include <string.h>

// ========== 헥스 문자열 → 바이트 배열 변환 ==========
static int hex_to_bytes(const char *hex, unsigned char *out, int max_len) {
    if (!hex) return 0;
    int len = strlen(hex);
    if (len % 2 != 0) {
        fprintf(stderr, "Invalid hex string length: %d\n", len);
        return -1;
    }
    int out_len = len / 2;
    if (out_len > max_len) {
        fprintf(stderr, "Output buffer too small: %d > %d\n", out_len, max_len);
        return -1;
    }
    for (int i = 0; i < out_len; i++) {
        if (sscanf(hex + 2*i, "%2hhx", &out[i]) != 1) {
            fprintf(stderr, "Failed to parse hex at index %d\n", i);
            return -1;
        }
    }
    return out_len;
}

// ========== ENT 벡터 테스트 ==========
int test_ent_vector(const char* filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("[-] Failed to open %s\n", filename);
        return 0;
    }
    printf("[*] Testing ENT with %s\n", filename);

    char line[4096], key[64], value[4096];
    Bignum n, e, m, c_expected, c_actual;
    bignum_init(&n); bignum_init(&e); bignum_init(&m);
    bignum_init(&c_expected); bignum_init(&c_actual);

    int test_passed = 1, test_count = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (sscanf(line, "%63[^=] = %4095s", key, value) == 2) {
            if (strcmp(key, "n ") == 0) {
                bignum_from_hex(&n, value);
            } else if (strcmp(key, "e ") == 0) {
                bignum_from_hex(&e, value);
            } else if (strcmp(key, "M ") == 0) {
                bignum_from_hex(&m, value);
            } else if (strcmp(key, "C ") == 0) {
                bignum_from_hex(&c_expected, value);

                RSA_PublicKey pub = { n, e };
                rsa_encrypt(&c_actual, &m, &pub);

                test_count++;
                if (bignum_compare(&c_actual, &c_expected) != 0) {
                    printf("[-] ENT Test %d failed!\n", test_count);
                    test_passed = 0;
                    print_bignum("n", &n);
                    print_bignum("e", &e);
                    print_bignum("M", &m);
                    print_bignum("c_ex", &c_expected);
                    print_bignum("c_ac", &c_actual);
                }
            }
        }
    }
    fclose(fp);

    if (test_passed) printf("[+] All %d ENT tests passed!\n", test_count);
    else printf("[-] Some ENT tests failed.\n");
    return test_passed;
}

// ========== DET 벡터 테스트 ==========
int test_det_vector(const char* filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("[-] Failed to open %s\n", filename);
        return 0;
    }
    printf("[*] Testing DET with %s\n", filename);

    char line[4096], key[64], value[4096];
    Bignum n, e, C, C_actual;
    unsigned char msg_bytes[1024], seed_bytes[32], em_bytes[RSA_KEY_BITS/8];
    Bignum M;
    uint8_t seed = 0;
    bignum_init(&M);
    int test_passed = 1, test_count = 0, test_ready = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (sscanf(line, "%63[^=] = %4095s", key, value) == 2) {
            if (strcmp(key, "n ") == 0) {
                bignum_init(&n); bignum_from_hex(&n, value);
            } else if (strcmp(key, "e ") == 0) {
                bignum_init(&e); bignum_from_hex(&e, value);
            } else if (strcmp(key, "M ") == 0) {
                bignum_from_hex(&M, value);
            } else if (strcmp(key, "Seed ") == 0) {
                seed = bignum_from_hex(value, seed_bytes, sizeof(seed_bytes));
            } else if (strcmp(key, "C ") == 0) {
                bignum_init(&C); bignum_from_hex(&C, value);
                test_ready = 1;
            }
        }
        if (test_ready) {
            test_count++;

            int rsa_size_bytes = RSA_KEY_BITS / 8;
            rsa_oaep_encode(em_bytes, "", &M, RSA_KEY_BITS, seed);

            Bignum em_bn;
            bignum_init(&em_bn);
            for (int i = 0; i < rsa_size_bytes; i++) {
                int limb_idx = i / 4;
                int byte_idx = i % 4;
                em_bn.limbs[limb_idx] |= ((uint32_t)em_bytes[i]) << (byte_idx * 8);
            }
            em_bn.size = rsa_size_bytes / 4;
            while (em_bn.size > 0 && em_bn.limbs[em_bn.size - 1] == 0) em_bn.size--;

            RSA_PublicKey pub_key = { n, e };
            rsa_encrypt(&C_actual, &em_bn, &pub_key);

            if (bignum_compare(&C_actual, &C) != 0) {
                printf("[-] DET Test %d failed!\n", test_count);
                test_passed = 0;
            }
            test_ready = 0;
        }
    }
    fclose(fp);

    if (test_passed) printf("[+] All %d DET tests passed!\n", test_count);
    else printf("[-] Some DET tests failed.\n");
    return test_passed;
}

// ========== KGT 벡터 테스트 ==========
int test_kgt_vector(const char* filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("[-] Failed to open %s\n", filename);
        return 0;
    }
    printf("[*] Testing KGT with %s\n", filename);

    char line[4096], key[64], value[4096];
    Bignum p_file, q_file, n_file, e_file, d_file, dP_file, dQ_file, qInv_file;
    bignum_init(&p_file); bignum_init(&q_file); bignum_init(&n_file); bignum_init(&e_file);
    bignum_init(&d_file); bignum_init(&dP_file); bignum_init(&dQ_file); bignum_init(&qInv_file);

    int test_count = 0, test_ready = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (sscanf(line, "%63[^=] = %4095s", key, value) == 2) {
            if (strcmp(key, "n ") == 0) bignum_from_hex(&n_file, value);
            else if (strcmp(key, "e ") == 0) bignum_from_hex(&e_file, value);
            else if (strcmp(key, "p ") == 0) bignum_from_hex(&p_file, value);
            else if (strcmp(key, "q ") == 0) bignum_from_hex(&q_file, value);
            else if (strcmp(key, "d ") == 0) bignum_from_hex(&d_file, value);
            else if (strcmp(key, "dP ") == 0) bignum_from_hex(&dP_file, value);
            else if (strcmp(key, "dQ ") == 0) bignum_from_hex(&dQ_file, value);
            else if (strcmp(key, "qInv ") == 0) {
                bignum_from_hex(&qInv_file, value);
                test_ready = 1;
            }
        }

        if (test_ready) {
            test_count++;
            RSA_PublicKey pub_gen;
            RSA_PrivateKey priv_gen;
            rsa_generate_keys(&pub_gen, &priv_gen, &p_file, &q_file);

            if (bignum_compare(&pub_gen.n, &n_file) != 0 ||
                bignum_compare(&pub_gen.e, &e_file) != 0 ||
                bignum_compare(&priv_gen.d, &d_file) != 0 ||
                bignum_compare(&priv_gen.p, &p_file) != 0 ||
                bignum_compare(&priv_gen.q, &q_file) != 0 ||
                bignum_compare(&priv_gen.dP, &dP_file) != 0 ||
                bignum_compare(&priv_gen.dQ, &dQ_file) != 0 ||
                bignum_compare(&priv_gen.qInv, &qInv_file) != 0) {
                printf("[-] KGT Test %d failed!\n", test_count);
                fclose(fp);
                return 0;
            } else {
                printf("[+] KGT Test %d passed.\n", test_count);
            }
            test_ready = 0;
        }
    }
    fclose(fp);
    printf("[+] All %d KGT tests completed successfully.\n", test_count);
    return 1;
}

// ========== 메인 ==========
int main() {
    printf("==== RSA-3072 Test Start ====\n");

    // test_bignum_divide();
    // test_bignum_modexp();

    int overall_ok = 1;

    if (!test_ent_vector("./test/RSAES_(3072)(65537)(SHA256)_ENT.txt")) overall_ok = 0;
    printf("\n");
    if (!test_det_vector("./test/RSAES_(3072)(65537)(SHA256)_DET.txt")) overall_ok = 0;
    printf("\n");
    if (!test_kgt_vector("./test/RSAES_(3072)(65537)(SHA256)_KGT.txt")) overall_ok = 0;

    printf("\n==== RSA-3072 Test Result ====\n");
    if (overall_ok) printf("[+] All tests passed successfully.\n");
    else printf("[-] Some tests failed.\n");

    return overall_ok ? 0 : 1;
}
