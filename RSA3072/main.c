#include "rsa.h"
#include <string.h>
#include <stdio.h>
#include <stdint.h>

// ========== 헥스 문자열 → 바이트 배열 변환 ==========
static int hex_to_bytes(const char *hex, uint8_t *out, int max_len) {
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

static void secure_zero(void* p, size_t n) {
    volatile uint8_t* v = (volatile uint8_t*)p;
    while (n--) *v++ = 0;
}

// ========== DET 벡터 테스트 (RSAES-OAEP with SHA-256, Deterministic Seed) ==========
int test_ent_vector(const char* filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("[-] Failed to open %s\n", filename);
        return 0;
    }
    printf("[*] Testing ENT with %s\n", filename);

    char line[4096], key[64], value[4096];
    Bignum n, e, C_expected, C_actual;
    uint8_t msg_bytes[1024];
    uint8_t em_bytes[384]; // RSA3072_K_BYTES -> 384
    int msg_len = 0;
    Bignum M;

    // Seed를 파일에서 읽는 대신 고정값으로 설정
    const uint8_t fixed_seed[32] = { // OAEP_HLEN -> 32
        0x94,0x11,0x58,0x6e,0x48,0x76,0x6f,0x3d,
        0x56,0x7b,0x14,0x98,0x77,0x05,0x77,0x9a,
        0x32,0x18,0x18,0x8a,0x47,0x47,0xa1,0x01,
        0x0c,0xf0,0x6f,0x56,0x90,0x24,0x18,0x86
    };

    int test_passed = 1, test_count = 0, test_ready = 0;

    bignum_init(&n); bignum_init(&e);
    bignum_init(&C_expected); bignum_init(&C_actual);

    while (fgets(line, sizeof(line), fp)) {
        if (sscanf(line, "%63[^=] = %4095s", key, value) == 2) {
            if (strcmp(key, "n ") == 0) {
                bignum_from_hex(&n, value);
                print_bignum("n", &n);
            } else if (strcmp(key, "e ") == 0) {
                bignum_from_hex(&e, value);
                print_bignum("e", &e);
                printf("\n");
            } else if (strcmp(key, "M ") == 0) {
                bignum_from_hex(&M, value);
            } else if (strcmp(key, "C ") == 0) {
                bignum_init(&C_expected);
                bignum_from_hex(&C_expected, value);
                test_ready = 1;
            }
        }

        if (test_ready) {
            test_count++;
            
            // OAEP 인코딩: 고정된 Seed 값을 사용
            if (rsa_oaep_encode(em_bytes, (const uint8_t*)"", &M, 384, fixed_seed) != 0) { // RSA3072_K_BYTES -> 384
                printf("[-] OAEP encoding failed for test %d\n", test_count);
                test_passed = 0;
                test_ready = 0;
                continue;
            }

            // EM 바이트를 Bignum으로 변환 및 RSA 암호화
            Bignum em_bn;
            bignum_init(&em_bn);
            int rsa_size_bytes = 384; // RSA3072_K_BYTES -> 384
            for (int i = 0; i < rsa_size_bytes; i++) {
                int limb_idx = i / 4;
                int byte_idx = i % 4;
                em_bn.limbs[limb_idx] |= ((uint32_t)em_bytes[i]) << (byte_idx * 8);
            }
            em_bn.size = rsa_size_bytes / 4;
            while (em_bn.size > 0 && em_bn.limbs[em_bn.size - 1] == 0) em_bn.size--;

            RSA_PublicKey pub_key = { n, e };
            rsa_encrypt(&C_actual, &em_bn, &pub_key);

            // 비교
            if (bignum_compare(&C_actual, &C_expected) != 0) {
                printf("[-] DET Test %d failed!\n", test_count);
                test_passed = 0;
                print_bignum("M", &M);
                print_bignum("c_ex", &C_expected);
                print_bignum("c_ac", &C_actual);
                printf("\n");
                break;
            } else {
                printf("[+] DET Test %d passed.\n", test_count);
            }

            secure_zero(em_bytes, sizeof(em_bytes));
            secure_zero(msg_bytes, sizeof(msg_bytes));
            bignum_init(&M);
            bignum_init(&C_expected);
            bignum_init(&C_actual);
            test_ready = 0;
        }
    }
    fclose(fp);

    if (test_passed) printf("[+] All %d DET tests passed!\n", test_count);
    else printf("[-] Some DET tests failed.\n");
    return test_passed;
}


// ========== 메인 함수 ==========
int main() {
    printf("==== RSA-3072 Test Start ====\n");

    // test_bignum_divide();
    test_bignum_modexp();

    int overall_ok = 1;

    if (!test_ent_vector("./test/RSAES_(3072)(65537)(SHA256)_ENT.txt")) overall_ok = 0;
    printf("\n");

    // if (!test_det_vector("./test/RSAES_(3072)(65537)(SHA256)_DET.txt")) overall_ok = 0;
    // printf("\n");
   
    // if (!test_kgt_vector("./test/RSAES_(3072)(65537)(SHA256)_KGT.txt")) overall_ok = 0;

    printf("\n==== RSA-3072 Test Result ====\n");
    if (overall_ok) printf("[+] All tests passed successfully.\n");
    else printf("[-] Some tests failed.\n");

    return overall_ok ? 0 : 1;
}