#include "rsa.h"
#include <stdbool.h>

int is_probably_prime(const Bignum* n, const int k) {//k=64
    // 변수 선언 및 초기화
    bool isprime=true;
    Bignum q, r, n_minus_one, result, n_minus_two, d, result_sq, result_sq_quotient, sudo_number;
    int s=0;
    const int size=sizeof(n); // 난수열 크기 지정
    const Bignum zero={.limbs={0}, .size=0};
    const Bignum one={.limbs={1}, .size=1};
    const Bignum two={.limbs={2}, .size=1};
    const Bignum three={.limbs={3}, .size=1};
    bignum_init(&q);
    bignum_init(&r);
    bignum_init(&n_minus_one);
    bignum_init(&result);
    bignum_init(&n_minus_two);
    bignum_init(&d);
    bignum_init(&result_sq);
    bignum_init(&result_sq_quotient);
    bignum_init(&sudo_number);

    // test_value = n-1;
    bignum_subtract(&n_minus_one,n,&one);
    const Bignum test_value=n_minus_one;

    // range = n-3;
    bignum_subtract(&n_minus_two, n, &two);
    const Bignum range=n_minus_two;

    // 2의 배수이면 바로 0 반환
    bignum_divide(&q, &r, n, &two);
    if(bignum_compare(&r, &zero)==0) return 0;
    if(bignum_compare(n, &one)==0) return 0;

    //2. n-1 = 2^k * d (d는 홀수)
    bignum_copy(&d, &test_value);
    while(1) {
        bignum_divide(&q, &r, &d, &two);
        if (bignum_compare(&r, &zero) != 0) break;
        bignum_copy(&d, &q);  
        s++;
    }

    for(int i=0; i<k; i++) {
        // 2이상 n-2이하 난수 뽑을때까지 반복문
        while(1){
            bignum_set_zero(&sudo_number);
            generate_secure_random(&sudo_number, size);
            if(bignum_compare(&sudo_number, &two)!=-1 && bignum_compare(&sudo_number,&range)!=1) {
                break;
            }
        }
        
        // 3. a^d mod n = 1 or n-1 n은 소수
        bignum_mod_exp(&result, &sudo_number, &d, n);
        if(bignum_compare(&result, &one)==0 || bignum_compare(&result, &test_value)==0) {
            continue;
        }

        // 4. r을 0~k-1까지 증가 ((a^d)^2) mod n이 n-1 인지 검사
        else {
            isprime=false;
            for (int m = 1; m < s; m++) {
                bignum_multiply(&result_sq, &result, &result);
                bignum_divide(&result_sq_quotient, &result, &result_sq, n);
                // n-1이면 계속
                if (bignum_compare(&result, &test_value) == 0) {
                    isprime = true;
                    break;
                }
                // 1이면 합성수
                if(bignum_compare(&result, &one)==0) {
                    return 0;
                } 
            }
            if (!isprime) {
                return 0; // 합성수
            }
        }
        //5. 만약 n-1이 한번이라도 나오면 n은 소수
        //6. 위의 모든 조건을 통과하지 못하면 n은 합성수
    }
    return 1;
}



void generate_prime(Bignum* prime, int bits) {
    Bignum candidate;

    while (1) {
        bignum_init(&candidate);

        // 1. 난수 생성에 필요한 바이트 수 계산
        int bytes = (bits + 7) / 8;
        unsigned char* buffer = (unsigned char*)malloc(bytes);
        if (buffer == NULL) continue;

        // 2. 안전한 난수 생성
        if (generate_secure_random(buffer, bytes) != 0) {
            free(buffer);
            continue;
        }

        // 3. 바이트를 Bignum으로 변환
        for (int i = 0; i < bytes; i++) {
            int limb_idx = i / 4; // 몇 번째 limb인지
            int byte_idx = i % 4; //그 limb 안에서 몇 번째 바이트 자리인지

            candidate.limbs[limb_idx] |= ((uint32_t)buffer[i]) << (byte_idx * 8);
            if (limb_idx >= candidate.size) {
                candidate.size = limb_idx + 1;
            }
        }

        free(buffer);

        // 4. 밀러-라빈 테스트로 소수인지 확인
        if (is_probably_prime(&candidate, MILLER_RABIN_ROUNDS)) {
            bignum_copy(prime, &candidate);
            return; // 소수이면 루프 종료
        }
    }
}
