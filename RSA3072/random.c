#include "rsa.h"
#include "bignum.h"
#include <stddef.h>
#include <windows.h>
#include <bcrypt.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#pragma comment(lib, "bcrypt.lib")

// 버퍼를 무작위 바이트로 채움
// 성공 시 0, 인자 유효하지 않을 경우 -1, 실패 시 -2
int fill_with_random_bytes(uint32_t* buffer, const size_t size) {
    if (buffer == NULL || size == 0) {
        return -1;
    }
    NTSTATUS st = BCryptGenRandom(
        NULL,
        buffer,
        (ULONG)size,
        BCRYPT_USE_SYSTEM_PREFERRED_RNG
    );

    if (!BCRYPT_SUCCESS(st)) {
        return -2;
    }

    return 0;
}

// 지정된 바이트 크기의 Bignum 생성
int generate_secure_random(Bignum* result, const size_t size) {
    if (result == NULL || size == 0) {
        return -1;
    }

    // 메모리 할당 실패
    uint32_t* temp_buffer = (uint32_t*)malloc(size);
    if (temp_buffer == NULL) {
        return -3; 
    }

    // 버퍼를 무작위 바이트로 채우기
    int rc = fill_with_random_bytes(temp_buffer, size);
    if (rc != 0) {
        free(temp_buffer);
        return rc;
    }

    //바이트 배열에서 Bignum을 생성하는함수 호출
    rc = bignum_from_binary(result, temp_buffer, size);

    // 임시 버퍼를 안전하게 삭제하고 해제
    SecureZeroMemory(temp_buffer, size);
    free(temp_buffer);

    return rc; //성공 시 0
}

// 3072비트 소수 후보 난수를 생성합니다.
// 생성된 후보는 반드시 3072비트 길이(MSB=1)이고 홀수(LSB=1)임이 보장됩니다.
// 성공 시 0, 실패 시 0이 아닌 값.
int random_3072_candidate(Bignum* out) {
    const size_t BYTES = 384;

    int rc = generate_secure_random(out, BYTES);
    if (rc != 0) {
        return rc;
    }

    /// 숫자가 정확히 3072비트이고 홀수가 되도록 비트 조정
    bn_set_bit_local(out, 3071); // 최상위 비트(MSB)를 1로 설정
    bn_set_bit_local(out, 0);    // 최하위 비트(LSB)를 1로 설정 (홀수로 만듦)

    return 0;
}

// OAEP/PSS 패딩을 위한 32바이트(256비트) 무작위 시드를 생성합니다.
// 성공 시 0, 실패 시 0이 아닌 값.
int random_oaep_seed32(Bignum* out) {
    return generate_secure_random(out, 32);
}
