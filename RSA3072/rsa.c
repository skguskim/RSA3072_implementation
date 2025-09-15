#include "rsa.h"
#include <string.h>

#define RSA_SUCCESS 0
#define RSA_ERROR_NULL_ARGUMENT -1
#define RSA_ERROR_MESSAGE_TOO_LARGE -2

// RSA Ecryption
int rsa_encrypt(Bignum* ciphertext, const Bignum* message, const RSA_PublicKey* pub_key) {
	// input validation
    if (ciphertext == NULL || message == NULL || pub_key == NULL) {
        return RSA_ERROR_NULL_ARGUMENT;
    }
	
	// RSA validation (message length < n)
    if (bignum_compare(message, &pub_key->n) >= 0) {
        return RSA_ERROR_MESSAGE_TOO_LARGE;
    }

	Bignum m_reduced;
	bignum_mod(&m_reduced, message, &pub_key->n);
	bignum_mod_exp(ciphertext, &m_reduced, &pub_key->e, &pub_key->n);
	return RSA_SUCCESS;
}

// RSA Decryption
int rsa_decrypt(Bignum* message, const Bignum* ciphertext, const RSA_PrivateKey* priv_key) {
	// input validation
    if (message == NULL || ciphertext == NULL || priv_key == NULL) {
        return RSA_ERROR_NULL_ARGUMENT;
    }
	
	Bignum m_p, m_q, h;
	// m_p = c^dP mod p
	bignum_mod_exp(&m_p, ciphertext, &priv_key->dP, &priv_key->p);
	// m_q = c^dQ mod q
	bignum_mod_exp(&m_q, ciphertext, &priv_key->dQ, &priv_key->q);

	// h = q^(-1) * (m_p - m_q) mod p
	if (bignum_compare(&m_p, &m_q) < 0) {
		bignum_add(&h, &m_p, &priv_key->p);
		bignum_subtract(&h, &h, &m_q);
	}
	else {
		bignum_subtract(&h, &m_p, &m_q);
	}
	bn_mod_mul(&h, &h, &priv_key->qInv, &priv_key->p);
	bignum_divide(NULL, &h, &h, &priv_key->p);

    // m = m_q + q * h
    Bignum temp;
	bignum_multiply(&temp, &h, &priv_key->q);
	bignum_add(message, &m_q, &temp);
	return RSA_SUCCESS;
}


/* ========================================================================
 * 내부 유틸
 * ===================================================================== */
// 0 초기화 + 작은 상수 세팅
static void bignum_set_u32(Bignum* r, uint32_t v) {
	bignum_init(r); // r->size = 0으로 초기화

	if (v == 0) {
		// v가 0이면, size가 0인 상태 그대로 반환
		return;
	}

	r->limbs[0] = v;
	r->size = 1;
}

// a = a - v (v는 작음)
static void bignum_sub_u32(Bignum* a, uint32_t v) {
	Bignum t; bignum_set_u32(&t, v);
	Bignum r; bignum_subtract(&r, a, &t);
	bignum_copy(a, &r);
}

// g = gcd(a,b)
static void bignum_gcd(Bignum* g, const Bignum* a, const Bignum* b) {
	Bignum A, B, R, zero; bignum_copy(&A, a); bignum_copy(&B, b); bignum_set_u32(&zero, 0);
	while (bignum_compare(&B, &zero) != 0) {
		Bignum Q;
		bignum_divide(&Q, &R, &A, &B);	// R = A % B
		bignum_copy(&A, &B);
		bignum_copy(&B, &R);
	}
	bignum_copy(g, &A);
}

// L = Lcm(a,b) = (a/gcd(a,b)) * b
static void bignum_lcm(Bignum* l, const Bignum* a, const Bignum* b) {
	Bignum g, q, r, t;
	bignum_gcd(&g, a, b);
	bignum_divide(&q, &r, a, &g);	// q = a / g
	bignum_multiply(&t, &q, b);	// t = q * b
	bignum_copy(l, &t);
}

// inv = a^{-1} mod m  (gcd(a,m)=1 가정, 반환 1=성공, 0=실패)
static int bignum_modinv(Bignum* inv, const Bignum* a, const Bignum* m) {
	Bignum r0, r1, t0, t1, zero;
	bignum_copy(&r0, m);
	bignum_copy(&r1, a);
	bignum_set_u32(&t0, 0);
	bignum_set_u32(&t1, 1);
	bignum_set_u32(&zero, 0);

	while (bignum_compare(&r1, &zero) != 0) {
		Bignum q, r;
		bignum_divide(&q, &r, &r0, &r1); // q = r0 / r1, r = r0 % r1
		bignum_copy(&r0, &r1);
		bignum_copy(&r1, &r);

		// ---- fix: (q*t1) % m 로 축소한 뒤 t0에서 빼기 ----
		Bignum q_t1, q_t1_mod, tmp;
		bignum_multiply(&q_t1, &q, &t1);
		bignum_mod(&q_t1_mod, &q_t1, m);    // q_t1_mod = (q*t1) % m

		if (bignum_compare(&t0, &q_t1_mod) >= 0) {
			bignum_subtract(&tmp, &t0, &q_t1_mod);
		}
		else {
			Bignum diff;                    // diff = q_t1_mod - t0
			bignum_subtract(&diff, &q_t1_mod, &t0);
			bignum_subtract(&tmp, m, &diff); // tmp = m - diff
		}
		bignum_copy(&t0, &t1);
		bignum_copy(&t1, &tmp);
	}

	Bignum one; bignum_set_u32(&one, 1);
	if (bignum_compare(&r0, &one) != 0) { // gcd != 1
		bignum_set_u32(inv, 0);
		return 0;
	}
	bignum_copy(inv, &t0); // 0..m-1
	return 1;
}

// e 선택: 65537 기본, gcd(e, λ(n))=1 보장되도록 증가
static void pick_public_exponent(Bignum* e, const Bignum* lambda_n) {
	Bignum g, one, two; bignum_set_u32(e, 65537); bignum_set_u32(&one, 1); bignum_set_u32(&two, 2);
	while (1) {
		bignum_gcd(&g, e, lambda_n);
		if (bignum_compare(&g, &one) == 0) break;
		Bignum t; bignum_add(&t, e, &two); bignum_copy(e, &t);	// 홀수만
	}
}


/* ========================================================================
 * 키 생성
 * ===================================================================== */
void rsa_generate_keys(RSA_PublicKey* pub_key, RSA_PrivateKey* priv_key, const Bignum* p, const Bignum* q) {
	// n = p*q
	Bignum n; bignum_multiply(&n, p, q);

	// p-1, q-1
	Bignum p1, q1; bignum_copy(&p1, p); bignum_sub_u32(&p1, 1);
	bignum_copy(&q1, q); bignum_sub_u32(&q1, 1);

	// λ(n) = lcm(p-1, q-1)
	Bignum lambda_n; bignum_lcm(&lambda_n, &p1, &q1);

	// e 선택
	Bignum e; pick_public_exponent(&e, &lambda_n);
	
	// d = e^(-1) mod λ(n)
	Bignum d;
	if (!bignum_modinv(&d, &e, &lambda_n)) {
		// 이론상 pick_public_exponent가 보장하지만, 방어적으로 처리
		// e를 65537+2k로 올리며 재시도해도 되지만, 여기서는 0 세팅 후 반환
		bignum_set_u32(&d, 0);
	}

	// CRT 파라미터
	Bignum dP, dQ, qInv;
	bignum_mod(&dP, &d, &p1);	// d mod (p-1)
	bignum_mod(&dQ, &d, &q1);	// d mod (q-1)
	bignum_modinv(&qInv, q, p);	// q^(-1) mod p

	// 결과 저장
	bignum_copy(&pub_key->n, &n);
	bignum_copy(&pub_key->e, &e);

	bignum_copy(&priv_key->n, &n);
	bignum_copy(&priv_key->d, &d);
	bignum_copy(&priv_key->p, p);
	bignum_copy(&priv_key->q, q);
	bignum_copy(&priv_key->dP, &dP);
	bignum_copy(&priv_key->dQ, &dQ);
	bignum_copy(&priv_key->qInv, &qInv);
}
