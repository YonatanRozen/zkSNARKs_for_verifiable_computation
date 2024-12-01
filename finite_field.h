/* A finite field element type.
- It is a small type, so I don't use pointer management. 
*/

#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <stdbool.h>
#include <immintrin.h>

typedef __uint128_t u128;
typedef __int128_t i128;
// typedef __m256_u u256;
//typedef __m256i i256;

typedef struct{
    u128 first;
    u128 second;
    u128 third;
}Triplet;

typedef struct{
    u128 value;
    u128 prime;
}Fpe;

u128 efficient_pow(u128 a, u128 d, u128 modulus);
bool miller_rabin_test(u128 n, unsigned int k);
u128 gcd(u128 a, u128 b);
Fpe Fp_Init(u128 value, u128 prime);
Fpe Fp_add(Fpe a, Fpe b);
Fpe Fp_neg(Fpe a);
Fpe Fp_sub(Fpe a, Fpe b);
Fpe Fp_mul(Fpe a, Fpe b);
Fpe Fp_inverse(Fpe a);

#endif

