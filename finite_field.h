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
#include <stdint.h>
#include <immintrin.h>

#define MILLER_RABIN_K 100 // Chances for false positive is 4^(-K) so if k=100 we get 6.22301528 * 10^(-61) :)
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

typedef struct{
    Fpe first;
    Fpe second;
}Pair;

u128 epow(u128 a, u128 d, u128 modulus);
bool miller_rabin_test(u128 n);
bool is_qr(Fpe a);
Pair tonelli_shanks(Fpe a);
u128 gcd(u128 a, u128 b);
Triplet ext_euclid(i128 a, i128 b);
Fpe fp_init(u128 value, u128 prime);
bool fp_equals(Fpe a, Fpe b);
Fpe fp_add(Fpe a, Fpe b);
Fpe fp_neg(Fpe a);
Fpe fp_sub(Fpe a, Fpe b);
Fpe fp_mul(Fpe a, Fpe b);
Fpe fp_smul(u128 a, Fpe b);
Fpe fp_pow(Fpe a, u128 b);
Fpe fp_inverse(Fpe a);
void print_int128(__int128_t value);


#endif
