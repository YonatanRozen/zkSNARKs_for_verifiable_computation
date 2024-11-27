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

typedef u_int64_t u64;

typedef struct{
    u64 first;
    u64 second;
}Pair;

typedef struct{
    u64 value;
    u64 prime;
}Fpe;

u64 efficient_pow(u64 a, u64 d, u64 modulus);
bool miller_rabin_test(u64 n, unsigned int k);
u64 gcd(u64 a, u64 b);
Fpe Fp_Init(u64 value, u64 prime);
Fpe Fp_add(Fpe a, Fpe b);
Fpe Fp_neg(Fpe a);
Fpe Fp_sub(Fpe a, Fpe b);
Fpe Fp_mul(Fpe a, Fpe b);
Fpe Fp_inverse(Fpe a);

#endif

