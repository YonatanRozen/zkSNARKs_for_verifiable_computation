/* An eliptic curve type.
*/

#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <stdbool.h>
#include <immintrin.h>
#include <finite_field.h>

typedef struct{
    u128 a;
    u128 b;
    u128 p;
}EC;

typedef struct{
    u128 x;
    u128 y;
    bool is_infinity;
}ECP;

EC ec_init(u128 a, u128 b, u128 p);
ECP ecp_init(u128 x, u128 y, bool is_infinity, EC e);
bool ecp_equals(ECP p, ECP q);
ECP ec_add(ECP p, ECP q);
ECP ec_neg(ECP p);
ECP ec_double(ECP p);

#endif