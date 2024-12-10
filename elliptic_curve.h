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
#include <pbc/pbc.h>

typedef struct{
    u128 p;
    u128 n; 
    Fpe a;
    Fpe b; 
    Fpe y; 
}EC;

typedef struct{
    Fpe x;
    Fpe y;
    bool is_infinity;
    EC e;
}ECP;

EC ec_init(u128 mbits);
ECP ecp_init(Fpe x, Fpe y, bool is_infinity, EC e);
bool ecp_equals(ECP p, ECP q);
ECP ecp_neg(ECP p);
ECP ecp_add(ECP p, ECP q);
ECP ecp_smul(u128 s, ECP p);
ECP pairing(ECP p, ECP q);

#endif