#include <elliptic_curve.h>

EC ec_init(u128 a, u128 b, u128 p){
    if (!miller_rabin_test(p)){
        printf("Error: ec_init: 3rd arg isn't a prime with high probability!\n"); 
        exit(EXIT_FAILURE);
    }else if (a >= p || b >= p){
        printf("Error: ec_init: 1st or 2nd arg isn't in range {0...p-1}!\n"); 
        exit(EXIT_FAILURE);
    }else if ((4 * epow(a, 3, p) + 27 * epow(b, 2, p)) % p == 0){
        printf("Error: the curve is singular for a and b values!\n");
        exit(EXIT_FAILURE);
    }

    EC elliptic_curve = {.a = a, .b = b, .p = p};
    return elliptic_curve;
}

ECP ecp_init(u128 x, u128 y, bool is_infinity, EC e){
    if (x >= e.p || y >= e.p){
        printf("Error: ecp_init: x or y aren't in the range {0...p(e)}!\n");
        exit(EXIT_FAILURE);
    }else if (epow(y, 2, e.p) != (epow(x, 3, e.p) + (e.a * x) % e.p + e.b % e.p) % e.p){
        printf("Error: ecp_init: (x,y) isn't on the eliptic curve e\n!");
        exit(EXIT_FAILURE);
    }

    ECP ecp = {.x = x, .y = y, .is_infinity = is_infinity};
    return ecp;
}

bool ecp_equals(ECP p, ECP q){
    return (p.x == q.x && p.y == q.y) || (p.is_infinity == q.is_infinity == true);
}

ECP ec_add(ECP p, ECP q){
    if (ecp_equals(p, q)){
        ECP res = {.x = 1, .y = 2, .is_infinity = true}; // garbage in a and b, won't use them anyway.
        return res;
    }else if (p.is_infinity){
        return q;
    }else if (q.is_infinity){
        return p;
    }

}

ECP ec_neg(ECP p){

}

ECP ec_double(ECP p){

}