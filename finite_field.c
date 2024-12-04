#include <finite_field.h>

u128 epow(u128 a, u128 d, u128 modulus){
    u128 res = 1;
    u128 mul = a; 

    while (a > 0){
        if (a & 1){
            res = (u128)(res * mul) % modulus;
        }
        mul = (u128)(mul * mul) % modulus;
        a >>= 1;
    }
    return res;
}

bool miller_rabin_test(u128 n){
    // definitely not a prime
    if (n <= 2 || n % 2 == 0)
        return false;

    // factor out powers of 2 from n to find s > 0 and d > 0 s.t n-1=2^(s)*d and d is odd. 
    u128 s = 0;
    u128 d = n-1;
    while (d % 2 == 0){
        s++;
        d /= 2;
    }

    for (size_t i = 0; i < MILLER_RABIN_K; i++)
    {
        u128 a = (u128)rand() % (n - 3) + 2;   // a should be in the range (2, n-2)
        u128 x = epow(a, d, n);
        if (x == 1){   // means a^d - 1 = 0 (mod n) so n divides one of the factors of a^(n-1) => \
                                            // n is a prime with high probability (satisfies fermats theorem for random base).
            return true;
        }
        u128 j = 1;
        while (j < s){
            x = (u128)(x * x) % n;
            if (x == n - 1){
                return true;
            }
        }
    }
    return false;
}

u128 gcd(u128 a, u128 b){
    u128 r;
    while (a % b){
        r = a % b;
        a = b;
        b = r;
    }

    return r;
}

Triplet ext_euclid(i128 a, i128 b){
    if (a == 0){
        Triplet res = {.first = b, .second = 0, .third = 1};
        return res;
    }

    i128 unPrev = 1;
    i128 vnPrev = 0;
    i128 unCurr = 0;
    i128 vnCurr = 1;

    while (b != 0){
        i128 qn = a / b;
        i128 rn = a % b;

        i128 unNew = unPrev - qn * unCurr;
        i128 vnNew = vnPrev - qn * vnCurr;

        a = b;
        b = rn;

        unPrev = unCurr;
        vnPrev = vnCurr;
        unCurr = unNew;
        vnCurr = vnNew;
    }

    Triplet res = {.first = a, .second = unPrev, .third = vnPrev};
    return res;
}

Fpe Fp_init(u128 value, u128 prime){
    if (!miller_rabin_test(prime)){
        printf("Error: Fp_init: 2nd arg isn't a prime with high probability!\n"); 
        exit(EXIT_FAILURE);
    }else if(value >= prime){
        printf("Error: Fp_init: 1st arg isn't in range {0,...,p-1}!\n");
        exit(EXIT_FAILURE);
    }
    Fpe elem;
    elem.value = value;
    elem.prime = prime;
    return elem;
}

Fpe fp_add(Fpe a, Fpe b){
    if (a.prime != b.prime){
        printf("Error: Fp_add: trying to add elements from different fields!\n"); 
        exit(EXIT_FAILURE);
    }

    Fpe elem;
    elem.value = (a.value + b.value) % a.prime;
    elem.prime = a.prime;
    return elem;
}

Fpe fp_neg(Fpe a){
    Fpe elem;
    elem.value = a.prime - a.value;     // For example in p=5 -2 = 3 = 5-2 
    elem.prime = a.prime;
    return elem;
}

Fpe Fp_sub(Fpe a, Fpe b){
    if (a.prime != b.prime){
        printf("Error: Fp_sub: trying to sub elements from different fields!\n"); 
        exit(EXIT_FAILURE);
    }
    Fpe elem = fp_neg(b);   // -b
    elem.value = (elem.value + a.value) % elem.prime;   // (-b + a) % p
    return elem;
}

Fpe fp_mul(Fpe a, Fpe b){
    if (a.prime != b.prime){
        printf("Error: Fp_mul: trying to mul elements from different fields!\n"); 
        exit(EXIT_FAILURE);
    }
    Fpe elem;
    elem.value = (a.value * b.value) % a.prime;   
    elem.prime = a.prime;
    return elem;
}

Fpe fp_inverse(Fpe a){
    if (a.value == 0){
        printf("Error: Fp_inverse: 0 has no inverse!\n"); 
        exit(EXIT_FAILURE);
    }
}