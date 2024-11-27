#include <finite_field.h>

#define MILLER_RABIN_K 100 // Chances for false positive are 4^(-K) so if k=100 we get 6.22301528 * 10^(-61) :)

u64 efficient_pow(u64 a, u64 d, u64 modulus){
    u64 res = 1;
    u64 mul = a; 

    while (a > 0){
        if (a & 1){
            res = (u64)(res * mul) % modulus;
        }
        mul = (u64)(mul * mul) % modulus;
        a >>= 1;
    }
    return res;
}

bool miller_rabin_test(u64 n, unsigned int k){
    // definitely not a prime
    if (n <= 2 || n % 2 == 0)
        return false;

    // factor out powers of 2 from n to find s > 0 and d > 0 s.t n-1=2^(s)*d and d is odd. 
    u64 s = 0;
    u64 d = n-1;
    while (d % 2 == 0){
        s++;
        d /= 2;
    }

    for (size_t i = 0; i < k; i++)
    {
        u64 a = (u64)rand() % (n - 3) + 2;   // a should be in the range (2, n-2)
        u64 x = efficient_pow(a, d, n);
        if (x == 1){   // means a^d - 1 = 0 (mod n) so n divides one of the factors of a^(n-1) => \
                                            // n is a prime with high probability (satisfies fermats theorem for random base).
            return true;
        }
        u64 j = 1;
        while (j < s){
            x = (u64)(x * x) % n;
            if (x == n - 1){
                return true;
            }
        }
    }
    return false;
}

Fpe Fp_init(u64 value, u64 prime){
    if (!miller_rabin_test(prime, MILLER_RABIN_K)){
        printf("Error: Fp_init: 2nd arg isn't a prime with high probability!\n"); 
        exit(EXIT_FAILURE);
    }else if(value > prime){
        printf("Error: Fp_init: 1st arg isn't in range {0,...,prime}!\n");
        exit(EXIT_FAILURE);
    }
    Fpe elem;
    elem.value = value;
    elem.prime = prime;
    return elem;
}

Fpe Fp_add(Fpe a, Fpe b){
    if (a.prime != b.prime){
        printf("Error: Fp_add: trying to add elements from different fields!\n"); 
        exit(EXIT_FAILURE);
    }

    Fpe elem;
    elem.value = (a.value + b.value) % a.prime;
    elem.prime = a.prime;
    return elem;
}

Fpe Fp_neg(Fpe a){
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
    Fpe elem = Fp_neg(b);   // -b
    elem.value = (elem.value + a.value) % elem.prime;   // (-b + a) % p
    return elem;
}

Fpe Fp_mul(Fpe a, Fpe b){
    if (a.prime != b.prime){
        printf("Error: Fp_mul: trying to mul elements from different fields!\n"); 
        exit(EXIT_FAILURE);
    }
    Fpe elem;
    elem.value = (a.value * b.value) % a.prime;   
    elem.prime = a.prime;
    return elem;
}

Fpe Fp_inverse(Fpe a){
    if (a.value == 0){
        printf("Error: Fp_inverse: 0 has no inverse!\n"); 
        exit(EXIT_FAILURE);
    }   
}