#include "finite_field.h"

// Prints a 128 bit integer.
// I didn't write it myself but it works :)
void print_int128(__int128_t value) {
    if (value < 0) {
        putchar('-');
        value = -value;
    }
    if (value > UINT64_MAX) {
        print_int128(value / 1000000000000000000);
        printf("%018lu\n", (uint64_t)(value % 1000000000000000000));
    } else {
        printf("%lu\n", (uint64_t)value);
    }
}

u128 epow(u128 a, u128 d, u128 modulus){
    u128 res = 1;
    u128 mul = a; 

    if (modulus != 0){
        while (d > 0){
            if (d & 1){
                res = (res * mul) % modulus;
            }
            mul = (mul * mul) % modulus;
            d >>= 1;
        }
    }else{
        while (d > 0){
            if (d & 1){
                res = (res * mul);
            }
            mul = (mul * mul);
            d >>= 1;
        }
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

Fpe fp_init(u128 value, u128 prime){
    if (!miller_rabin_test(prime)){
        printf("Error: fp_init: 2nd arg isn't a prime with high probability!\n"); 
        exit(EXIT_FAILURE);
    }else if(value >= prime){
        printf("Error: fp_init: 1st arg isn't in range {0,...,p-1}!\n");
        exit(EXIT_FAILURE);
    }
    
    Fpe elem;
    elem.value = value;
    elem.prime = prime;
    return elem;
}

bool fp_equals(Fpe a, Fpe b){
    if (a.prime != b.prime){
        printf("Error: fp_equals: a and b aren't in the same field!\n"); 
        exit(EXIT_FAILURE);
    }

    return a.value == b.value;
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
        printf("Error: fp_mul: a or b aren't in F[p]!\n");
        exit(EXIT_FAILURE);
    }

    Fpe elem;
    elem.value = (a.value * b.value) % a.prime;
    elem.prime = a.prime;
    return elem;
}

Fpe fp_smul(u128 a, Fpe b){
    Fpe elem;
    elem.value = (a * b.value) % b.prime;   
    elem.prime = b.prime;
    return elem;
}

Fpe fp_pow(Fpe a, u128 b){
    Fpe elem;
    elem.value = epow(a.value, b, a.prime);   
    elem.prime = a.prime;
    return elem;
}

// The algorithm uses euler's criterion internally to check if a is a QR mod p.
bool is_qr(Fpe a){
    if (a.value == 0){
        return false;
    }
    Fpe one = fp_init(1, a.prime);
    return fp_equals(fp_pow(a, (a.prime - 1) / 2), one);
}

Fpe fp_inverse(Fpe a){
    if (a.value == 0){
        printf("Error: Fp_inverse: 0 has no inverse!\n"); 
        exit(EXIT_FAILURE);
    }

    Triplet extEuclidRes = ext_euclid(a.value, a.prime); // find k, l s.t a * k + p * l = 1
    Fpe res = {.value = extEuclidRes.second, .prime = a.prime};
    return res; 
}

Pair tonelli_shanks(Fpe a){
    if (!is_qr(a)){
        printf("Error: fp_qr: a isn't a qr!\n");
        exit(EXIT_FAILURE);
    }  

    Pair res;
    u128 p = a.prime;
    Fpe one = fp_init(1, p);
    // Case 1: p=2 => root is a
    if (p == 2){
        res.first = a;
        res.second = a;
        return res;
    }else if (p % 4 == 3){    // Case 2: p=4n+3 for some n => root is +-a^((p+1)/4)
        Fpe r = fp_pow(a, (p + 1) / 4);
        res.first = r;
        res.second = fp_neg(r);
        return res;
    }else if (p % 8 == 5){    // Case 3: p=8n+5 for some n => root is +-av(i-1) where v=(2a)^((p-5)/8), i=2av^2
        Fpe v = fp_pow(fp_smul(2, a), (p - 5) / 8);
        Fpe i = fp_smul(2, fp_mul(a, fp_pow(v, 2)));
        res.first = fp_mul(a, fp_mul(v, fp_add(i, fp_neg(one))));
        res.second = fp_neg(res.first);
        return res;
    }

    // Case 4: p=8n+1 for some n => very long computation. 
    // 1st step: factor p to 2^e*q + 1 where q is odd
    u128 e = 0;
    u128 q = p - 1;
    while (q % 2 == 0){
        e++;
        q /= 2;
    }
    // 2nd step: generate random x s.t 1 < x < p until x^(q*(2^(e-1)))=1 (mod p) satisfies 
    Fpe x = fp_init((u128)rand() % (p - 2) + 2, p);   // x should be in the range [2, p-1]
    u128 expo = (u128)epow(2, e-1, 0);
    Fpe z = fp_pow(x, q);
    Fpe z2e = fp_pow(z, expo);
    while (fp_equals(z, one)){
        x = fp_init((u128)rand() % (p - 2) + 2, p);
        z2e = fp_pow(z, expo);
        
    }

    // 3rd step: set y <- z, r <- e, x <- a^((q-1)/2) (mod p), v <- ax (mod p), w <- vx (mod p)
    Fpe y = z;
    u128 r = e;
    x = fp_pow(a, (q - 1) / 2);
    Fpe v = fp_mul(a, x);
    Fpe w = fp_mul(v, x);
    // 4th step: if w = 1 (mod p) return +-v as the root. else enter \ continue loop
    while (!fp_equals(w, one)){
        // 5th step: find smallest k s.t w^(2^k)=1 (mod p)
        u128 k = 1; // We already know that w = w^(2^0) != 1 (mod p)
        u128 temp = epow(2, k, 0);
        Fpe w_2k = fp_pow(w, temp);
        while (!fp_equals(w_2k, one)){
            k++;
            temp = epow(2, k, 0);
            w_2k = fp_pow(w, temp);
        }
        // 6th step: set d <- y^(2^(r-k-1)) (mod p), y <- d^2 (mod p), r <- k, v <- dv (mod p), w <- wy (mod p)
        temp = epow(2, r - k - 1, 0);
        Fpe d = fp_pow(y, temp);
        y = fp_pow(d, 2);
        r = k;
        v = fp_mul(d, v);
        w = fp_mul(w, y);

        //7th step: Back to step 4...
    } 
    res.first = v;
    res.second = fp_neg(v);
    return res; 
}