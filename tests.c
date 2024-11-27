#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <assert.h>
#include <stdbool.h>

typedef u_int64_t u64;

typedef struct{
    u64 first;
    u64 second;
}Pair;

u64 efficient_pow(u64 a, u64 d, u64 modulus){
    u64 res = 1;
    u64 mul = a; 

    while (d > 0){
        if (d & 1){
            res = (res * mul) % modulus;
        }
        mul = (mul * mul) % modulus;
        d >>= 1;
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

u64 gcd(u64 a, u64 b){ 
    u64 r;
    while (a % b){
        r = a % b;
        a = b;
        b = r;
    }
    return r;
}

Pair ext_euclid(u64 a, u64 b){ 
    u64 unPrev = 1;
    u64 vnPrev = 0;
    u64 unCur = 0;
    u64 vnCur = 1;

    while (b != 0){
        u64 bn = a; // b
        u64 newB = a % b;
        a = b;
        b = newB;

        // Update coefficients
        u64 unNew = unPrev - bn * unCur;
        u64 vnNew = vnPrev - bn * vnCur;

        // Shift coefficients
        unPrev = unCur;
        vnPrev = vnCur;
        unCur = unNew;
        vnCur = vnNew;
    }
    
    Pair p;
    p.first = unPrev;
    p.second = vnPrev;
    return p;
}

void test_efficient_pow(){
    u64 a = 2;
    u64 d = 3;
    u64 n = 5;
    assert(efficient_pow(a, d, n) == 3);
}

void test_miller_rabin(){
    assert(miller_rabin_test(104729, 100));
    assert(!miller_rabin_test(2, 100));
    assert(!miller_rabin_test(20, 100));
    assert(!miller_rabin_test(40124014202, 100));
    assert(!miller_rabin_test(12419082402, 100));
    assert(miller_rabin_test(11, 100));
    assert(miller_rabin_test(79, 100));
    assert(!miller_rabin_test(10, 100));
    assert(miller_rabin_test(1319, 100));
}

void test_gcd(){
    assert(gcd(3,5)==1);
    assert(gcd(3,9)==3);
    assert(gcd(14,70)==14);
    assert(gcd(128,48)==16);
    assert(gcd(12,55)==1);
    assert(gcd(2,5)==1);
    assert(gcd(7,5)==1);
    assert(gcd(7,90)==1);
    assert(gcd(2,8)==2);
}

void test_ext_euclid(){
    Pair res1 = ext_euclid(23, 70);
    printf("Ext euclid 1: (%lu,%lu), %lu\n", res1.first, res1.second, (23*res1.first + 70*res1.second));
}

int main(){
    test_efficient_pow();
    test_miller_rabin();
    test_gcd();
    test_ext_euclid();
    return 0;
}