#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>

#define MILLER_RABIN_K 100 // Chances for false positive is 4^(-K) so if k=100 we get 6.22301528 * 10^(-61) :)
typedef __uint128_t u128;
typedef __int128_t i128;

typedef struct{
    i128 first;
    i128 second;
    i128 third;
}Triplet;

u128 epow(u128 a, u128 d, u128 modulus){
    u128 res = 1;
    u128 mul = a; 

    while (d > 0){
        if (d & 1){
            res = (res * mul) % modulus;
        }
        mul = (mul * mul) % modulus;
        d >>= 1;
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

i128 gcd(i128 a, i128 b){ 
    i128 r;
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

void test_efficient_pow(){
    u128 a = 2;
    u128 d = 3;
    u128 n = 5;
    assert(epow(a, d, n) == 3);
}

void test_miller_rabin(){
    assert(miller_rabin_test(104729));
    assert(!miller_rabin_test(2));
    assert(!miller_rabin_test(20));
    assert(!miller_rabin_test(40124014202));
    assert(!miller_rabin_test(12419082402));
    assert(miller_rabin_test(11));
    assert(miller_rabin_test(79));
    assert(!miller_rabin_test(10));
    assert(miller_rabin_test(1319));
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

// Prints a 128 bit integer.
// I didn't write it myself but it works :)
void print_int128(__int128_t value) {
    if (value < 0) {
        putchar('-');
        value = -value;
    }
    if (value > UINT64_MAX) {
        print_int128(value / 1000000000000000000);
        printf("%018lu", (uint64_t)(value % 1000000000000000000));
    } else {
        printf("%lu", (uint64_t)value);
    }
}

void test_ext_euclid() {
    __int128_t a1 = 104; // Example values
    __int128_t b1 = 47;

    Triplet result1 = ext_euclid(a1, b1);
    assert(a1 * result1.second + b1 * result1.third == result1.first);

    __int128_t a2 = 23423526534523; // Example values
    __int128_t b2 = 477961987446346346;

    Triplet result2 = ext_euclid(a2, b2);
    assert(a2 * result2.second + b2 * result2.third == result2.first);

    __int128_t a3 = 104123598246245; // Example values
    __int128_t b3 = 7879532459834234385;

    Triplet result3 = ext_euclid(a3, b3);
    assert(a3 * result3.second + b3 * result3.third == result3.first);

    __int128_t a4 = 1; // Example values
    __int128_t b4 = 2;

    Triplet result4 = ext_euclid(a4, b4);
    assert(a4 * result4.second + b4 * result4.third == result4.first);

    __int128_t a5 = 5; // Example values
    __int128_t b5 = 3;

    Triplet result5 = ext_euclid(a5, b5);
    assert(a5 * result5.second + b5 * result5.third == result5.first);
    // Loooooooooooooong print for sanity check
    // printf("GCD(");
    // print_int128(a);
    // printf(",");
    // print_int128(b);
    // printf(") = ");
    // print_int128(result.first);
    // printf(" = ");
    // print_int128(result.second);
    // printf(" * ");
    // print_int128(a);
    // printf(" + ");
    // print_int128(result.third);
    // printf(" * ");
    // print_int128(b);
    // printf("\n");
}


int main(){
    test_efficient_pow();
    test_miller_rabin();
    test_gcd();
    test_ext_euclid();
    return 0;
}