#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include "finite_field.h"
#include "elliptic_curve.h"

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

void test_tonelli_shanks(){
    Fpe a1 = fp_init(58, 101);
    Fpe r1_1 =  fp_init(82, 101);
    Fpe r1_2 = fp_init(19, 101);
    Pair res1 = tonelli_shanks(a1);
    // printf("1st root: ");
    // print_int128(res1.first.value);
    // printf("2nd root: ");
    // print_int128(res1.second.value);
    assert((fp_equals(res1.first, r1_1) && fp_equals(res1.second, r1_2)) ||
           ((fp_equals(res1.first, r1_2) && fp_equals(res1.second, r1_1))));
    
    Fpe a2 = fp_init(111, 113);
    Fpe r2_1 =  fp_init(87, 113);
    Fpe r2_2 = fp_init(26, 113);
    Pair res2 = tonelli_shanks(a2);
    // printf("1st root: ");
    // print_int128(res2.first.value);
    // printf("2nd root: ");
    // print_int128(res2.second.value);
    assert((fp_equals(res2.first, r2_1) && fp_equals(res2.second, r2_2)) ||
           ((fp_equals(res2.first, r2_2) && fp_equals(res2.second, r2_1))));
}

void test_elliptic_curves_initialization(){
    
}

void test_other_things(){
    // Test if conversion to double happens in pow
    u128 m = 41;
    u128 y = epow(2, m / 4, 0);
    //print_int128(y);
    assert(y == 1024);
}

int main(){
    test_efficient_pow();
    test_miller_rabin();
    test_gcd();
    test_ext_euclid();
    test_tonelli_shanks();
    test_elliptic_curves_initialization();
    test_other_things();
    return 0;
}