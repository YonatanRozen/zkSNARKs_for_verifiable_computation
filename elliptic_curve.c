#include "elliptic_curve.h"

u128 calc_p(i128 z){
    return (u128)(36*pow(z,4) + 36*pow(z,3) + 24*pow(z,2) + 6*z + 1);
}

u128 calc_t(i128 z){
    return (u128)(6*pow(z,2) + 1);
}

EC ec_init(u128 mbits){

    // Start with z ≈ 2^(mbits/4) 
    i128 z = pow(2, mbits / 4);
    // Loop until ceil(log_2(p(-z)))=m
    u128 p, n;
    while (true){
        u128 t = calc_t(z);
        p = calc_p(-z);
        n = p + 1 - t;
        if (miller_rabin_test(p) && miller_rabin_test(n))   // Modulus p is prime and curver order n is prime => break. 
            break;
        p = calc_p(z);
        n = p + 1 - t;
        if (miller_rabin_test(p) && miller_rabin_test(n))   // Try again with +z. 
            break;
        // If didn't succeed, increase z; 
        z++;
    }

    // Initialize BN-curve params: we know a=1, need to find b.
    Fpe a = fp_init(1, p);
    Fpe b = fp_init(0, p);
    Fpe one = fp_init(1, p);
    EC ec = {.p = p, .n = n, .a = a, .b = b, .y = b};   // y = b doesn't matter (will be changed anyway).
    ECP g;
    do{
        do{
            b = fp_add(b, one);
        }while(!is_qr(fp_add(b, one)));
        // TODO: 
        // 1. Find y ∈ F[p] s.t y^2=b+1 (mod p)
        // 2. G <- (1,y) on the curve E: y^2=x^3+b
        // 3. If nG = Infinity point "O" then break.
        Pair roots = tonelli_shanks(fp_add(b, one));
        ec.b = b;
        ec.y = roots.first;
        g = ecp_init(one, ec.y, false, ec);
    }while (!ecp_smul(n, g).is_infinity);

    // Return an EC defined by (p,n,b,y).
    return ec;
}

ECP ecp_init(Fpe x, Fpe y, bool is_infinity, EC e){
    if (x.prime != e.p || y.prime != e.p || x.value >= e.p || y.value >= e.p){
        printf("Error: ecp_init: x or y aren't in F[p(e)]!\n");
        exit(EXIT_FAILURE);
    }else if (fp_pow(y, 2).value != fp_add(fp_pow(x, 3), fp_add(fp_mul(e.a, x), e.b)).value){
        printf("Error: ecp_init: (x,y) isn't on the eliptic curve e\n!");
        exit(EXIT_FAILURE);
    }

    ECP ecp = {.x = x, .y = y, .is_infinity = is_infinity, .e = e};
    return ecp;
}

bool ecp_equals(ECP p, ECP q){
    return (p.x.value == q.x.value && p.y.value == q.y.value) || 
           (p.is_infinity == q.is_infinity == true);
}

ECP ecp_neg(ECP p){
    ECP res = {.x = p.x, .y = fp_neg(p.y), .is_infinity = p.is_infinity};
    return res; 
}

ECP ecp_add(ECP p, ECP q){
    if (ecp_equals(p, ecp_neg(q))){ // Q = -P => Q + (-P) = P + (-P) = O
        ECP res = {.x = fp_init(1,2), .y = fp_init(3,4), .is_infinity = true}; // garbage in a and b, won't use them anyway.
        return res;
    }else if (p.is_infinity){
        return q;
    }else if (q.is_infinity){
        return p;
    }else if (ecp_equals(p, q)){
        // Here we assume that P = Q
        // lambda = (3x^2 + a) / (2y)
        Fpe lambda = fp_mul(fp_add(fp_smul(3, fp_pow(p.x, 2)), p.e.a), fp_inverse(fp_smul(2, p.y))); 
        
        Fpe x_r = fp_add(fp_pow(lambda, 2), fp_neg(fp_smul(2, p.x)));
        Fpe y_r = fp_add(fp_mul(lambda, fp_add(p.x, fp_neg(x_r))), fp_neg(p.y));
        ECP res = {.x = x_r, .y = y_r, .is_infinity = false, .e = p.e};
        return res; 
    }

    // Here we assume that P != +-Q and P != O and Q != O
    // lambda = (y_p - y_q) / (x_p - x_q)
    Fpe lambda = fp_mul(fp_add(p.y, fp_neg(q.y)), fp_inverse(fp_add(p.x, fp_neg(q.x))));
    Fpe x_r = fp_add(fp_pow(lambda, 2), fp_neg(fp_add(p.x, q.x))); // x = lambda^2 - x_p - x_q
    Fpe y_r = fp_add(fp_mul(lambda, fp_add(p.x, fp_neg(x_r))), fp_neg(p.y)); // y = lambda * (x_p - x_r) - y_p
    ECP res = {.x = x_r, .y = y_r, .is_infinity = false, .e = p.e};
    return res; 
}

ECP ecp_smul(u128 s, ECP p){
    if (p.is_infinity){
        return p;
    }

    ECP mul = p;
    ECP res = {.x = fp_init(1,2), .y = fp_init(3,4), .is_infinity = true, .e = p.e};

    while (s > 0){
        if (s & 1){
            // Sanity check: in the 1st time we encounter a '1' bit then 
            // res = O but mul != O (since mul=p and from 1st cond p != O). 
            // This means that ecp_add(res,mul) is like O + mul = mul and that's what we want!
            res = ecp_add(res, mul);
        }
        mul = ecp_add(mul, mul);
        s >>= 1;
    }

    return res;
}

ECP pairing(ECP p, ECP q){
    
}