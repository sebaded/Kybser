#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include "inverse.h"

#define PKBYTES 800
#define SKBYTES 1632

#define CTBYTES 1024
#define SSBYTES 32

#define KYBER_N 256
#define KYBER_Q 3329
#define KYBER_K 2

#define KYBER_SYMBYTES 32


typedef struct{
  poly vec[KYBER_K];
} polyvec;


// kyber functions
extern void pqcrystals_kyber512_ref_keypair(uint8_t* pk, uint8_t* sk);
extern void pqcrystals_kyber512_ref_enc(uint8_t* ct, uint8_t* ss, uint8_t* pk);
extern void pqcrystals_kyber512_ref_dec(uint8_t* ss, uint8_t* ct, uint8_t* sk);
extern void pqcrystals_kyber512_ref_poly_frombytes(poly* r, const uint8_t a[CTBYTES-640]);
extern void pqcrystals_kyber512_ref_poly_tobytes(const uint8_t a[CTBYTES-640], poly* r);
extern void pqcrystals_kyber512_ref_enc_derand(uint8_t* ct, uint8_t* ss, uint8_t* pk, uint8_t* coins);
extern void randombytes(uint8_t* out, size_t outlen);
extern int16_t pqcrystals_kyber512_ref_barrett_reduce(int16_t a);
extern int16_t pqcrystals_kyber512_ref_montgomery_reduce(int32_t a);
extern void pqcrystals_kyber512_ref_poly_sub(poly *r, const poly* a, const poly* b);
extern void pqcrystals_kyber512_ref_poly_reduce(poly* r);
extern void pqcrystals_kyber512_ref_polyvec_ntt(polyvec* r);
extern void pqcrystals_kyber512_ref_poly_ntt(poly* r);
extern void pqcrystals_kyber512_ref_poly_basemul_montgomery(poly* r, const poly* a, const poly* b);
extern void pqcrystals_kyber512_ref_poly_tomont(poly* r);
extern void pqcrystals_kyber512_ref_poly_invntt_tomont(poly* r);
extern void pqcrystals_kyber512_ref_polyvec_invntt_tomont(polyvec* r);
extern void pqcrystals_kyber512_ref_poly_frommsg(poly *r, const uint8_t msg[KYBER_SYMBYTES]);
extern void pqcrystals_kyber512_ref_polyvec_decompress(polyvec *r, const uint8_t a[640]);
extern void pqcrystals_kyber512_ref_polyvec_frombytes(polyvec *r, const uint8_t a[768]);

// modified encryption to get msg
void enc_mod(uint8_t* ct, uint8_t* ss, uint8_t* pk, uint8_t* msg){
    uint8_t coins[KYBER_SYMBYTES];
    randombytes(coins, KYBER_SYMBYTES);
    memcpy(msg, &coins, KYBER_SYMBYTES);
    pqcrystals_kyber512_ref_enc_derand(ct, ss, pk, coins);
}


long int time_enc(uint8_t* ss, uint8_t* ct, uint8_t* sk, int loops) {
    long int sum = 0;
    struct timespec start, end;

    for(int i = 0; i < loops; i++) {
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);   // CLOCK_PROCESS_CPUTIME_ID , CLOCK_MONOTONIC, CLOCK_THREAD_CPUTIME_ID  /** MONOTONIC does not work on Lin's computer but THREAD_CPUTIME could work perfectly **/
        pqcrystals_kyber512_ref_dec(ss, ct, sk);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);

        sum += (end.tv_nsec - start.tv_nsec + (end.tv_sec - start.tv_sec)*1000000000);
    }
    return sum/loops;
}



int time_difference(uint8_t* ct, uint8_t* compare_ct, uint8_t* sk) {
    uint8_t ss[SSBYTES];
    long int t1, t2;
    int correct = 0;
    int loop = 100;

    for(int i = 0; i < loop; i++){
        time_enc(ss, ct, sk, 1); // warmup  helps for some reason

        t2 = time_enc(ss, ct, sk, 1);
        t1 = time_enc(ss, compare_ct, sk, 1);

        if(t2 < t1)
            correct++;
    }
 //     printf("%lu %lu %u \n", t1, t2, (((100*correct)/loop)));

    if(((100*correct)/loop) >= 80){
        return 1;
    }
    else if (((100*correct)/loop) <= 55){
        return 0;
    }
    else{
        return time_difference(ct, compare_ct, sk);
    }
}




int16_t recover_error_coeff(uint8_t* ct, uint8_t m[KYBER_SYMBYTES], uint8_t* sk, size_t i){
    poly v, v_ch;
    uint8_t cct[CTBYTES];
    uint8_t compare_ct[CTBYTES];

    memcpy(&cct, ct, CTBYTES);
    memcpy(&compare_ct, ct, CTBYTES);
    compare_ct[640] += 1;

    pqcrystals_kyber512_ref_poly_frombytes(&v, &ct[640]);
    memcpy(&v_ch, &v, sizeof(poly));

    int lower = 0;
    int upper = 1200;

    while(lower <= upper) {

        int middle = (int)((lower + upper)/2);
        v_ch.coeffs[i] = (v.coeffs[i] + middle) % KYBER_Q;
   //     printf("v: %u v_: %u ", v.coeffs[i], v_ch.coeffs[i]);
        pqcrystals_kyber512_ref_poly_tobytes(&cct[640], &v_ch);

        if(time_difference(cct, compare_ct, sk)) {

            v_ch.coeffs[i] = (v_ch.coeffs[i] - 1) % KYBER_Q;
     //       printf("v: %u v_: %u ", v.coeffs[i], v_ch.coeffs[i]);
            pqcrystals_kyber512_ref_poly_tobytes(&cct[640], &v_ch);

            if(!time_difference(cct, compare_ct, sk)){

                int b = m[i/8] >> (i%8) & 1;
                return 833 - b - middle;
            }
            else
                upper = middle - 1;
        }
        else
            lower = middle + 1;
    }
}




void revmont(poly* r){
    for(int i = 0; i < KYBER_N; i++)
        r->coeffs[i] = pqcrystals_kyber512_ref_barrett_reduce(pqcrystals_kyber512_ref_montgomery_reduce(r->coeffs[i]));
}




void correct_repr(poly* r){
    for(int i = 0; i < KYBER_N; i++)
        r->coeffs[i] = (r->coeffs[i] + KYBER_Q) % KYBER_Q;
}





int main(){
    char command[3];

    uint8_t pk[PKBYTES];
    uint8_t sk[SKBYTES];
    uint8_t ss[SSBYTES];

    uint8_t ct1[CTBYTES];
    uint8_t ct2[CTBYTES];

    uint8_t hm1[KYBER_SYMBYTES];
    uint8_t hm2[KYBER_SYMBYTES];

    pqcrystals_kyber512_ref_keypair(pk, sk);

    enc_mod(ct1, ss, pk, hm1);
    enc_mod(ct2, ss, pk, hm2);


    poly e1;
    poly e2;

    clock_t start = clock();

    for(int i = 0; i < KYBER_N; i++){
        e1.coeffs[i] = recover_error_coeff(ct1, hm1, sk, i);
        printf("recovered e1[%i]: %i\n", i, e1.coeffs[i]);
    }

    for(int i = 0; i < KYBER_N; i++){
        e2.coeffs[i] = recover_error_coeff(ct2, hm2, sk, i);
        printf("recovered e2[%i]: %i\n", i, e2.coeffs[i]);
    }

    printf("Time to recover Error_polynomials: %fs", (double)(clock()-start)/CLOCKS_PER_SEC);
    printf("\n");


    polyvec s;
    polyvec u1, u2;

    poly s0_calc, s1_calc, v1, v2, m1, m2;
    poly y1, y2;

    poly u1_inv;

    // get v1 and v2 from respective ct
    pqcrystals_kyber512_ref_poly_frombytes(&v1, ct1+640);
    pqcrystals_kyber512_ref_poly_frombytes(&v2, ct2+640);

    // get m1 and m2 from respective hm
    pqcrystals_kyber512_ref_poly_frommsg(&m1, hm1);
    pqcrystals_kyber512_ref_poly_frommsg(&m2, hm2);

    // get u1 and u2 from respective ct
    pqcrystals_kyber512_ref_polyvec_decompress(&u1, ct1);
    pqcrystals_kyber512_ref_polyvec_decompress(&u2, ct2);

    // get s from sk
    pqcrystals_kyber512_ref_polyvec_frombytes(&s, sk);


    // y1 = v1 - E1 - m1
    pqcrystals_kyber512_ref_poly_sub(&y1, &v1, &e1);
    pqcrystals_kyber512_ref_poly_sub(&y1, &y1, &m1);
    pqcrystals_kyber512_ref_poly_reduce(&y1);

    // y2 = v2 - E2 - m2
    pqcrystals_kyber512_ref_poly_sub(&y2, &v2, &e2);
    pqcrystals_kyber512_ref_poly_sub(&y2, &y2, &m2);
    pqcrystals_kyber512_ref_poly_reduce(&y2);


    // u1_inv = u_10^-1
    inverse(&u1_inv, &(u1.vec[0]));


    // s1 = (y2 - y1 * u_10^-1 * u_20) * (u_21 - u_10^-1 * u_11 * u_20)^-1
    // s1 = lefty * rightu^-1
    poly left, lefty;

    pqcrystals_kyber512_ref_polyvec_ntt(&u1);
    pqcrystals_kyber512_ref_polyvec_ntt(&u2);
    pqcrystals_kyber512_ref_poly_ntt(&u1_inv);
    pqcrystals_kyber512_ref_poly_ntt(&y1);


    // left = (y2 - y1 * u_10^-1 * u_20)
    // lefty = (y1 * left)  left = (u_10^-1 * u_20)
    // lefty = (y2 - lefty)
    pqcrystals_kyber512_ref_poly_basemul_montgomery(&left, &u1_inv, &u2.vec[0]); // left = u_10^-1 * u_20
    pqcrystals_kyber512_ref_poly_tomont(&left);
    pqcrystals_kyber512_ref_poly_basemul_montgomery(&lefty, &left, &y1);  // lefty = y1 * left
    pqcrystals_kyber512_ref_poly_invntt_tomont(&lefty);
    pqcrystals_kyber512_ref_poly_sub(&lefty, &y2, &lefty);  // left = y2 - left
    pqcrystals_kyber512_ref_poly_reduce(&lefty);


    // right = (u_21 - u_10^-1 * u_11 * u_20)
    // rightu = (u_10^-1 * right)  right = (u_11 * u_20)
    // rightu = (u_21 - rightu)
    poly right, rightu;

    pqcrystals_kyber512_ref_poly_basemul_montgomery(&right, &u1.vec[1], &u2.vec[0]);  // right = u_11 * u_20
    pqcrystals_kyber512_ref_poly_tomont(&right);
    pqcrystals_kyber512_ref_poly_basemul_montgomery(&rightu, &right, &u1_inv);  // rightu = right * u_10^-1
    pqcrystals_kyber512_ref_poly_invntt_tomont(&rightu);
    pqcrystals_kyber512_ref_polyvec_invntt_tomont(&u2);
    revmont(&u2.vec[0]);
    revmont(&u2.vec[1]);
    pqcrystals_kyber512_ref_poly_sub(&rightu, &u2.vec[1], &rightu);  // right = u_21 - right
    pqcrystals_kyber512_ref_poly_reduce(&rightu);


    // right_inv = rightu^-1
    poly right_inv;
    inverse(&right_inv, &rightu);

    // s1 = lefty * right_inv
    pqcrystals_kyber512_ref_poly_ntt(&lefty);
    pqcrystals_kyber512_ref_poly_ntt(&right_inv);
    pqcrystals_kyber512_ref_poly_basemul_montgomery(&s1_calc, &lefty, &right_inv);
    pqcrystals_kyber512_ref_poly_tomont(&s1_calc);


    // s0 = (y1 - s1 * u_11 ) * u_10^-1
    // s0 = (y1 - su1) * u_10^-1  su1 = (s1 * u_11)
    poly su1;

    pqcrystals_kyber512_ref_poly_basemul_montgomery(&su1, &s1_calc, &u1.vec[1]);  // su1 = s1 * u_11
    pqcrystals_kyber512_ref_poly_invntt_tomont(&su1);
    pqcrystals_kyber512_ref_poly_invntt_tomont(&y1);
    revmont(&y1);

    pqcrystals_kyber512_ref_poly_sub(&su1, &y1, &su1);   // su1 = (y1 - s1 * u_11)
    pqcrystals_kyber512_ref_poly_reduce(&su1);

    pqcrystals_kyber512_ref_poly_ntt(&su1);
    pqcrystals_kyber512_ref_poly_basemul_montgomery(&s0_calc, &su1, &u1_inv); // s0 = su * u1_inv

    pqcrystals_kyber512_ref_poly_tomont(&s0_calc);

    // convert negative coefficients to positive ones
    correct_repr(&s0_calc);
    correct_repr(&s1_calc);


    printf("s0_rec:\n");
    for(int i = 0; i < KYBER_N; i++) {
        printf("%4d ", s0_calc.coeffs[i]);
    }
    printf("\n");

    printf("s1_rec:\n");
    for(int i = 0; i < KYBER_N; i++) {
        printf("%4d ", s1_calc.coeffs[i]);
    }
    printf("\n");

    printf("s0:\n");
    for(int i = 0; i < KYBER_N; i++) {
        printf("%4d ", s.vec[0].coeffs[i]);
    }
    printf("\n");

    printf("s1:\n");
    for(int i = 0; i < KYBER_N; i++) {
        printf("%4d ", s.vec[1].coeffs[i]);
    }
    printf("\n");
    printf("\n");

    int correct = 1;
    for(int i = 0; i < KYBER_N; i++) {
        if(s0_calc.coeffs[i] != s.vec[0].coeffs[i]){
            correct = 0;
        }
    }
    if(correct == 1)
        printf("s0 correctly recovered!\n");
    else
        printf("s0 recover incorrect.\n");

    correct = 1;
    for(int i = 0; i < KYBER_N; i++) {
        if(s1_calc.coeffs[i] != s.vec[1].coeffs[i]){
            correct = 0;
        }
    }
    if(correct == 1)
        printf("s1 correctly recovered!\n");
    else
        printf("s1 recover incorrect.\n");

    return 0;
}






