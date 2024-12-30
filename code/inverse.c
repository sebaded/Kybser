#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "inverse.h"


int16_t find_quot(int num, int d){
    int i = 0;
    while(!((num + (KYBER_Q * i)) % d == 0))
        i++;
    return (int16_t)(((num + (KYBER_Q * i)) / d) % KYBER_Q) ;
}





int degree(dpoly num) {
    for(int i = POLYBYTES; i >= 0; i--)
        if(num.coeffs[i] != 0)
            return i;
    return 0;
}





void poly_div(dpoly qr[2], dpoly a, dpoly b){

    for(int i = 0; i <= POLYBYTES; i++){
        qr[0].coeffs[i] = 0;
        qr[1].coeffs[i] = 0;
    }

    while(degree(a) >= degree(b)) {
        int i = degree(a);
        int j = degree(b);

        qr[0].coeffs[i-j] = find_quot(a.coeffs[i], b.coeffs[j]);

        for(int k = j; k >= 0; k--) {
            a.coeffs[k+(i-j)] -= (b.coeffs[k] * qr[0].coeffs[i-j]) % KYBER_Q;
            a.coeffs[k+(i-j)] %= KYBER_Q;
        }
        if(i+j == 0)
            break;
    }
    qr[1] = a;
}





int is_empty(dpoly r) {
    int z = 1;
    for(int i = 0; i < POLYBYTES+1; i++){
        if(r.coeffs[i] != 0 )
            z = 0;
    }
    return z;
}





dpoly mult(dpoly a, dpoly b) {
    dpoly m;
    for(int i = 0; i < POLYBYTES+1; i++)
        m.coeffs[i] = 0;

    for(int i = 0; i < POLYBYTES+1; i++){
        for(int j = 0; j < POLYBYTES+1; j++){
            if((i+j) < POLYBYTES+1){
                m.coeffs[i+j] += ((a.coeffs[i] * b.coeffs[j]) % KYBER_Q);
                m.coeffs[i+j] %= KYBER_Q;
            }
        }
    }
    return m;
}




dpoly egcd(dpoly a, dpoly b, dpoly* s, dpoly* t) {
    dpoly qr[2];
    dpoly q;
    dpoly old_s;

    if(is_empty(b)){
        s->coeffs[0] = 1;
        t->coeffs[0] = 0;
        return a;
    }
    else {
        poly_div(qr, a, b);
        memcpy(&q, &qr[0], sizeof(dpoly));
        dpoly d = egcd(b, qr[1], s, t);
        memcpy(&old_s, s, sizeof(dpoly));
        memcpy(s, t, sizeof(dpoly));
        dpoly qs = mult(q, *s);

        for(int i = 0; i < POLYBYTES+1; i++) {
            t->coeffs[i] = (old_s.coeffs[i] - qs.coeffs[i]);
            t->coeffs[i] %= KYBER_Q;
        }
        return d;
    }
}

void inverse(poly* r, poly* a) {
    dpoly a1, mod, s, t;
    memcpy(&a1, a, sizeof(poly));
    a1.coeffs[POLYBYTES] = 0;

    for(int i = 0; i <= POLYBYTES; i++){
        mod.coeffs[i] = 0;
        s.coeffs[i] = 0;
        t.coeffs[i] = 0;
    }

    //mod = x^256 + 1
    mod.coeffs[POLYBYTES] = 1;
    mod.coeffs[0] = 1;

    dpoly d = egcd(mod, a1, &s, &t);

    if(degree(d) == 0 && d.coeffs[0] != 0) {
        for(int i = 0; i < POLYBYTES; i++)
            r->coeffs[i] = find_quot(t.coeffs[i], d.coeffs[0]) % KYBER_Q;
    }
    else {
        printf("Couldn't find inverse. Try again!\n");
    }
}

