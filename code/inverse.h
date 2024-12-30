#ifndef INVERSE_H
#define INVERSE_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define KYBER_Q 3329
#define POLYBYTES 256


typedef struct{
  int16_t coeffs[POLYBYTES];
} poly;

typedef struct{
  int16_t coeffs[POLYBYTES+1];
} dpoly;

void inverse(poly* r, poly* a);

#endif
