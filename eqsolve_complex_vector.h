#ifndef _EQSOLVE_COMPLEX_VECTOR_H_
#define _EQSOLVE_COMPLEX_VECTOR_H_
#include "complex.h"
#include "fftw3.h"

typedef struct eqsolve_vector_c {
  size_t size;
  double complex *vector;
  size_t stride;
} eqsolve_vector_c;
extern const struct eqsolve_vector_c_default eqsolve_vector_c_default;

eqsolve_vector_c *eqsolve_vector_c_alloc(size_t N);
int eqsolve_vector_c_set(eqsolve_vector_c *vector, int i, double complex value);
int eqsolve_vector_c_set_all(eqsolve_vector_c *vector, double complex value);
int eqsolve_vector_c_free(eqsolve_vector_c *vector);

#include "eqsolve_complex_vector.c"
#endif