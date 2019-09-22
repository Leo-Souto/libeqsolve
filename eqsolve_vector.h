#ifndef _EQSOLVE__VECTOR_H_
#define _EQSOLVE__VECTOR_H_

typedef struct eqsolve_vector {
  size_t size;
  double *vector;
  size_t stride;
} eqsolve_vector;
extern const struct eqsolve_vector_default eqsolve_vector_default;

eqsolve_vector *eqsolve_vector_alloc(size_t N);
int eqsolve_vector_set(eqsolve_vector *vector, int i, double value);
int eqsolve_vector_set_all(eqsolve_vector *vector, double value);
int eqsolve_vector_free(eqsolve_vector *vector);
double eqsolve_vector_abs_max(eqsolve_vector *vector);

#include "eqsolve_vector.c"

#endif