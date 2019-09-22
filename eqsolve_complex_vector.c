#ifndef _EQSOLVE_COMPLEX_VECTOR_C_
#define _EQSOLVE_COMPLEX_VECTOR_C_
#include "eqsolve_complex_vector.h"
#include <stdlib.h>

eqsolve_vector_c *eqsolve_vector_c_alloc(size_t N) {
  eqsolve_vector_c *object = malloc(sizeof(eqsolve_vector_c));
  object->size = N;
  object->vector = malloc(N * sizeof(double complex));
  return object;
};

int eqsolve_vector_c_set(eqsolve_vector_c *vector, int i,
                         double complex value) {
  if (i >= 0 && i < vector->size)
    vector->vector[i] = value;
  return 0;
};

int eqsolve_vector_c_set_all(eqsolve_vector_c *vector, double complex value) {
  if (vector->size > 0) {
    for (int i = 0; i < vector->size; i++) {
      vector->vector[i] = value;
    }
  }
  return 0;
};

int eqsolve_vector_c_free(eqsolve_vector_c *vector) {
  free(vector->vector);
  free(vector);
  return EXIT_SUCCESS;
};

#endif