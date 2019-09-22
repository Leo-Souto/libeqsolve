#ifndef _EQSOLVE_VECTOR_C_
#define _EQSOLVE_VECTOR_C_
#include "eqsolve_vector.h"
#include <stdlib.h>

eqsolve_vector *eqsolve_vector_alloc(size_t N) {
  eqsolve_vector *object = malloc(sizeof(eqsolve_vector));
  object->size = N;
  object->vector = malloc(N * sizeof(double));
  return object;
};

int eqsolve_vector_set(eqsolve_vector *vector, int i, double value) {
  if (i >= 0 && i < vector->size)
    vector->vector[i] = value;
  return 0;
};

int eqsolve_vector_set_all(eqsolve_vector *vector, double value) {
  if (vector->size > 0) {
    for (int i = 0; i < vector->size; i++) {
      vector->vector[i] = value;
    }
  }
  return 0;
};

int eqsolve_vector_free(eqsolve_vector *vector) {
  free(vector->vector);
  free(vector);
  return EXIT_SUCCESS;
};

double eqsolve_vector_abs_max(eqsolve_vector *vector) {
  double max = 0;
  for (int i = 0; i < vector->size; i++)
    if (vector->vector[i] > max)
      max = vector->vector[i];
  return max;
};

#endif