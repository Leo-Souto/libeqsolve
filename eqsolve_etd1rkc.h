#ifndef _EQSOLVE_PDE_ETD1RKC_H_
#define _EQSOLVE_PDE_ETD1RKC_H_
#include "complex.h"
#include "eqsolve_complex_vector.h"
#include "math.h"

#ifndef pi
#define pi M_PI
#endif
typedef struct eqsolve_etd1rkc {
  double dt;
  double time;
  int time_it;
  int numel_U;
  int numel_M;
  void *params;

  eqsolve_vector_c *U;
  eqsolve_vector_c *U_source;
  eqsolve_vector_c *N_source;
  eqsolve_vector_c *L;
  eqsolve_vector_c *N;
  eqsolve_vector_c *Ech;
  eqsolve_vector_c *Echm1;

  int (*sf)(eqsolve_vector_c *source, eqsolve_vector_c *destination,
            void *sf_params, double time_1);
  int (*nf)(eqsolve_vector_c *source, eqsolve_vector_c *destination,
            void *nf_params, double time_2);
  int (*ff)(eqsolve_vector_c *source, eqsolve_vector_c *destination,
            void *ff_params, double time_3);
} eqsolve_etd1rkc;

eqsolve_etd1rkc *eqsolve_etd1rkc_alloc(int n);
int eqsolve_etd1rkc_init(
    eqsolve_etd1rkc *solver, eqsolve_vector_c *U_source_,
    eqsolve_vector_c *N_source_, eqsolve_vector_c *L_, void *params_,
    double time_step,
    int (*sf_)(eqsolve_vector_c *U_source_1, eqsolve_vector_c *N_source_1,
               void *sf_params, double time_1),
    int (*nf_)(eqsolve_vector_c *U_source_2, eqsolve_vector_c *N_source_2,
               void *nf_params, double time_2),
    int (*ff_)(eqsolve_vector_c *U_source_3, eqsolve_vector_c *N_source_3,
               void *ff_params, double time_3));
int eqsolve_etd1rkc_get_solution(eqsolve_etd1rkc *solver,
                                 eqsolve_vector_c *destination);
int eqsolve_etd1rkc_compute_terms(eqsolve_etd1rkc *solver);

int eqsolve_etd1rkc_free(eqsolve_etd1rkc *solver);
int eqsolve_etd1rkc_update(eqsolve_etd1rkc *solver);

#include "eqsolve_etd1rkc.c"
#endif