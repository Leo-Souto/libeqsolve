#ifndef _EQSOLVE_PDE_ETD4RKC_H_
#define _EQSOLVE_PDE_ETD4RKC_H_
#include "complex.h"
#include "eqsolve_complex_vector.h"
#include "math.h"

#ifndef pi
#define pi M_PI
#endif
typedef struct eqsolve_etd4rkc {
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
  eqsolve_vector_c *A;
  eqsolve_vector_c *B;
  eqsolve_vector_c *C;
  eqsolve_vector_c *N;
  eqsolve_vector_c *Ech;
  eqsolve_vector_c *Ech2;
  eqsolve_vector_c *Echm1;
  eqsolve_vector_c *Ech2m1;
  eqsolve_vector_c *T1;
  eqsolve_vector_c *T2;
  eqsolve_vector_c *T3;
  eqsolve_vector_c *U_N;
  eqsolve_vector_c *A_N;
  eqsolve_vector_c *B_N;

  int (*sf)(eqsolve_vector_c *source, eqsolve_vector_c *destination,
            void *sf_params, double time);
  int (*nf)(eqsolve_vector_c *source, eqsolve_vector_c *destination,
            void *nf_params, double time);
  int (*ff)(eqsolve_vector_c *source, eqsolve_vector_c *destination,
            void *ff_params, double time);
} eqsolve_etd4rkc;

eqsolve_etd4rkc *eqsolve_etd4rkc_alloc(int n);
int eqsolve_etd4rkc_init(
    eqsolve_etd4rkc *solver, eqsolve_vector_c *U_source_,
    eqsolve_vector_c *N_source_, eqsolve_vector_c *L_, void *params_,
    double time_step,
    int (*sf_)(eqsolve_vector_c *U_source_1, eqsolve_vector_c *N_source_1,
               void *sf_params, double time),
    int (*nf_)(eqsolve_vector_c *U_source_2, eqsolve_vector_c *N_source_2,
               void *nf_params, double time),
    int (*ff_)(eqsolve_vector_c *U_source_3, eqsolve_vector_c *N_source_3,
               void *ff_params, double time));

int eqsolve_etd4rkc_get_solution(eqsolve_etd4rkc *solver,
                                 eqsolve_vector_c *destination);
int eqsolve_etd4rkc_compute_terms(eqsolve_etd4rkc *solver);

int eqsolve_etd4rkc_free(eqsolve_etd4rkc *solver);
int eqsolve_etd4rkc_update(eqsolve_etd4rkc *solver);

#include "eqsolve_etd4rkc.c"
#endif