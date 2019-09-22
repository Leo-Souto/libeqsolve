#ifndef _EQSOLVE_PDE_ETD1RKC_C_
#define _EQSOLVE_PDE_ETD1RKC_C_

#include "eqsolve_etd1rkc.h"
#include "complex.h"
#include "math.h"
#include "stdlib.h"
#ifndef pi
#define pi M_PI
#endif

/*Allocates memory for complex ETD1RK.*/
eqsolve_etd1rkc *eqsolve_etd1rkc_alloc(int N) {
  eqsolve_etd1rkc *object = malloc(sizeof(eqsolve_etd1rkc));
  object->U = eqsolve_vector_c_alloc(N);
  object->L = eqsolve_vector_c_alloc(N);
  object->N = eqsolve_vector_c_alloc(N);
  object->Ech = eqsolve_vector_c_alloc(N);
  object->Echm1 = eqsolve_vector_c_alloc(N);
  object->numel_U = N;
  object->numel_M = 32;
  object->time = 0;
  object->time_it = 0;

  return object;
};

/*Initiate the ETD1RK complex solver and compute the terms.*/
int eqsolve_etd1rkc_init(
    eqsolve_etd1rkc *solver, eqsolve_vector_c *U_source_,
    eqsolve_vector_c *N_source_, eqsolve_vector_c *L_, void *params_,
    double time_step,
    int (*sf_)(eqsolve_vector_c *U_source_1, eqsolve_vector_c *N_source_1,
               void *sf_params, double time_1),
    int (*nf_)(eqsolve_vector_c *U_source_2, eqsolve_vector_c *N_source_2,
               void *nf_params, double time_2),
    int (*ff_)(eqsolve_vector_c *U_source_3, eqsolve_vector_c *N_source_3,
               void *ff_params, double time_3)) {
  if (U_source_->size == solver->numel_U && solver->numel_U != 0) {
    if (L_->size == U_source_->size) {
      for (int i = 0; i < solver->numel_U; i++) {
        solver->U->vector[i] = U_source_->vector[i];
        solver->L->vector[i] = L_->vector[i];
      }
      solver->dt = time_step;
      solver->U_source = U_source_;
      solver->N_source = N_source_;
      solver->sf = sf_;
      solver->nf = nf_;
      solver->ff = ff_;
      solver->params = params_;
      eqsolve_etd1rkc_compute_terms(solver);
    } else {
      perror("eqsolve_vector_c sizes not compatible\n");
      exit(EXIT_FAILURE);
    }
  } else {
    perror("U_source_ size different from solver size, or solver not allocated "
           "properly\n");
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
};

/*Get the current complex ETD1RK solution.*/
int eqsolve_etd1rkc_get_solution(eqsolve_etd1rkc *solver,
                                 eqsolve_vector_c *destination) {
  if (destination->size == solver->numel_U && solver->numel_U != 0) {
    for (int i = 0; i < destination->size; i++)
      destination->vector[i] = solver->U->vector[i];
  } else {
    perror("not compatible sizes between solver and source\n");
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
};

/*Set the current complex ETD1RK solution.*/
int eqsolve_etd1rkc_set_solution(eqsolve_etd1rkc *solver,
                                 eqsolve_vector_c *source) {
  if (source->size == solver->numel_U && solver->numel_U != 0) {
    for (int i = 0; i < source->size; i++)
      solver->U->vector[i] = source->vector[i];
  } else {
    perror("not compatible sizes between solver and source\n");
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
};
/*Compute the omplex ETD1RK constant terms.*/
int eqsolve_etd1rkc_compute_terms(eqsolve_etd1rkc *solver) {
  if (solver->numel_U != 0) {
    eqsolve_vector_c *pts = eqsolve_vector_c_alloc(solver->numel_M);
    eqsolve_vector_c *LR =
        eqsolve_vector_c_alloc((solver->numel_M) * (solver->numel_U));
    for (int i = 0; i < solver->numel_M; i++) {
      pts->vector[i] = cexp(I * pi * ((i + 1) - 0.5) / solver->numel_U);
    }
    for (int i = 0; i < solver->numel_U; i++) {
      for (int j = 0; j < solver->numel_M; j++) {
        LR->vector[j + solver->numel_M * i] =
            solver->dt * solver->L->vector[i] + pts->vector[j];
      }
    }
    for (int i = 0; i < solver->numel_U; i++) {
      solver->Ech->vector[i] = cexp(solver->L->vector[i] * solver->dt);
    }

    double complex sum;
    for (int i = 0; i < solver->numel_U; i++) {
      sum = 0;
      for (int j = 0; j < solver->numel_M; j++) {
        sum += (cexp(LR->vector[j + solver->numel_M * i]) - 1.0) /
               LR->vector[j + solver->numel_M * i];
      }
      solver->Echm1->vector[i] = solver->dt * creal(sum) / solver->numel_M;
    }
    eqsolve_vector_c_free(pts);
    eqsolve_vector_c_free(LR);
  } else {
    perror("solver not allocated\n");
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
};

/*Compute the next iteration of the solution.*/
int eqsolve_etd1rkc_update(eqsolve_etd1rkc *solver) {
  if (solver->sf != NULL) {
    solver->sf(solver->U_source, solver->N_source, solver->params,
               solver->time);
    for (int i = 0; i < solver->numel_U; i++)
      solver->U->vector[i] = solver->U_source->vector[i];
  }
  solver->nf(solver->U_source, solver->N_source, solver->params, solver->time);
  for (int i = 0; i < solver->numel_U; i++) {
    solver->U->vector[i] = solver->U_source->vector[i] =
        solver->U->vector[i] * solver->Ech->vector[i] +
        solver->Echm1->vector[i] * solver->N_source->vector[i];
  }
  if (solver->ff != NULL) {
    solver->ff(solver->U_source, solver->N_source, solver->params,
               solver->time);
    for (int i = 0; i < solver->numel_U; i++)
      solver->U->vector[i] = solver->U_source->vector[i];
  }
  solver->time_it++;
  solver->time = (solver->time_it) * solver->dt;
  return EXIT_SUCCESS;
};
/*Free the memory allocated for ETD1RK solver.*/
int eqsolve_etd1rkc_free(eqsolve_etd1rkc *solver) {
  eqsolve_vector_c_free(solver->U);
  eqsolve_vector_c_free(solver->L);
  eqsolve_vector_c_free(solver->N);
  eqsolve_vector_c_free(solver->Ech);
  eqsolve_vector_c_free(solver->Echm1);

  free(solver);

  return EXIT_SUCCESS;
};

#endif