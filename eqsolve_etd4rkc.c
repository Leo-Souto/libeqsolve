#ifndef _EQSOLVE_PDE_ETD4RKC_C_
#define _EQSOLVE_PDE_ETD4RKC_C_

#include "eqsolve_etd4rkc.h"
#include "complex.h"
#include "math.h"
#include "stdlib.h"
#include <omp.h>
#ifndef pi
#define pi M_PI
#endif

/*Allocates memory for complex ETD4RK.*/
eqsolve_etd4rkc *eqsolve_etd4rkc_alloc(int N) {
  eqsolve_etd4rkc *object = malloc(sizeof(eqsolve_etd4rkc));
  object->U = eqsolve_vector_c_alloc(N);
  object->L = eqsolve_vector_c_alloc(N);
  object->N = eqsolve_vector_c_alloc(N);
  object->A = eqsolve_vector_c_alloc(N);
  object->B = eqsolve_vector_c_alloc(N);
  object->C = eqsolve_vector_c_alloc(N);
  object->T1 = eqsolve_vector_c_alloc(N);
  object->T2 = eqsolve_vector_c_alloc(N);
  object->T3 = eqsolve_vector_c_alloc(N);
  object->Ech = eqsolve_vector_c_alloc(N);
  object->Ech2 = eqsolve_vector_c_alloc(N);
  object->Echm1 = eqsolve_vector_c_alloc(N);
  object->Ech2m1 = eqsolve_vector_c_alloc(N);
  object->U_N = eqsolve_vector_c_alloc(N);
  object->A_N = eqsolve_vector_c_alloc(N);
  object->B_N = eqsolve_vector_c_alloc(N);
  object->numel_U = N;
  object->numel_M = 32;
  return object;
};

/*Initiate the ETD4RK complex solver and compute the terms.*/
int eqsolve_etd4rkc_init(
    eqsolve_etd4rkc *solver, eqsolve_vector_c *U_source_,
    eqsolve_vector_c *N_source_, eqsolve_vector_c *L_, void *params_,
    double time_step,
    int (*sf_)(eqsolve_vector_c *U_source_1, eqsolve_vector_c *N_source_1,
               void *sf_params, double time),
    int (*nf_)(eqsolve_vector_c *U_source_2, eqsolve_vector_c *N_source_2,
               void *nf_params, double time),
    int (*ff_)(eqsolve_vector_c *U_source_3, eqsolve_vector_c *N_source_3,
               void *ff_params, double time)) {
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
      solver->time = 0;
      solver->time_it = 0;
      eqsolve_etd4rkc_compute_terms(solver);
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

/*Get the current complex ETD4RK solution.*/
int eqsolve_etd4rkc_get_solution(eqsolve_etd4rkc *solver,
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

/*Set the current complex ETD4RK solution.*/
int eqsolve_etd4rkc_set_solution(eqsolve_etd4rkc *solver,
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
/*Compute the complex ETD4RK constant terms.*/
int eqsolve_etd4rkc_compute_terms(eqsolve_etd4rkc *solver) {
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
      solver->Ech2->vector[i] = cexp(solver->L->vector[i] * solver->dt / 2.0);
    }

    double complex sum, sum2;
    for (int i = 0; i < solver->numel_U; i++) {
      sum = 0;
      sum2 = 0;
      for (int j = 0; j < solver->numel_M; j++) {
        sum += (cexp(LR->vector[j + solver->numel_M * i]) - 1.0) /
               LR->vector[j + solver->numel_M * i];
        sum2 += (cexp(LR->vector[j + solver->numel_M * i] / 2.0) - 1.0) /
                LR->vector[j + solver->numel_M * i];
      }
      solver->Echm1->vector[i] = solver->dt * creal(sum) / solver->numel_M;
      solver->Ech2m1->vector[i] = solver->dt * creal(sum2) / solver->numel_M;
    }

    for (int i = 0; i < solver->numel_U; i++) {
      sum = 0;
      for (int j = 0; j < solver->numel_M; j++) {
        // [−4 − hc + e^ch(4 − 3hc + h^2c^2)]
        sum += (-4 - LR->vector[j + solver->numel_M * i] +
                cexp(LR->vector[j + solver->numel_M * i]) *
                    (4 - 3 * LR->vector[j + solver->numel_M * i] +
                     cpow(LR->vector[j + solver->numel_M * i], 2))) /
               cpow(LR->vector[j + solver->numel_M * i], 3);
      }
      solver->T1->vector[i] = solver->dt * creal(sum) / solver->numel_M;
    }
    for (int i = 0; i < solver->numel_U; i++) {
      sum = 0;
      for (int j = 0; j < solver->numel_M; j++) {
        //[2 + hc + e^ch(−2 + hc)]
        sum += (2 + LR->vector[j + solver->numel_M * i] +
                cexp(LR->vector[j + solver->numel_M * i]) *
                    (-2 + LR->vector[j + solver->numel_M * i])) /
               cpow(LR->vector[j + solver->numel_M * i], 3);
      }
      solver->T2->vector[i] = solver->dt * creal(sum) / solver->numel_M;
    }
    for (int i = 0; i < solver->numel_U; i++) {
      sum = 0;
      for (int j = 0; j < solver->numel_M; j++) {
        //[−4 − 3hc − h^2c^2 + e^ch(4 − hc)]
        sum += (-4 - 3 * LR->vector[j + solver->numel_M * i] -
                cpow(LR->vector[j + solver->numel_M * i], 2) +
                cexp(LR->vector[j + solver->numel_M * i]) *
                    (4 - LR->vector[j + solver->numel_M * i])) /
               cpow(LR->vector[j + solver->numel_M * i], 3);
      }
      solver->T3->vector[i] = solver->dt * creal(sum) / solver->numel_M;
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
int eqsolve_etd4rkc_update(eqsolve_etd4rkc *solver) {
  if (solver->sf != NULL) {
    solver->sf(solver->U_source, solver->N_source, solver->params,
               solver->time);
    for (int i = 0; i < solver->numel_U; i++)
      solver->U->vector[i] = solver->U_source->vector[i];
  }
  solver->nf(solver->U_source, solver->N_source, solver->params, solver->time);
#pragma omp parallel for
  for (int i = 0; i < solver->numel_U; i++) {
    // Compute A(u);
    solver->A->vector[i] = solver->U_source->vector[i] =
        solver->U->vector[i] * solver->Ech2->vector[i] +
        solver->Ech2m1->vector[i] * solver->N_source->vector[i];
    // copies N(u) to another vector
    solver->U_N->vector[i] = solver->N_source->vector[i];
  }

  solver->nf(solver->U_source, solver->N_source, solver->params,
             (solver->time_it + 0.5) * solver->dt);
#pragma omp parallel for
  for (int i = 0; i < solver->numel_U; i++) {
    solver->B->vector[i] = solver->U_source->vector[i] =
        solver->U->vector[i] * solver->Ech2->vector[i] +
        solver->Ech2m1->vector[i] * solver->N_source->vector[i];
    solver->A_N->vector[i] = solver->N_source->vector[i];
  }
  solver->nf(solver->U_source, solver->N_source, solver->params,
             (solver->time_it + 0.5) * solver->dt);
#pragma omp parallel for
  for (int i = 0; i < solver->numel_U; i++) {
    solver->C->vector[i] = solver->U_source->vector[i] =
        solver->A->vector[i] * solver->Ech2->vector[i] +
        solver->Ech2m1->vector[i] *
            (2 * solver->N_source->vector[i] - solver->U_N->vector[i]);
    solver->B_N->vector[i] = solver->N_source->vector[i];
  }

  solver->time_it++;
  solver->time = (solver->time_it) * solver->dt;
  solver->nf(solver->U_source, solver->N_source, solver->params, solver->time);
#pragma omp parallel for
  for (int i = 0; i < solver->numel_U; i++) {
    solver->U->vector[i] = solver->U_source->vector[i] =
        solver->U->vector[i] * solver->Ech->vector[i] +
        solver->U_N->vector[i] * solver->T1->vector[i] +
        2 * (solver->A_N->vector[i] + solver->B_N->vector[i]) *
            solver->T2->vector[i] +
        solver->N_source->vector[i] * solver->T3->vector[i];
  }
  if (solver->ff != NULL) {
    solver->ff(solver->U_source, solver->N_source, solver->params,
               solver->time);
    for (int i = 0; i < solver->numel_U; i++)
      solver->U->vector[i] = solver->U_source->vector[i];
  }
  return EXIT_SUCCESS;
};

/*Free the memory allocated for ETD4RK solver.*/
int eqsolve_etd4rkc_free(eqsolve_etd4rkc *solver) {
  eqsolve_vector_c_free(solver->U);
  eqsolve_vector_c_free(solver->L);
  eqsolve_vector_c_free(solver->N);
  eqsolve_vector_c_free(solver->A);
  eqsolve_vector_c_free(solver->B);
  eqsolve_vector_c_free(solver->C);
  eqsolve_vector_c_free(solver->T1);
  eqsolve_vector_c_free(solver->T2);
  eqsolve_vector_c_free(solver->T3);
  eqsolve_vector_c_free(solver->Ech);
  eqsolve_vector_c_free(solver->Ech2);
  eqsolve_vector_c_free(solver->Echm1);
  eqsolve_vector_c_free(solver->Ech2m1);
  eqsolve_vector_c_free(solver->U_N);
  eqsolve_vector_c_free(solver->A_N);
  eqsolve_vector_c_free(solver->B_N);

  free(solver);

  return EXIT_SUCCESS;
};

#endif