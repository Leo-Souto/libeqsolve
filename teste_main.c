#include "complex.h"
#include "eqsolve_complex_vector.h"
#include "eqsolve_etd4rkc.h"
#include "eqsolve_vector.h"

#include "fftw3.h"
#include "math.h"
#include "stdio.h"
#include "time.h"
#define LL 100.0
#define N 1000
#define Nk (N / 2 + 1)
#define alpha 0.0

struct params {
  int N_;
  fftw_plan Ur2c;
  fftw_plan Uc2r;
  fftw_plan U2r2c;
  fftw_plan U2c2r;
  eqsolve_vector *Un;
  eqsolve_vector *Un2;
  eqsolve_vector_c *Unk;
  eqsolve_vector_c *Un2k;
} params;

int nonlinear(eqsolve_vector_c *U_source, eqsolve_vector_c *N_source,
              void *params, double time) {

  fftw_execute(((struct params *)params)->Uc2r);

  for (int i = 0; i < ((struct params *)params)->Un->size; i++)
    ((struct params *)params)->Un2->vector[i] =
        ((struct params *)params)->Un->vector[i] *
        ((struct params *)params)->Un->vector[i] / N / N;

  fftw_execute(((struct params *)params)->U2r2c);

  for (int i = 0; i < N_source->size; i++)
    N_source->vector[i] =
        I * pi * i / LL * ((struct params *)params)->Un2k->vector[i];
  return EXIT_SUCCESS;
};

int main(int argc, char const *argv[]) {

  struct params *parameters = malloc(sizeof(struct params));
  parameters->Un = eqsolve_vector_alloc(N);
  parameters->Un2 = eqsolve_vector_alloc(N);
  parameters->Unk = eqsolve_vector_c_alloc(Nk);
  parameters->Un2k = eqsolve_vector_c_alloc(Nk);
  parameters->Ur2c = fftw_plan_dft_r2c_1d(
      N, parameters->Un->vector, parameters->Unk->vector, FFTW_PATIENT);
  parameters->Uc2r = fftw_plan_dft_c2r_1d(N, parameters->Unk->vector,
                                          parameters->Un->vector, FFTW_PATIENT);
  parameters->U2r2c = fftw_plan_dft_r2c_1d(
      N, parameters->Un2->vector, parameters->Un2k->vector, FFTW_PATIENT);
  parameters->U2c2r = fftw_plan_dft_c2r_1d(
      N, parameters->Un2k->vector, parameters->Un2->vector, FFTW_PATIENT);
  parameters->N_ = N;
  srand(time(NULL));
  eqsolve_vector_c *L_ = eqsolve_vector_c_alloc(Nk);
  for (int i = 0; i < Nk; i++) {
    L_->vector[i] =
        cpow(2.0 * pi * i / LL, 2) - cpow(2.0 * pi * i / LL, 4) - alpha;
  }

  for (int i = 0; i < N; i++)
    parameters->Un->vector[i] = 2.0 * rand() / RAND_MAX - 1.0;
  fftw_execute(parameters->Ur2c);

  eqsolve_etd4rkc *solver = eqsolve_etd4rkc_alloc(Nk);

  eqsolve_etd4rkc_init(solver, parameters->Unk, parameters->Un2k, L_,
                       parameters, 0.001, NULL, nonlinear, NULL);

  for (int i = 0; i < 1000; i++) {
    for (int j = 0; j < 500; j++)
      eqsolve_etd4rkc_update(solver);
    fftw_execute(parameters->Uc2r);
    double max = eqsolve_vector_abs_max(parameters->Un);
    for (int j = 0; j < N; j++) {
      printf("%f\n", parameters->Un->vector[j] / max);
    }
  }
  eqsolve_etd4rkc_free(solver);
  return 0;
}
