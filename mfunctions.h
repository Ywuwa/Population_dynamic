#pragma once
#include <string>
#include "omp.h"
#include "args.h"

double alpha(const double* x, const double* y, const float sig, const int dim);
double alpha_m(const double* x, const double* y, const float* m, const float sig, const int dim);
double Kore(const double* x, const float sig, const int k0, const int dim);
double integrate(const double *R, const double *x, const float sig_a, const float eps, const int dim);

double res_dense(const double* y, const std::vector<population> &poplist,
                 const double R, Model &model, size_t l);
double runge_compute_res(const double* y, const std::vector<population> &poplist,
                         const double R, Model &model, size_t l);

double pop_dense(const std::vector<population> &poplist, population &cur_p,
                 const double *R, const float *m, const double cur_N,
                 const float xi, const float delta, const float eps, const float sig_a, const float sig_g,
                 const int dim, size_t l);

double runge_compute_pop(const std::vector<population> &poplist, population &cur_p,
                 const double *R, const float *m, const double cur_N,
                 const float xi, const float delta, const float eps, const float sig_a, const float sig_g,
                 const int dim, size_t l);

void paral_runge_pop(std::vector<population> &poplist,
                     const double *R, const float *m,
                     const float xi, const float delta, const float eps, const float sig_a, const float sig_g,
                     const int dim, size_t l);
