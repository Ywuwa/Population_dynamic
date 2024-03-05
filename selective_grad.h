#include "mfunctions.h"
# define M_PI           3.14159265358979323846  /* pi */
double selective_grad(Model &model, std::vector<population> &poplist);
double len(const double *v, const int dim);

double ifunc(const double *x, const double *y, const float sig_a, const int dim, const int pos);
typedef double(* integrateFunc)(const double *x, const double *y, const float sig_a, const int dim, const int pos);
double integrate_derivative(const double *R, const double *x, const float sig_a, const float eps, const int dim, const int pos, integrateFunc func);
