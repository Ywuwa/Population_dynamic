#include "mfunctions.h"
# define M_PI           3.14159265358979323846  /* pi */
//================================= CONSUMPTION RATE FUNC =========================================
double alpha(const double* x, const double* y, const float sig, const int dim)
{
    double res = 0.0;
    for (int i=0; i<dim; i++){
        const double buf = x[i] - y[i];
        res -= pow(buf, 2);
    }
    double norm = 2*pow(sig, 2);
    res = exp(res/norm);
    norm = pow(norm*M_PI, dim/2.0);
    res /= norm;
    return res;
}
//================================ ATTACK EFFICIENCY FUNC =========================================
double alpha_m(const double* x, const double* y, const float* m, const float sig, const int dim)
{
    double res = 0.0;
    for (int i=0; i<dim; i++){
        const double buf = x[i] - (y[i]+m[i]);
        res -= pow(buf, 2);
    }
    double norm = 2*pow(sig, 2);
    res = exp(res/norm);
    norm = pow(norm*M_PI, dim/2.0);
    res /= norm;
    return res;
}
//=============================== MODEL CARRYING CAPACITY =========================================
double Kore(const double* x, const float sig, const int k0, const int dim)
{
    double res = 0.0;
    for (int i=0; i<dim; i++){
        res -= pow(x[i], 4);
    }
    double norm = 4.0; // until sig_k=1, otherwise norm=4*pow(sig,4);
    res = exp(res/norm);
    res *= k0;
    return res;
}
//====================================== INTEGRATION =============================================
double integrate(const double *R, const double *x, const float sig_a, const float eps, const int dim)
{
    double res = 0.0;
    std::unique_ptr<double []> y = std::make_unique<double []> (dim);
    const double b = 3.0;
    unsigned int n = 1;
    for (int i=0; i<dim; i++) n *= 60;
    switch(dim)
    {
        case 1:
            {
            const double dx = 0.1;
            for (size_t i=0; i<n; i++){
                y[0] = -b + i*dx;
                const double f = alpha(x, y.get(), sig_a, dim);
                res += R[i]*eps*f;
            }
            res *= dx;
            break;
            }
        case 2:
            {
            const double dx1 = 0.1;
            const double dx2 = 0.1;
            for (unsigned int j=0; j<n; j++){
                y[0] = -b + (j/60)*dx1;
                y[1] = -b + (j%60)*dx2;
                const double f = alpha(x, y.get(), sig_a, dim);
                res += R[j]*f;
            }
            res *= dx1*dx2*eps;
            break;
            }
        default:
            break;
    }
    return res;
}
//=================================== RESOURSE DENSITY ============================================
double res_dense(const double* y, const std::vector<population> &poplist, const double R, Model &model, size_t l)
{
    double res = 0.0;
    double vec_sum = 0.0;
    for (size_t i=0; i<l; i++){
        const double dist = alpha(poplist[i].x.get(), y, model.sig_a, model.dim);
        vec_sum += poplist[i].r*poplist[i].N*dist;
    }
    const float sig_k = 1.0;
    const double k = Kore(y, sig_k, model.k0, model.dim);
    if (model.rtype==1) res = R * (model.birth - (model.birth*R)/k - vec_sum);
    else res = k - R * (vec_sum + model.dissip);
    return res;
}
//================================= RUNGE KUTTA RESOURSE ===========================================
double runge_compute_res(const double* y, const std::vector<population> &poplist, const double R, Model &model, size_t l)
{
    const double T = 0.05; // time step for runge_kutta
    const double k1 = res_dense(y, poplist, R, model, l);
    const double k2 = res_dense(y, poplist, R+T*k1/2.0, model, l);
    const double k3 = res_dense(y, poplist, R+T*k2/2.0, model, l);
    const double k4 = res_dense(y, poplist, R+T*k3, model, l);

    double res = R + (T/6.0) * (k1 + 2*k2 + 2*k3 + k4);
    return res;
}
//================================== POPULATION DENSITY ============================================
double pop_dense(const std::vector<population> &poplist, population &cur_p,
                 const double *R, const float *m, const double cur_N,
                 const float xi, const float delta, const float eps, const float sig_a, const float sig_g,
                 const int dim, size_t l)
{
    double vec_sum = 0.0;
    for (size_t i=0; i<l; i++){
        const double buf1 = alpha_m(cur_p.x.get(), poplist[i].x.get(), m, sig_g, dim);
        const double buf2 = alpha_m(poplist[i].x.get(), cur_p.x.get(), m, sig_g, dim);
        vec_sum += poplist[i].N * (cur_p.p*xi*buf1 - poplist[i].p*buf2);
    }
    const double buf = cur_p.r * integrate(R, cur_p.x.get(), sig_a, eps, dim) - delta;
    double res = cur_p.N * (buf + vec_sum);

    return res;
}
//================================ RUNGE KUTTA POPULATION ===========================================
double runge_compute_pop(const std::vector<population> &poplist, population &cur_p,
                 const double *R, const float *m, const double cur_N,
                 const float xi, const float delta, const float eps, const float sig_a, const float sig_g,
                 const int dim, size_t l)
{
    const double T = 0.05; // time step for runge_kutta
    const double k1 = pop_dense(poplist, cur_p,
                                    R, m, cur_N,
                                    xi, delta, eps, sig_a, sig_g,
                                    dim, l);
    const double k2 = pop_dense(poplist, cur_p,
                                    R, m, cur_N+T*k1/2.0,
                                    xi, delta, eps, sig_a, sig_g, dim, l);
    const double k3 = pop_dense(poplist, cur_p,
                                    R, m, cur_N+T*k2/2.0,
                                    xi, delta, eps, sig_a, sig_g, dim, l);
    const double k4 = pop_dense(poplist, cur_p,
                                    R, m, cur_N+T*k3,
                                    xi, delta, eps, sig_a, sig_g, dim, l);

    double res = cur_N + (T/6.0) * (k1 + 2*k2 + 2*k3 + k4);
    //printf("res: %10.5e, R: %10.5e, k1: %10.5e, k2: %10.5e, k3: %10.5e, k4: %10.5e\n", res, R, k1, k2, k3, k4);
    return res;
}
//================================ PARALLEL RUNGE KUTTA POPULATION ===========================================
void paral_runge_pop(std::vector<population> &poplist,
                 const double *R, const float *m,
                 const float xi, const float delta, const float eps, const float sig_a, const float sig_g,
                 const int dim, size_t l)
{
    std::unique_ptr<double []> k1 = std::make_unique<double []> (poplist.size());
    std::unique_ptr<double []> k2 = std::make_unique<double []> (poplist.size());
    std::unique_ptr<double []> k3 = std::make_unique<double []> (poplist.size());
    std::unique_ptr<double []> k4 = std::make_unique<double []> (poplist.size());
    std::unique_ptr<double []> N = std::make_unique<double []> (poplist.size());
    const double T = 0.05; // time step for runge_kutta
    int num = omp_get_num_procs();
    omp_set_num_threads(num);
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++)
        {
            N[i] = poplist[i].N;
            k1[i] = pop_dense(poplist, poplist[i], R, m, poplist[i].N, xi, delta, eps, sig_a, sig_g, dim, l);
        }
        #pragma omp barrier
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++)
        {
            poplist[i].N = N[i] + T*k1[i]*0.5;
        }
        #pragma omp barrier
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++)
        {
            k2[i] = pop_dense(poplist, poplist[i], R, m, poplist[i].N, xi, delta, eps, sig_a, sig_g, dim, l);
        }
        #pragma omp barrier
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++)
        {
            poplist[i].N = N[i] + T*k2[i]*0.5;
        }
        #pragma omp barrier
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++)
        {
            k3[i] = pop_dense(poplist, poplist[i], R, m, poplist[i].N, xi, delta, eps, sig_a, sig_g, dim, l);
        }
        #pragma omp barrier
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++)
        {
            poplist[i].N = N[i] + T*k3[i];
        }
        #pragma omp barrier
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++)
        {
            k4[i] = pop_dense(poplist, poplist[i], R, m, poplist[i].N, xi, delta, eps, sig_a, sig_g, dim, l);
        }
        #pragma omp barrier
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++)
        {
            poplist[i].N = N[i] + (T/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        }
    }
}
