#include "selective_grad.h"
//==================================== SELECTIVE GRADIENT ========================================

double selective_grad(Model &model, std::vector<population> &poplist)
{
    std::vector<double> sgrad; // compute gradient length for every population gradient, construct this vector
    std::unique_ptr<double []> v = std::make_unique<double[]>(model.dim);
    double x1 = 0.0;
    double x2 = 0.0;
    int num = omp_get_num_procs();
    omp_set_num_threads(num);
    #pragma omp parallel
    {
        std::unique_ptr<double []> grad = std::make_unique<double[]>(model.dim); // gradient for resident population
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++){
            for(int j=0; j<model.dim; j++){
                const double buf = poplist[i].r * integrate_derivative(model.R.get(), poplist[i].x.get(), model.sig_a, model.eps, model.dim, j, (integrateFunc)ifunc);
                double vec_sum = 0.0;
                for (size_t k=0; k<poplist.size(); k++){
                    const double mult1 = -(poplist[i].x[j]-poplist[k].x[j]-model.m[j]);
                    const double mult2 = (poplist[k].x[j]-poplist[i].x[j]-model.m[j]);
                    const double del = pow(model.sig_g, 2);
                    const double buf1 = mult1*alpha_m(poplist[i].x.get(), poplist[k].x.get(), model.m.get(), model.sig_g, model.dim)/del;
                    const double buf2 = mult2*alpha_m(poplist[k].x.get(), poplist[i].x.get(), model.m.get(), model.sig_g, model.dim)/del;
                    vec_sum += poplist[k].N * (poplist[i].p*model.xi*buf1 - poplist[k].p*buf2);
                }
                grad[j] = (buf + vec_sum) * poplist[i].N;
            }
            const double l = len(grad.get(), model.dim);
            #pragma omp critical
            {
                sgrad.push_back(l);
                x1 += grad[0];
                x2 += grad[1];
            }
        }
        /*v[0] = x1; v[1] = x2;
        const double l = len(v.get(), model.dim);
        sgrad.push_back(l); x1 = 0.0; x2 = 0.0;*/
    }
    double res = 0.0;
    for (size_t i=0; i<sgrad.size(); i++){
        res += sgrad[i];
    }
    /*std::unique_ptr<double []> v = std::make_unique<double[]>(model.dim);
    v[0] = x1; v[1] = x2;
    res = len(v.get(), model.dim);*/
    return res;
}
double len(const double *v, const int dim)
{
    double res = 0.0;
    for (int i=0; i<dim; i++){
        res += pow(v[i],2);
    }
    res = pow(res, 0.5);
    return res;
}
//================================== DERIVATIVE INTEGRATION =====================================
double ifunc(const double *x, const double *y, const float sig_a, const int dim, const int pos)
{
    const double a = alpha(x, y, sig_a, dim);
    const double norm = pow(sig_a, 2);
    double res = -(x[pos]-y[pos])*a/norm;

    return res;
}
double integrate_derivative(const double *R, const double *x, const float sig_a, const float eps, const int dim, const int pos, integrateFunc func)
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
                const double f = func(x, y.get(), sig_a, dim, pos);//alpha(x, y.get(), sig_a, dim);
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
                const double f = func(x, y.get(), sig_a, dim, pos);
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
//==================================== SELECTIVE GRADIENT BACKUP ========================================

double selective_grad_back(Model &model, std::vector<population> &poplist)
{
    std::vector<double> sgrad; // compute gradient length for every population gradient, construct this vector
    double x1 = 0.0;
    double x2 = 0.0;
    #pragma omp parallel
    {
        std::unique_ptr<double []> grad = std::make_unique<double[]>(model.dim); // gradient for resident population
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<poplist.size(); i++){
            for(int j=0; j<model.dim; j++){
                const double buf = poplist[i].r * integrate_derivative(model.R.get(), poplist[i].x.get(), model.sig_a, model.eps, model.dim, j, (integrateFunc)ifunc) - model.death;
                double vec_sum = 0.0;
                for (size_t k=0; k<poplist.size(); k++){
                    //if (k!=i){
                        const double mult1 = -(poplist[i].x[j]-poplist[k].x[j]-model.m[j]);
                        const double mult2 = (poplist[k].x[j]-poplist[i].x[j]-model.m[j]);
                        const double del = pow(model.sig_g, 2);
                        const double buf1 = mult1*alpha_m(poplist[i].x.get(), poplist[k].x.get(), model.m.get(), model.sig_g, model.dim)/del;
                        const double buf2 = mult2*alpha_m(poplist[k].x.get(), poplist[i].x.get(), model.m.get(), model.sig_g, model.dim)/del;
                        vec_sum += poplist[k].N * (poplist[i].p*model.xi*buf1 - poplist[k].p*buf2);
                    //}
                }
                grad[j] = buf + vec_sum;
            }
            const double l = len(grad.get(), model.dim);
            #pragma omp critical
            {
                sgrad.push_back(l);
                x1 += grad[0];
                x2 += grad[1];
            }
        }
    }

    double res = 0.0;
    for (size_t i=0; i<sgrad.size(); i++){
        res += sgrad[i];
    }
    /*std::unique_ptr<double []> v = std::make_unique<double[]>(model.dim);
    v[0] = x1; v[1] = x2;
    res = len(v.get(), model.dim);*/
    return res;
}
