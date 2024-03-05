#include "args.h"

#include <random>
#include <string>

population::population(const population &other) {
    dim = other.dim;
    N = other.N;
    p = other.p;
    r = other.r;

    x = std::make_unique<double[]>(dim);
    std::copy_n(other.x.get(), dim, x.get());
}

population::population(population &&other) noexcept {
    dim = other.dim;
    N = other.N;
    p = other.p;
    r = other.r;
    x = std::move(other.x);
}

void population::init(int idim, double iN, double ip, double* ix) {
    dim = idim;
    N = iN;
    p = ip;
    r = pow(1 - pow(p, lambda), 1.0 / lambda);
    x = std::make_unique<double []>(dim);

    for (int i = 0; i < dim; i++) {
        x[i] = ix[i];
    }
}

void population::mutate(const population&P) {
    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    const double mu = 0.01;
    std::normal_distribution<double> distribution(0.0, mu);

    dim = P.dim;
    N = 2e-06;
    p = P.p + distribution(generator);
    if (p < 0) { p = 0.0; }
    if (p > 1) { p = 1.0; }
    r = pow(1 - pow(p, lambda), 1.0 / lambda);
    x = std::make_unique<double []>(dim);
    for (int i = 0; i < dim; i++) {
        x[i] = P.x[i] + distribution(generator);
    }
}

void population::descrip() const {
    printf("pop density: %10.3e, ", N);
    printf("predation/consumption coef: %10.3e/%10.3e, ", p, r);
    for (int i = 0; i < dim; i++) {
        printf("x[%d]: %10.3e, ", i, x[i]);
    }
    printf("\n");
}
