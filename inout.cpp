#include "inout.h"

#include <string.h>

//===================================== MODEL INPUT =================================================
int model_input(Model &model, const std::string &filename) {

    FILE* IN;
    char name[100];
    int dim, restype, k0;
    float sig_a, sig_g, eps, xi, birth, death, dissip;
    float m;

    if ((IN = fopen(filename.c_str(), "r")) == NULL) {
        printf("Can't open file\n");
        return -2;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("The first line is missing\n");
        fclose(IN);
        return -1;
    }
    // 2 DIM
    if (fscanf(IN, "%d", &dim) != 1) {
        printf("Wrong input dim\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("dim line is missing\n");
        fclose(IN);
        return -1;
    }
    // RESource TYPE
    if (fscanf(IN, "%d", &restype) != 1) {
        printf("Wrong input restype\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("restype line is missing\n");
        fclose(IN);
        return -1;
    }
    // 3 K0
    if (fscanf(IN, "%d", &k0) != 1) {
        printf("Wrong input dim\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("K0 line is missing\n");
        fclose(IN);
        return -1;
    }
    // 4 SIG_A
    if (fscanf(IN, "%f", &sig_a) != 1) {
        printf("Wrong input sig_a\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("sig_a line is missing\n");
        fclose(IN);
        return -1;
    }
    // 5 SIG_G
    if (fscanf(IN, "%f", &sig_g) != 1) {
        printf("Wrong input sig_g\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("sig_g line is missing\n");
        fclose(IN);
        return -1;
    }
    // 6 EPS
    if (fscanf(IN, "%f", &eps) != 1) {
        printf("Wrong input eps\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("eps line is missing\n");
        fclose(IN);
        return -1;
    }
    // 7 XI
    if (fscanf(IN, "%f", &xi) != 1) {
        printf("Wrong input xi\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("xi line is missing\n");
        fclose(IN);
        return -1;
    }
    // 8 BIRTH
    if (fscanf(IN, "%f", &birth) != 1) {
        printf("Wrong input birth\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("birth line is missing\n");
        fclose(IN);
        return -1;
    }
    // 9 DEATH
    if (fscanf(IN, "%f", &death) != 1) {
        printf("Wrong input death\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("death line is missing\n");
        fclose(IN);
        return -1;
    }
    // DISSIPATION
    if (fscanf(IN, "%f", &dissip) != 1) {
        printf("Wrong input dissipation\n");
        fclose(IN);
        return -1;
    }
    if (fgets(name, 100, IN) == 0) {
        printf("dissipation line is missing\n");
        fclose(IN);
        return -1;
    }
    // 10 M
    if (fgets(name, 100, IN) == 0) {
        printf("m line is missing\n");
        fclose(IN);
        return -1;
    }

    model.m = std::make_unique<float []> (dim);
    model.y = std::make_unique<double []>(dim);
    for (int i = 0; i < dim; i++) {
        if (fscanf(IN, "%f", &m) != 1) {
            printf("Wrong input m, %d\n", i);
            fclose(IN);
            return -1;
        }
        model.m[i] = m;
        model.y[i] = 0.0;
    }

    model.dim = dim;
    model.rtype = restype;
    model.k0 = k0;
    model.sig_a = sig_a;
    model.sig_g = sig_g;
    model.eps = eps;
    model.xi = xi;
    model.birth = birth;
    model.death = death;
    model.dissip = dissip;

    fclose(IN);
    return 0;
}

//===================================== RESOURSE DENSITY DISTRIBUTION INPUT =================================================
int dense_input(Model &model, const std::string &filename) {
    FILE* IN;
    if ((IN = fopen(filename.c_str(), "r")) == NULL) {
        printf("Can't open file2\n");
        return -2;
    }
    int n = 1;
    for (int i = 0; i < model.dim; i++) { n *= 60; } // 60 grid knot over each dimension

    model.R = std::make_unique<double []>(n);
    for (int i = 0; i < n; i++) {
        if (fscanf(IN, "%lf ", &model.R[i]) != 1) { return -1; }
    }
    fclose(IN);
    return 0;
}

//===================================== RESOURSE DENSITY DISTRIBUTION OUTPOT ================================================
int dense_output(const Model &model, const std::string &filename) {
    FILE* OUT;
    if ((OUT = fopen(filename.c_str(), "w")) == NULL) {
        printf("Can't open file\n");
        return -2;
    }
    int n = 1;
    for (int i = 0; i < model.dim; i++) { n *= 60; } // 60 grid knot over each dimension
    for (int i = 0; i < n; i++) {
        fprintf(OUT, "%10.5e ", model.R[i]);
    }
    fclose(OUT);
    return 0;
}

//============================================== POPULATIONZ INPUT ==========================================================
int populationz_input(std::vector<population> &poplist, const int dim, const std::string &filename) {
    FILE* IN;

    if ((IN = fopen(filename.c_str(), "r")) == NULL) {
        printf("Can't open file\n");
        return -2;
    }

    double N;
    while (true) {
        if (fscanf(IN, "%lf", &N) != 1) {
            printf("density wrong input\n");
            fclose(IN);
            return -1;
        }

        std::unique_ptr<double[]> x = std::make_unique<double[]> (dim);
        for (int i = 0; i < dim; i++) {
            if (fscanf(IN, "%lf", &x[i]) != 1) {
                printf("phenotype component wrong input, %10.5e\n", x[i]);
                fclose(IN);
                return -1;
            }
        }

        double p;
        if (fscanf(IN, "%lf", &p) != 1) {
            printf("p-coef wrong input\n");
            fclose(IN);
            return -1;
        }

        population P;
        P.init(dim, N, p, x.get());
        poplist.push_back(std::move(P));

        if (fgetc(IN) == EOF) {
            break;
        }
    }
    fclose(IN);
    return 0;
}

//============================================== POPULATIONZ OUTPUT =========================================================
int populationz_output(const std::vector<population> &poplist, const int dim, const std::string &filename) {
    FILE* OUT;
    if ((OUT = fopen(filename.c_str(), "w")) == NULL) {
        printf("Can't open file\n");
        return -2;
    }

    for (size_t i = 0; i < poplist.size()-1; i++) {
        fprintf(OUT, "%10.5e ", poplist[i].N);
        for (int j = 0; j < dim; j++) {
            fprintf(OUT, "%10.5e ", poplist[i].x[j]);
        }
        fprintf(OUT, "%10.5e ", poplist[i].p);
        fprintf(OUT, "\n");
    }

    fprintf(OUT, "%10.5e ", poplist[poplist.size()-1].N);
    for (int j = 0; j < dim; j++) {
        fprintf(OUT, "%10.5e ", poplist[poplist.size()-1].x[j]);
    }
    fprintf(OUT, "%10.5e", poplist[poplist.size()-1].p);
    fclose(OUT);
    return 0;
}
void dense_fill(const Model &model, double value)
{
    int n = 1;
    for (int i = 0; i < model.dim; i++) { n *= 60; } // 60 grid knot over each dimension
    //model.R = std::make_unique<double []>(n);
    for (int i = 0; i < n; i++) {
        model.R[i] = value;
    }
}

//============================ SELECTIVE GRADIENT METRIC OUTPOT =================================
int grad_output(const std::vector<double> &grad, const std::string &filename) {
    FILE* OUT;
    if ((OUT = fopen(filename.c_str(), "w")) == NULL) {
        printf("Can't open file\n");
        return -2;
    }
    for (size_t i = 0; i < grad.size(); i++) {
        fprintf(OUT, "%10.5e ", grad[i]);
    }
    fclose(OUT);
    return 0;
}
//====================================== MINMAX INFO ===========================================
void info(std::vector<population> &poplist)
{
    double max_n = 0.0;
    double max_nx = 0.0;
    double max_ny = 0.0;
    double max_x = -4.0;
    double max_y = -4.0;
    double min_y = 4.0;
    double min_x = 4.0;

    double max_p = 0.0;
    double max_px = 0.0;
    double max_py = 0.0;

    for (size_t i=0; i<poplist.size(); i++){
        if (poplist[i].N > max_n) {max_n = poplist[i].N; max_nx=poplist[i].x[0]; max_ny=poplist[i].x[1];}
        if (poplist[i].x[0] > max_x) max_x = poplist[i].x[0];
        if (poplist[i].x[0] < min_x) min_x = poplist[i].x[0];
        if (poplist[i].x[1] > max_y) max_y = poplist[i].x[1];
        if (poplist[i].x[1] < min_y) min_y = poplist[i].x[1];

        if (poplist[i].p > max_p) {max_p = poplist[i].p; max_px=poplist[i].x[0]; max_py=poplist[i].x[1];}
    }
    printf("Situation: \n");
    printf("//------- max_N: %10.3e, (%10.3e, %10.3e)\n", max_n, max_nx, max_ny);
    printf("//------- max_p: %10.3e, (%10.3e, %10.3e)\n", max_p, max_px, max_py);
    printf("//------- max_x, max_y: %10.3e, %10.3e\n", max_x, max_y);
    printf("//------- min_x, min_y: %10.3e, %10.3e\n", min_x, min_y);
    printf("//-------------------------------------------------------\n");
}
