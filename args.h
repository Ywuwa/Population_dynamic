#pragma once

#include <pthread.h>

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

// {} is brace-initializer https://en.cppreference.com/w/cpp/language/aggregate_initialization,
// It's used to initialize object and prevent using non initialized variables.
struct args {
    int n                 {};
    int m                 {};
    int r                 {};
    int k                 {};
    int p                 {};
    int s                 {};
    int* res              {};
    double* normA         {};
    double* mas_norm      {};
    int* mas_q            {};
    int* mas_err          {};
    int* c                {};
    double* A             {};
    double* E             {};
    char* filename        {};
    int* mas_num          {};
    pthread_barrier_t* b  {};
    pthread_mutex_t* mut  {};
    double elapsed_thread {};
    double elapsed_total  {};
};

struct Model
{
    int   dim    {};   // dimensions of the problem
    int   rtype  {};   // resource type: 1 logistic, 0 chemostat
    int   k0     {};
    float sig_a  {};
    float sig_g  {};
    float eps    {};
    float xi     {};
    float birth  {};
    float death  {};
    float dissip {};

    std::unique_ptr<float  []> m;
    std::unique_ptr<double []> R; // resource dense
    std::unique_ptr<double []> y;

    ~Model() {
        printf("deleting the model...\n");
    };
};

struct population
{
    int dim  {};
    float lambda = 1.0;
    double N {}; // population density
    double p {}; // predation coef
    double r {}; // consumption coef

    std::unique_ptr<double []> x; // phenotype vector

    population () = default;
    population (const population &other);
    population (population &&other) noexcept;

    void init (int idim, double iN, double ip, double *ix);
    void mutate (const population &P);
    void descrip () const;
};
