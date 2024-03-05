#pragma once

#include "args.h"

double dist(const double* x, const double* y, const double xp, const double yp, int dim);

void weights(double N, const std::vector<population> &poplist, population &new_p);

void pop_clust(std::vector<population> &poplist, double eps);

int pop_choice(std::vector<population> &poplist);

void pop_clust2(std::vector<population> &poplist, double eps);
