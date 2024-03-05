#pragma once

#include "args.h"

int model_input(Model &model, const std::string &filename);

int dense_input(Model &model, const std::string &filename);
int dense_output(const Model &model, const std::string &filename);

int populationz_input(std::vector<population> &poplist, int dim, const std::string &filename);
int populationz_output(const std::vector<population> &poplist, int dim, const std::string &filename);

void dense_fill(const Model &model, double value);
void info(std::vector<population> &poplist);

int grad_output(const std::vector<double> &grad, const std::string &filename);
