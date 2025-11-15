#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <iostream>
#include <vector>
#include <algorithm>

double compute_first_derivative(const std::vector<double>& x, const std::vector<double>& y, double x_point);
double compute_second_derivative(const std::vector<double>& x, const std::vector<double>& y, double x_point);

#endif //NUMERIC_METHODS_MAIN_H