#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <cmath>
#include <iostream>
#include <vector>
#include <format>

using func = double(*)(double);

double function(double x);

std::vector<double> Lagrange_polynomial(func y, const std::vector<double>& x_points, double x_precision);
std::vector<double> Newton_polynomial(const func y, const std::vector<double>& x_points, const double x_precision);

#endif //NUMERIC_METHODS_MAIN_H