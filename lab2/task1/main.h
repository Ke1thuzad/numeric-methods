#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <iostream>
#include <cmath>

using func = double(*)(double);

double f1(double x);
double f1_dx(double x);
double phi(double x);
double Newton_method(double x0, func f, func f_dx, double eps);
double simple_iterations_method(double x0, func phi, double q, double eps);

#endif //NUMERIC_METHODS_MAIN_H