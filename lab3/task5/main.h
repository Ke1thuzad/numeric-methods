#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <cmath>
#include <iostream>

using func = double(*)(double);

double rectangle_method(func f, double x_beg, double x_end, double h);
double trapezoid_method(func f, double x_beg, double x_end, double h);
double simpson_method(func f, double x_beg, double x_end, double h);
double runge_romberg_richardson_error_method(double F, double Fh1, double Fh2, double k, double p);

#endif //NUMERIC_METHODS_MAIN_H