#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <iostream>
#include <cmath>
#include <vector>
#include <format>

double find_h(std::vector<double> x, int index);

double cubic_spline(const std::vector<double>& x, const std::vector<double>& f, double point_x);

#endif //NUMERIC_METHODS_MAIN_H