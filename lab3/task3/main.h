#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include "../../lab1/task1/lup.h"

#include <vector>
#include <iostream>
#include <format>
#include <algorithm>
#include <numeric>
#include <cmath>

std::vector<double> approximation_polynomial_first_degree(const std::vector<double>& x, const std::vector<double>& y);
std::vector<double> approximation_polynomial_second_degree(const std::vector<double>& x, const std::vector<double>& y);

#endif //NUMERIC_METHODS_MAIN_H