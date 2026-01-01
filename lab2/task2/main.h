#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include "../../lab1/task1/lup.h"

#include <matrix.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <format>

using func = double(*)(const Matrix<double>&);
using mx_func = Matrix<double>(*)(const Matrix<double>&);

Matrix<double> Jacobi(const Matrix<double>& x);
Matrix<double> A1(const Matrix<double>& x);
Matrix<double> A2(const Matrix<double>& x);

Matrix<double> simple_iterations_method_system(const Matrix<double>& starting_vector, double q, const std::vector<func>& phi, double eps);
Matrix<double> Newton_method_system(const Matrix<double>& starting_vector, const std::vector<mx_func>& A, mx_func J, const double eps);

#endif //NUMERIC_METHODS_MAIN_H