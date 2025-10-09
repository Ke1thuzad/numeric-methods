#ifndef NUMERAL_METHODS_MAIN_H
#define NUMERAL_METHODS_MAIN_H

#include <bits/stdc++.h>
#include "matrix.h"

std::pair<Matrix<double>, double*> read_tridiagonal_matrix_from_file(const char *path);

double* thomas_algorithm(const Matrix<double>& mat, const double* d);
void thomas_propagation(const Matrix<double>& mat, double* P, double* Q, const double* d);

#endif //NUMERAL_METHODS_MAIN_H