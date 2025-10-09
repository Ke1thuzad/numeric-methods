#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <matrix.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <limits>

template<class T> Matrix<T> read_matrix_from_file(const char *path);

template<class T>
std::pair<Matrix<T>, Matrix<T>> QR_Decomposition(Matrix<T>& A);

Matrix<double> QR_eigenvalues(const Matrix<double>& A, const double eps);

#endif //NUMERIC_METHODS_MAIN_H