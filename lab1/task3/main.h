#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include "../task1/lup.h"

#include <matrix.h>
#include <iostream>
#include <fstream>
#include <vector>

LinearEquation<double> read_equation_from_file(const char *path);

template<class T>
T norm(Matrix<T> A);

template<class T>
Matrix<T> simple_iterations_method(const LinearEquation<T> &equation, T eps);

template<class T>
Matrix<T> seidel_method(const LinearEquation<T> &equation, T eps);

template<class T>
void get_alpha_beta_components(LinearEquation<T> equation, Matrix<T> &alpha, Matrix<T> &beta);

#endif //NUMERIC_METHODS_MAIN_H