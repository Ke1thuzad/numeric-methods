#ifndef NUMERAL_METHODS_MAIN_H
#define NUMERAL_METHODS_MAIN_H

#include <iostream>
#include <fstream>

#include "matrix.h"

template <class T>
struct LUPResult {
    Matrix<T> L;
    Matrix<T> U;
    Matrix<T> P;
    int swaps;

    explicit LUPResult(int n) : L(n, n), U(n, n), P(n, n), swaps(0) {}
};

LinearEquation<double> read_equation_from_file(const char *path);

template <class T>
LUPResult<T> LUPDecomposition(const Matrix<T>&);
template <class T>
T* solveLUP(const LUPResult<T>&, const T*);
template <class T>
T determinantLUP(const LUPResult<T>&);
template <class T>
Matrix<T> inverseMatrix(const Matrix<T>&);

#endif //NUMERAL_METHODS_MAIN_H