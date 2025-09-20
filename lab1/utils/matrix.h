#ifndef NUMERAL_METHODS_MATRIX_H
#define NUMERAL_METHODS_MATRIX_H

#include <initializer_list>

template<class T>
class Matrix {
public:
    int n, m;
    T **coefficients;

    Matrix(const int n = 2, const int m = 2) : n(n), m(m) {
        coefficients = new T*[n];
        for (int i = 0; i < n; i++) {
            coefficients[i] = new T[m]();
        }
    }

    Matrix(std::initializer_list<std::initializer_list<T>> list) : n(list.size()), m(list.begin()->size()) {
        coefficients = new T*[n];
        size_t i = 0;
        for (const auto& row : list) {
            if (row.size() != m) {
                throw std::invalid_argument("All rows must have the same length");
            }
            coefficients[i] = new T[m];
            std::copy(row.begin(), row.end(), coefficients[i]);
            i++;
        }
    }

    Matrix(const Matrix& other) : n(other.n), m(other.m) {
        coefficients = new T*[n];
        for (int i = 0; i < n; i++) {
            coefficients[i] = new T[m];
            for (int j = 0; j < m; j++) {
                coefficients[i][j] = other.coefficients[i][j];
            }
        }
    }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            for (int i = 0; i < n; i++) {
                delete[] coefficients[i];
            }
            delete[] coefficients;

            n = other.n;
            m = other.m;
            coefficients = new T*[n];
            for (int i = 0; i < n; i++) {
                coefficients[i] = new T[m];
                for (int j = 0; j < m; j++) {
                    coefficients[i][j] = other.coefficients[i][j];
                }
            }
        }
        return *this;
    }

    ~Matrix() {
        for (int i = 0; i < n; ++i) {
            delete[] coefficients[i];
        }

        delete[] coefficients;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Matrix& mat) {
        stream << "Matrix (N: " << mat.n << ", M: " << mat.m << ")" << std::endl;

        for (int i = 0; i < mat.n; i++) {
            for (int j = 0; j < mat.m; j++) {
                stream << mat.coefficients[i][j] << ' ';
            }
            stream << std::endl;
        }
        return stream;
    }
};

#endif //NUMERAL_METHODS_MATRIX_H