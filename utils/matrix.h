#ifndef NUMERAL_METHODS_MATRIX_H
#define NUMERAL_METHODS_MATRIX_H

#include <cfloat>
#include <initializer_list>
#include <iomanip>

template<class T>
class Matrix {
public:
    int n, m;
    T **coefficients;

    explicit Matrix(const int n) {
        *this = Matrix(n, n);
    }

    Matrix(const int n = 2, const int m = 2) : n(n), m(m) {
        coefficients = new T*[n];
        for (int i = 0; i < n; i++) {
            coefficients[i] = new T[m](0);
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

    Matrix operator*(const Matrix& other) const {
        if (m != other.n) {
            throw std::invalid_argument("Matrix dimensions do not match for multiplication");
        }

        Matrix result(n, other.m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < other.m; j++) {
                result.coefficients[i][j] = 0;
                for (int k = 0; k < m; k++) {
                    result.coefficients[i][j] += coefficients[i][k] * other.coefficients[k][j];
                    if (std::abs(result.coefficients[i][j]) < FLT_EPSILON) {
                        result.coefficients[i][j] = 0;
                    }
                }
            }
        }
        return result;
    }

    Matrix operator*(const T other) const {
        Matrix result(n, m);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result.coefficients[i][j] = coefficients[i][j] * other;
            }
        }

        return result;
    }

    Matrix operator+(const Matrix& other) const {
        if (n != other.n && m != other.m)
            throw std::invalid_argument("Matrix dimensions do not match for summation");

        Matrix result(n, m);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result.coefficients[i][j] = coefficients[i][j] + other.coefficients[i][j];
            }
        }

        return result;
    }

    Matrix operator+(const T other) const {
        Matrix result(n, m);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result.coefficients[i][j] = coefficients[i][j] + other;
            }
        }

        return result;
    }

    Matrix operator-() const {
        Matrix result(n, m);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result.coefficients[i][j] = -coefficients[i][j];
            }
        }

        return result;
    }

    Matrix operator-(const Matrix& other) const {
        return *this + (-other);
    }

    Matrix operator-(const T other) const {
        return *this + (-other);
    }

    Matrix& operator*=(const Matrix& other) {
        *this = *this * other;
        return *this;
    }

    Matrix& operator+=(const Matrix& other) {
        *this = *this + other;
        return *this;
    }

    Matrix& operator-=(const Matrix& other) {
        *this = *this - other;
        return *this;
    }

    static Matrix Identity(int size) {
        if (size <= 0) {
            throw std::invalid_argument("Identity matrix size must be positive");
        }

        Matrix result(size, size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                result.coefficients[i][j] = (i == j) ? T(1) : T(0);
            }
        }
        return result;
    }

    Matrix transpose() const {
        Matrix result(m, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                result.coefficients[j][i] = coefficients[i][j];
            }
        }

        return result;
    }

    static Matrix Householder(Matrix v) {
        Matrix vT = v.transpose();

        Matrix squared = vT * v;

        squared.coefficients[0][0] = 2 / squared.coefficients[0][0];

        Matrix result = Matrix::Identity(v.n) - (v * vT) * squared.coefficients[0][0];

        return result;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Matrix& mat) {
        stream << "Matrix (N: " << mat.n << ", M: " << mat.m << ")" << std::endl;

        std::ios_base::fmtflags original_flags = stream.flags();
        std::streamsize original_precision = stream.precision();

        stream << std::fixed << std::setprecision(6);

        for (int i = 0; i < mat.n; i++) {
            for (int j = 0; j < mat.m; j++) {
                stream << std::setw(12) << mat.coefficients[i][j] << " ";
            }
            stream << std::endl;
        }

        stream.flags(original_flags);
        stream.precision(original_precision);

        return stream;
    }
};

template <class T>
class LinearEquation {
public:
    int n;
    Matrix<T> matrix;
    T* b;

    LinearEquation(int n = 4) : n(n), matrix(n, n), b(new T[n]()) {}

    LinearEquation(const LinearEquation& other) : n(other.n), matrix(other.matrix), b(new T[other.n]) {
        std::copy(other.b, other.b + n, b);
    }

    LinearEquation& operator=(const LinearEquation& other) {
        if (this != &other) {
            delete[] b;
            n = other.n;
            matrix = other.matrix;
            b = new T[n];
            std::copy(other.b, other.b + n, b);
        }
        return *this;
    }

    ~LinearEquation() {
        delete[] b;
    }
};

#endif //NUMERAL_METHODS_MATRIX_H