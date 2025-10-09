#include "main.h"

int main(const int argc, char **argv) {
    if (argc == 2) {
        const LinearEquation<double> equation = read_equation_from_file(argv[1]);

        constexpr double eps = 10e-32;

        std::cout << "-----------Simple Iterations Method-----------" << std::endl;

        Matrix<double> result = simple_iterations_method(equation, eps);

        std::cout << result << equation.matrix * result;

        std::cout << "------------Seidel Method------------" << std::endl;

        result = seidel_method(equation, eps);

        std::cout << result << equation.matrix * result;
    } else {
        throw std::invalid_argument("Incorrect arguments. Add path to your test file");
    }

    return 0;
}

LinearEquation<double> read_equation_from_file(const char *path) {
    std::ifstream stream(path);
    if (!stream.is_open()) {
        throw std::runtime_error("Cannot open file");
    }

    int n;
    stream >> n;

    if (n <= 0) {
        throw std::domain_error("Matrix size cannot be less than 1");
    }

    LinearEquation<double> equation(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (!(stream >> equation.matrix.coefficients[i][j])) {
                throw std::runtime_error("File format is not correct");
            }
        }
    }

    for (int i = 0; i < n; i++) {
        if (!(stream >> equation.b[i])) {
            throw std::runtime_error("File format is not correct");
        }
    }

    return equation;
}

template<class T>
T norm(Matrix<T> A) {
    T max_row_sum = 0;

    for (int i = 0; i < A.n; ++i) {
        T row_sum = 0;
        for (int j = 0; j < A.n; ++j) {
            row_sum += std::abs(A.coefficients[i][j]);
        }

        max_row_sum = std::max(max_row_sum, row_sum);
    }

    return max_row_sum;
}

template <class T>
Matrix<T> simple_iterations_loop(const Matrix<T>& alpha, const Matrix<T>& beta, const T eps) {
    Matrix<T> x_k, x_k_prev = beta;

    T eps_k;

    int i = 0;

    do {
        ++i;

        x_k = beta + alpha * x_k_prev;

        T norm_alpha = norm(alpha);

        T norm_x_k = norm(x_k - x_k_prev);

        x_k_prev = x_k;

        if (norm_alpha >= 1) {
            eps_k = norm_x_k;
        } else {
            eps_k = norm_alpha / (1 - norm_alpha) * norm_x_k;
        }
    } while (eps_k > eps);

    std::cout << i << std::endl;

    return x_k_prev;
}

template<class T>
Matrix<T> simple_iterations_method(const LinearEquation<T> &equation, const T eps) {
    const int n = equation.n;

    if (determinant(equation.matrix) == 0)
        throw std::invalid_argument("Matrix is singular");

    Matrix<T> beta(n, 1);
    Matrix<T> alpha(n, n);

    get_alpha_beta_components(equation, alpha, beta);

    return simple_iterations_loop<T>(alpha, beta, eps);
}

template<class T>
Matrix<T> seidel_method(const LinearEquation<T> &equation, const T eps) {
    const int n = equation.n;

    if (determinant(equation.matrix) == 0)
        throw std::invalid_argument("Matrix is singular");

    Matrix<T> beta(n, 1);
    Matrix<T> alpha(n, n);

    get_alpha_beta_components(equation, alpha, beta);

    Matrix<T> B(n, n);
    Matrix<T> C(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i <= j) {
                B.coefficients[i][j] = alpha.coefficients[i][j];
                C.coefficients[i][j] = T(0);
            } else {
                B.coefficients[i][j] = T(0);
                C.coefficients[i][j] = alpha.coefficients[i][j];
            }
        }
    }

    Matrix<T> inv_B = inverseMatrix(Matrix<T>::Identity(n) - B);

    return simple_iterations_loop<T>(inv_B * C, inv_B * beta, eps);
}

template<class T>
void get_alpha_beta_components(LinearEquation<T> equation, Matrix<T> &alpha, Matrix<T> &beta) {
    const int n = equation.n;

    for (int i = 0; i < n; ++i) {
        const T diag_elem = equation.matrix.coefficients[i][i];
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                alpha.coefficients[i][i] = 0;
                continue;
            }

            alpha.coefficients[i][j] = -equation.matrix.coefficients[i][j] / diag_elem;
        }
        beta.coefficients[i][0] = equation.b[i] / diag_elem;
    }
}
