#include "main.h"


int main() {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    constexpr double q = 0.9;

    std::cout << "Current epsilon: " << std::scientific << eps << std::endl << std::endl;

    const std::vector<func> phi = {
        [](const Matrix<double>& x) {  // 1 + cos(x2)
            return 1 + std::cos(x.coefficients[1][0]);
        },
        [](const Matrix<double>& x) {  // 2 + lg(x1 + 1)
            return 2 + std::log10(x.coefficients[0][0] + 1);
        }
    };

    const Matrix starting_vector {{0.45}, {2.15}};  // x1_0 = 0.45, x2_0 = 2.15

    std::cout << std::format("{:-^27}", "Newton Method") << std::endl;

    Matrix<double> root = Newton_method_system(starting_vector, {A1, A2}, Jacobi, eps);

    std::cout << "\nRoots:\n" << root;
    std::cout << std::format("{:-^27}", "") << std::endl;

    std::cout << std::format("{:-^27}", "Simple Iterations Method") << std::endl;

    root = simple_iterations_method_system(starting_vector, q, phi, eps);

    std::cout << "\nRoots:\n" << root;
    std::cout << std::format("{:-^27}", "") << std::endl;

    return 0;
}

Matrix<double> Jacobi(const Matrix<double>& x) {
    return {
        {1, std::sin(x.coefficients[1][0])},
        {-1 / ((x.coefficients[0][0] + 1) * std::log(10)), 1}
    };
}

Matrix<double> A1(const Matrix<double>& x) {
    double** xc = x.coefficients;
    return {
        {xc[0][0] - cos(xc[1][0]) - 1, sin(xc[1][0])},
        {xc[1][0] - log10(xc[0][0] + 1) - 2, 1}
    };
}

Matrix<double> A2(const Matrix<double>& x) {

    double** xc = x.coefficients;
    return {
        {1, xc[0][0] - cos(xc[1][0]) - 1},
        {-1 / ((xc[0][0] + 1) * log(10)), xc[1][0] - log10(xc[0][0] + 1) - 2}
    };
}

double norm(const Matrix<double>& m) {
    double mat_max = 0;

    for (int i = 0; i < m.n; ++i) {
        for (int j = 0; j < m.m; ++j) {
            mat_max = std::max(mat_max, std::abs(m.coefficients[i][j]));
        }
    }

    return mat_max;
}

Matrix<double> simple_iterations_method_system(const Matrix<double>& starting_vector, double q, const std::vector<func>& phi, const double eps) {
    const int n = starting_vector.n;

    Matrix<double> x_k(n, 1);

    for (int i = 0; i < n; ++i) {
        x_k.coefficients[i][0] = phi[i](starting_vector);
    }

    double delta = 1;

    int i = 0;

    std::cout << std::format("{:^16} {:^10}", "Error", "Iteration") << std::endl;

    while (delta > eps) {
        ++i;

        const Matrix<double> x_k_prev = x_k;

        for (int j = 0; j < n; ++j) {
            x_k.coefficients[j][0] = phi[j](x_k);
        }

        delta = q / (1 - q) * norm(x_k - x_k_prev);

        std::cout << std::format("{:^16.10f} {:^10}", delta, i) << std::endl;

    }

    return x_k;
}

Matrix<double> Newton_method_system(const Matrix<double>& starting_vector, const std::vector<mx_func>& A, const mx_func J, const double eps) {
    Matrix<double> x_k = starting_vector;

    double delta = 1;

    int i = 0;

    std::cout << std::format("{:^16} {:^10}", "Error", "Iteration") << std::endl;

    while (delta >= eps) {
        ++i;

        const Matrix<double> x_k_prev = x_k;

        for (int j = 0; j < A.size(); ++j) {
            x_k.coefficients[j][0] = x_k_prev.coefficients[j][0] - determinant(A[j](x_k_prev)) / determinant(J(x_k_prev));
        }


        delta = norm(x_k - x_k_prev);

        std::cout << std::format("{:^16.10f} {:^10}", delta, i) << std::endl;
    }

    return x_k;
}