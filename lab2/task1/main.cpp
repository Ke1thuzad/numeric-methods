#include "main.h"

#include <format>

int main() {
    constexpr double x0 = 1.5;
    constexpr double q = 0.5;
    constexpr double eps = std::numeric_limits<double>::epsilon();

    std::cout << std::format("{:-^27}", "Newton Method") << std::endl;

    double root = Newton_method(x0, f1, f1_dx, eps);

    std::cout << std::format("Root: {:.6f}", root) << std::endl;
    std::cout << std::format("{:-^27}", "") << std::endl;

    std::cout << std::format("{:-^27}", "Simple Iterations Method") << std::endl;

    root = simple_iterations_method(x0, phi, q, eps);

    std::cout << std::format("Root: {:.6f}", root) << std::endl;
    std::cout << std::format("{:-^27}", "") << std::endl;

    return 0;
}

double f1(const double x) {
    return std::cos(x) + 0.25 * x - 0.5;
}

double f1_dx(const double x) {
    return -std::sin(x) + 0.25;
}

double phi(const double x) {
    return std::acos(0.5 - 0.25 * x);
}

double Newton_method(const double x0, const func f, const func f_dx, const double eps) {
    double x_k = x0;
    double delta = 1;

    int i = 0;

    std::cout << std::format("{:^16} {:^10}", "Error", "Iteration") << std::endl;

    while (delta >= eps) {
        ++i;

        const double x_k_prev = x_k;

        x_k = x_k_prev - f(x_k_prev) / f_dx(x_k_prev);


        delta = std::abs(x_k - x_k_prev);

        std::cout << std::format("{:^16.10f} {:^10}", delta, i) << std::endl;
    }

    std::cout << std::endl;

    return x_k;
}

double simple_iterations_method(const double x0, const func phi, const double q, const double eps) {
    double x_k = phi(x0);
    double delta = 1;

    int i = 0;

    std::cout << std::format("{:^16} {:^10}", "Error", "Iteration") << std::endl;

    while (delta > eps) {
        ++i;

        const double x_k_prev = x_k;

        x_k = phi(x_k);

        delta = q / (1 - q) * std::abs(x_k - x_k_prev);

        std::cout << std::format("{:^16.10f} {:^10}", delta, i) << std::endl;

    }

    return x_k;
}