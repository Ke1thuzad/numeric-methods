#include "main.h"

double separated_difference(double f_i, double f_j, double x_i, double x_j) {
    return (f_i - f_j) / (x_i - x_j);
}

double compute_first_derivative(const std::vector<double>& x, const std::vector<double>& y, double x_point) {
    int x_target_i = -1;

    for (int i = 0; i < x.size(); ++i) {
        if (x[i] == x_point) {
            x_target_i = i;
            break;
        }
    }

    if (x_target_i == -1)
        throw std::exception();

    if (x_target_i == 0) {
        return separated_difference(y[x_target_i + 1], y[x_target_i], x[x_target_i + 1], x[x_target_i]);
    }

    if (x_target_i == (x.size() - 1)) {
        return separated_difference(y[x_target_i], y[x_target_i - 1], x[x_target_i], x[x_target_i - 1]);
    }

    double left_sided_derivative = separated_difference(y[x_target_i], y[x_target_i - 1], x[x_target_i], x[x_target_i - 1]);
    double right_sided_derivative = separated_difference(y[x_target_i + 1], y[x_target_i], x[x_target_i + 1], x[x_target_i]);

    double frac = separated_difference(right_sided_derivative, left_sided_derivative, x[x_target_i + 1], x[x_target_i - 1]);

    return left_sided_derivative + frac * (2 * x_point - x[x_target_i - 1] - x[x_target_i]);
}

double compute_second_derivative(const std::vector<double>& x, const std::vector<double>& y, double x_point) {
    int x_target_i = -1;

    for (int i = 0; i < x.size(); ++i) {
        if (x[i] == x_point) {
            x_target_i = i;
            break;
        }
    }

    double left_sided_derivative = separated_difference(y[x_target_i], y[x_target_i - 1], x[x_target_i], x[x_target_i - 1]);
    double right_sided_derivative = separated_difference(y[x_target_i + 1], y[x_target_i], x[x_target_i + 1], x[x_target_i]);

    double frac = separated_difference(right_sided_derivative, left_sided_derivative, x[x_target_i + 1], x[x_target_i - 1]);

    return 2 * frac;
}
