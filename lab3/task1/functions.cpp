#include "main.h"

double function(const double x) {
    return std::log(x);
}

std::vector<double> Lagrange_polynomial(const func y, const std::vector<double>& x_points, const double x_precision) {
    const int n = x_points.size();

    std::vector<double> coefficients(n);

    for (int i = 0; i < n; ++i) {
        double w = 1;

        for (int j = 0; j < n; ++j) {
            if (i == j)
                continue;

            w *= (x_points[i] - x_points[j]);
        }

        coefficients[i] = y(x_points[i]) / w;
    }

    double approximation = 0;

    std::string out{};

    for (int i = 0; i < n; ++i) {
        double x = 1;
        std::string x_out{};

        for (int j = 0; j < n; ++j) {
            if (i == j)
                continue;

            x *= (x_precision - x_points[j]);

            x_out += std::format("(x - {:.1f})", x_points[j]);

        }

        approximation += coefficients[i] * x;

        if (coefficients[i] == 0 || x == 0) continue;

        if (coefficients[i] < 0)
            out += " - ";
        else if (i > 0)
            out += " + ";

        out += std::format("{:.5f}{}", std::abs(coefficients[i]), x_out);
    }

    std::cout << out << std::endl << "Result in X*: " << approximation << std::endl;

    return coefficients;
}

double separated_difference(double f_i, double f_j, double x_i, double x_j) {
    return (f_i - f_j) / (x_i - x_j);
}

std::vector<double> Newton_polynomial(const func y, const std::vector<double>& x_points, const double x_precision) {
    const int n = x_points.size();

    std::vector<double> f_cur(x_points.size());

    for (int i = 0; i < n; ++i) {
        f_cur[i] = y(x_points[i]);
    }


    for (int i = 0; i < n; ++i) {
        for (int j = n - 1; j > i; --j) {
            f_cur[j] = separated_difference(f_cur[j - 1], f_cur[j], x_points[j - 1 - i], x_points[j]);
        }
    }

    double approximation = 0;
    std::string out{};

    for (int i = 0; i < n; ++i) {
        double x = 1;
        std::string x_out{};

        for (int j = 0; j < i; ++j) {
            x *= x_precision - x_points[j];

            x_out += std::format("(x - {:.1f})", x_points[j]);
        }

        approximation += f_cur[i] * x;

        if (f_cur[i] == 0 || x == 0) continue;

        if (f_cur[i] < 0)
            out += " - ";
        else if (i > 0)
            out += " + ";

        out += std::format("{:.5f}{}", std::abs(f_cur[i]), x_out);
    }

    std::cout << out << std::endl << "Result in X*: " << approximation << std::endl;

    return f_cur;
}