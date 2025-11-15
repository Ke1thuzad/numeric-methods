#include "main.h"

double rectangle_method(func f, double x_beg, double x_end, double h) {
    double prev_x = x_beg;
    double result = 0;

    for (double x = x_beg + h; x <= x_end; x += h) {
        result += f((prev_x + x) / 2);

        prev_x = x;
    }

    return result * h;
}

double trapezoid_method(func f, double x_beg, double x_end, double h) {
    double result = (f(x_beg) + f(x_end)) / 2;

    for (double x = x_beg + h; x < x_end; x += h) {
        result += f(x);
    }

    return h * result;
}

double simpson_method(func f, double x_beg, double x_end, double h) {
    double result = f(x_beg) + f(x_end);

    int i = 1;

    for (double x = x_beg + h; x < x_end; x += h) {
        double multiplier = 2;
        if (i % 2)
            multiplier = 4;

        result += multiplier * f(x);

        i++;
    }

    return (h / 3) * result;
}

double runge_romberg_richardson_error_method(double F, double Fh1, double Fh2, double k, double p) {
    return std::abs(F - Fh2 - (Fh2 - Fh1) / (std::pow(k, p) - 1));
}