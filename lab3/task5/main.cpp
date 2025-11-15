#include "main.h"
#include <iostream>
#include <cmath>

double metod_func(double x) {
    return x / std::pow((3 * x + 4), 2);
}

double var_func(double x) {
    return 1 / ((2 * x + 7) * (3 * x + 4));
}

int main() {
    double x_beg = -1;
    double x_end = 1;
    double h1 = 0.5;
    double h2 = 0.25;
    double F_exact = 0.104471;

    std::cout << "=== h1 = " << h1 << " ===" << std::endl;
    double Fh1_rect = rectangle_method(var_func, x_beg, x_end, h1);
    double Fh1_trap = trapezoid_method(var_func, x_beg, x_end, h1);
    double Fh1_simp = simpson_method(var_func, x_beg, x_end, h1);

    std::cout << "Rectangle: " << Fh1_rect << std::endl;
    std::cout << "Trapezoid: " << Fh1_trap << std::endl;
    std::cout << "Simpson: " << Fh1_simp << std::endl << std::endl;

    std::cout << "=== h2 = " << h2 << " ===" << std::endl;
    double Fh2_rect = rectangle_method(var_func, x_beg, x_end, h2);
    double Fh2_trap = trapezoid_method(var_func, x_beg, x_end, h2);
    double Fh2_simp = simpson_method(var_func, x_beg, x_end, h2);

    std::cout << "Rectangle: " << Fh2_rect << std::endl;
    std::cout << "Trapezoid: " << Fh2_trap << std::endl;
    std::cout << "Simpson: " << Fh2_simp << std::endl << std::endl;

    std::cout << "=== Runge-Romberg Error ===" << std::endl;
    double error_rect = runge_romberg_richardson_error_method(F_exact, Fh1_rect, Fh2_rect, h1 / h2, 2);
    double error_trap = runge_romberg_richardson_error_method(F_exact, Fh1_trap, Fh2_trap, h1 / h2, 2);
    double error_simp = runge_romberg_richardson_error_method(F_exact, Fh1_simp, Fh2_simp, h1 / h2, 4);

    std::cout << "Rectangle: " << error_rect << std::endl;
    std::cout << "Trapezoid: " << error_trap << std::endl;
    std::cout << "Simpson: " << error_simp << std::endl;

    return 0;
}