#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <vector>
#include <functional>
#include <string>
#include <format>


struct Solution {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

double f(double x, double y, double z);
double g(double x, double y, double z);

double exact_y(double x);
double exact_z(double x);

Solution euler_explicit(double x0, double y0, double z0, double h, int steps);
Solution euler_cauchy(double x0, double y0, double z0, double h, int steps);
Solution euler_improved(double x0, double y0, double z0, double h, int steps);
Solution runge_kutta_4(double x0, double y0, double z0, double h, int steps);
Solution adams_4(double x0, double y0, double z0, double h, int steps);

void print_solution(const Solution& sol, const std::string& method_name);
void save_to_file(const Solution& sol, const std::string& filename);
double calculate_error(const Solution& num, const Solution& exact);
double runge_romberg_error(const Solution& sol_h, const Solution& sol_2h, int p);

#endif //NUMERIC_METHODS_MAIN_H