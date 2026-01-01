#ifndef NUMERIC_METHODS_MAIN_H
#define NUMERIC_METHODS_MAIN_H

#include <vector>
#include <functional>
#include <string>

struct Solution {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

struct BoundarySolution {
    std::vector<double> x;
    std::vector<double> y;
};

double f_shooting(double x, double y, double z);
double g_shooting(double x, double y, double z);
double exact_y(double x);

Solution runge_kutta_4_system(double x0, double y0, double z0, double h, int steps);
double shooting_residual(double eta, double a, double b, double h, double target_value);
double shooting_method(double a, double b, double h, double initial_slope, double target_value, double eps);

double f_diff_eq(double x, double y);
double p_diff_eq(double x);
double q_diff_eq(double x);
BoundarySolution finite_difference_method(double a, double b, int N);

double runge_romberg_error(const Solution& sol_h, const Solution& sol_2h, int p);
double runge_romberg_error(const BoundarySolution& sol_h, const BoundarySolution& sol_2h, int p);
double compute_absolute_error(const std::vector<double>& x_vals, const std::vector<double>& y_vals);

void save_solution_to_file(const BoundarySolution& sol, const std::string& filename);

#endif //NUMERIC_METHODS_MAIN_H