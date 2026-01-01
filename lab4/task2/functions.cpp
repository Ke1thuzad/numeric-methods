#include "main.h"
#include "matrix.h"
#include "../../lab1/task2/main.h"
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>


double f_shooting(double x, double y, double z) {
    return z;
}

double g_shooting(double x, double y, double z) {
    return 2.0 * y / (x*x + 1.0);
}

double exact_y(double x) {
    return x * x + 1;
}

Solution runge_kutta_4_system(double x0, double y0, double z0, double h, int steps) {
    Solution sol;
    sol.x.resize(steps + 1);
    sol.y.resize(steps + 1);
    sol.z.resize(steps + 1);
    
    sol.x[0] = x0;
    sol.y[0] = y0;
    sol.z[0] = z0;
    
    for (int i = 0; i < steps; i++) {
        double x = sol.x[i];
        double y = sol.y[i];
        double z = sol.z[i];
        
        double k1_y = h * f_shooting(x, y, z);
        double k1_z = h * g_shooting(x, y, z);
        
        double k2_y = h * f_shooting(x + h/2, y + k1_y/2, z + k1_z/2);
        double k2_z = h * g_shooting(x + h/2, y + k1_y/2, z + k1_z/2);
        
        double k3_y = h * f_shooting(x + h/2, y + k2_y/2, z + k2_z/2);
        double k3_z = h * g_shooting(x + h/2, y + k2_y/2, z + k2_z/2);
        
        double k4_y = h * f_shooting(x + h, y + k3_y, z + k3_z);
        double k4_z = h * g_shooting(x + h, y + k3_y, z + k3_z);
        
        sol.x[i+1] = x + h;
        sol.y[i+1] = y + (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6;
        sol.z[i+1] = z + (k1_z + 2*k2_z + 2*k3_z + k4_z) / 6;
    }
    
    return sol;
}

double shooting_residual(double eta, double a, double b, double h, double bc_right) {
    int steps = static_cast<int>((b - a) / h);

    Solution sol = runge_kutta_4_system(a, eta, 0.0, h, steps);

    double y_b = sol.y.back();
    double z_b = sol.z.back();

    return (y_b - z_b) - bc_right;
}

double shooting_method(double a, double b, double h, double bc_left, double bc_right, double eps) {
    double eta0 = 0.5;
    double eta1 = 1.5;

    double phi0 = shooting_residual(eta0, a, b, h, bc_right);
    double phi1 = shooting_residual(eta1, a, b, h, bc_right);

    int max_iter = 100;
    double eta_current = eta1;
    double phi_current = phi1;

    for (int iter = 0; iter < max_iter; iter++) {
        if (fabs(phi_current) < eps) {
            std::cout << "Метод стрельбы сошелся за " << iter << " итераций\n";
            return eta_current;
        }

        double eta_next = eta_current - phi_current * (eta_current - eta0) / (phi_current - phi0);

        eta0 = eta_current;
        phi0 = phi_current;
        eta_current = eta_next;
        phi_current = shooting_residual(eta_current, a, b, h, bc_right);
    }

    return eta_current;
}

double f_diff_eq(double x, double y) {
    return 0.0;
}

double p_diff_eq(double x) {
    return 0.0;
}

double q_diff_eq(double x) {
    double tan_x = tan(x);
    return -2.0 * (1.0 + tan_x * tan_x);
}

BoundarySolution finite_difference_method(double a, double b, int N) {
    BoundarySolution solution;
    double h = (b - a) / N;

    solution.x.resize(N + 1);
    solution.y.resize(N + 1);

    for (int i = 0; i <= N; i++) {
        solution.x[i] = a + i * h;
    }

    // Для этого уравнения: (x^2+1)y'' - 2y = 0
    // Разностная схема: (x_i^2 + 1) * (y_{i-1} - 2y_i + y_{i+1})/h^2 - 2y_i = 0
    // Умножаем на h^2: (x_i^2 + 1)*(y_{i-1} - 2y_i + y_{i+1}) - 2h^2*y_i = 0

    Matrix<double> mat(N + 1, N + 1);
    double* d = new double[N + 1];

    for (int i = 0; i <= N; i++) {
        d[i] = 0.0;
    }

    for (int i = 1; i < N; i++) {
        double x = solution.x[i];
        double A = x * x + 1.0;

        // Коэффициенты для разностного уравнения:
        // A*(y_{i-1} - 2y_i + y_{i+1}) - 2h^2*y_i = 0
        // => A*y_{i-1} - (2A + 2h^2)*y_i + A*y_{i+1} = 0

        mat.coefficients[i][i-1] = A;
        mat.coefficients[i][i] = -2.0 * A - 2.0 * h * h;
        mat.coefficients[i][i+1] = A;
        d[i] = 0.0;
    }

    // Граничное условие в точке x = 0: y'(0) = bc_left = 0
    // Аппроксимация: (y_1 - y_{-1})/(2h) = 0 => y_1 = y_{-1}
    // Разностное уравнение в точке i=0:
    // A0*(y_{-1} - 2y_0 + y_1) - 2h^2*y_0 = 0
    // Подставляем y_{-1} = y_1:
    // A0*(y_1 - 2y_0 + y_1) - 2h^2*y_0 = 0
    // 2A0*y_1 - (2A0 + 2h^2)*y_0 = 0

    double x0 = solution.x[0];
    double A0 = x0 * x0 + 1.0;
    mat.coefficients[0][0] = -(2.0 * A0 + 2.0 * h * h);
    mat.coefficients[0][1] = 2.0 * A0;
    d[0] = 0.0;

    // Граничное условие в точке x = 2: y(2) - y'(2) = bc_right = 1
    // Аппроксимация: y_N - (y_{N+1} - y_{N-1})/(2h) = 1
    // Выражаем y_{N+1} из граничного условия:
    // y_N - (y_{N+1} - y_{N-1})/(2h) = 1
    // y_{N+1} = y_{N-1} + 2h*(y_N - 1)

    // Разностное уравнение в точке i=N:
    // AN*(y_{N-1} - 2y_N + y_{N+1}) - 2h^2*y_N = 0
    // Подставляем выражение для y_{N+1}:
    // AN*(y_{N-1} - 2y_N + y_{N-1} + 2h*(y_N - 1)) - 2h^2*y_N = 0
    // AN*(2y_{N-1} - 2y_N + 2h*y_N - 2h) - 2h^2*y_N = 0
    // 2AN*y_{N-1} + (-2AN + 2h*AN - 2h^2)*y_N - 2h*AN = 0
    // 2AN*y_{N-1} + (-2AN + 2h*AN - 2h^2)*y_N = 2h*AN

    double xN = solution.x[N];
    double AN = xN * xN + 1.0;
    mat.coefficients[N][N-1] = 2.0 * AN;
    mat.coefficients[N][N] = -2.0 * AN + 2.0 * h * AN - 2.0 * h * h;
    d[N] = 2.0 * h * AN;

    double* y_result = thomas_algorithm(mat, d);

    for (int i = 0; i <= N; i++) {
        solution.y[i] = y_result[i];
    }

    delete[] y_result;
    delete[] d;

    return solution;
}

void save_solution_to_file(const BoundarySolution& sol, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "x,y,y_exact\n";
    for (size_t i = 0; i < sol.x.size(); i++) {
        file << sol.x[i] << "," << sol.y[i] << "," << exact_y(sol.x[i]) << "\n";
    }
    
    file.close();
}

double runge_romberg_error(const Solution& sol_h, const Solution& sol_2h, int p) {
    double max_error_est = 0;
    for (size_t i = 0; i < sol_2h.x.size(); i++) {
        for (size_t j = 0; j < sol_h.x.size(); j++) {
            if (fabs(sol_h.x[j] - sol_2h.x[i]) < 1e-10) {
                double error_est = fabs(sol_h.y[j] - sol_2h.y[i]) / (pow(2, p) - 1);
                if (error_est > max_error_est) {
                    max_error_est = error_est;
                }
                break;
            }
        }
    }
    return max_error_est;
}

double runge_romberg_error(const BoundarySolution& sol_h, const BoundarySolution& sol_2h, int p) {
    double max_error_est = 0;
    for (size_t i = 0; i < sol_2h.x.size(); i++) {
        for (size_t j = 0; j < sol_h.x.size(); j++) {
            if (fabs(sol_h.x[j] - sol_2h.x[i]) < 1e-10) {
                double error_est = fabs(sol_h.y[j] - sol_2h.y[i]) / (pow(2, p) - 1);
                if (error_est > max_error_est) {
                    max_error_est = error_est;
                }
                break;
            }
        }
    }
    return max_error_est;
}

double compute_absolute_error(const std::vector<double>& x_vals, const std::vector<double>& y_vals) {
    double max_error = 0.0;
    for (size_t i = 0; i < x_vals.size(); i++) {
        double exact = exact_y(x_vals[i]);
        double error = fabs(y_vals[i] - exact);
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}
