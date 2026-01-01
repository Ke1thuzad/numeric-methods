#include "main.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x, double y, double z) {
    return z;
}

double g(double x, double y, double z) {
    double tanx = tan(x);
    return (1.0 + 2.0 * tanx * tanx) * y;
}

double exact_y(double x) {
    return 1.0/cos(x) + sin(x) + x/cos(x);
}

double exact_z(double x) {
    const double cosx = cos(x);

    return sin(x)/(cosx*cosx) + cos(x) + (cosx + x*sin(x)/(cosx*cosx));
}

Solution euler_explicit(double x0, double y0, double z0, double h, int steps) {
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

        sol.x[i+1] = x + h;
        sol.y[i+1] = y + h * f(x, y, z);
        sol.z[i+1] = z + h * g(x, y, z);
    }

    return sol;
}

Solution euler_cauchy(double x0, double y0, double z0, double h, int steps) {
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

        double y_pred = y + h * f(x, y, z);
        double z_pred = z + h * g(x, y, z);
        double x_next = x + h;

        sol.x[i+1] = x_next;
        sol.y[i+1] = y + h/2.0 * (f(x, y, z) + f(x_next, y_pred, z_pred));
        sol.z[i+1] = z + h/2.0 * (g(x, y, z) + g(x_next, y_pred, z_pred));
    }

    return sol;
}

Solution euler_improved(double x0, double y0, double z0, double h, int steps) {
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

        double h2 = h / 2.0;
        double y_half = y + h2 * f(x, y, z);
        double z_half = z + h2 * g(x, y, z);
        double x_half = x + h2;

        sol.x[i+1] = x + h;
        sol.y[i+1] = y + h * f(x_half, y_half, z_half);
        sol.z[i+1] = z + h * g(x_half, y_half, z_half);
    }

    return sol;
}

Solution runge_kutta_4(double x0, double y0, double z0, double h, int steps) {
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

        double k1_y = h * f(x, y, z);
        double k1_z = h * g(x, y, z);

        double k2_y = h * f(x + h/2, y + k1_y/2, z + k1_z/2);
        double k2_z = h * g(x + h/2, y + k1_y/2, z + k1_z/2);

        double k3_y = h * f(x + h/2, y + k2_y/2, z + k2_z/2);
        double k3_z = h * g(x + h/2, y + k2_y/2, z + k2_z/2);

        double k4_y = h * f(x + h, y + k3_y, z + k3_z);
        double k4_z = h * g(x + h, y + k3_y, z + k3_z);

        sol.x[i+1] = x + h;
        sol.y[i+1] = y + (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6.0;
        sol.z[i+1] = z + (k1_z + 2*k2_z + 2*k3_z + k4_z) / 6.0;
    }

    return sol;
}

Solution adams_4(double x0, double y0, double z0, double h, int steps) {
    Solution start = runge_kutta_4(x0, y0, z0, h, 3);

    Solution sol;
    sol.x = start.x;
    sol.y = start.y;
    sol.z = start.z;

    vector<double> f_vals, g_vals;
    for (int i = 0; i <= 3; i++) {
        f_vals.push_back(f(sol.x[i], sol.y[i], sol.z[i]));
        g_vals.push_back(g(sol.x[i], sol.y[i], sol.z[i]));
    }

    for (int i = 3; i < steps; i++) {
        double x = sol.x[i];
        double y = sol.y[i];
        double z = sol.z[i];

        double y_next = y + h/24.0 * (55*f_vals[i] - 59*f_vals[i-1] + 37*f_vals[i-2] - 9*f_vals[i-3]);
        double z_next = z + h/24.0 * (55*g_vals[i] - 59*g_vals[i-1] + 37*g_vals[i-2] - 9*g_vals[i-3]);

        double x_next = x + h;

        sol.x.push_back(x_next);
        sol.y.push_back(y_next);
        sol.z.push_back(z_next);

        f_vals.push_back(f(x_next, y_next, z_next));
        g_vals.push_back(g(x_next, y_next, z_next));
    }

    return sol;
}

void print_solution(const Solution& sol, const std::string& method_name) {
    cout << "\n" << method_name << ":\n";
    cout << format("{:^10} {:^15} {:^15} {:^15}\n", "x", "y_num", "y_exact", "Error");
    cout << string(55, '-') << "\n";

    for (size_t i = 0; i < sol.x.size(); i++) {
        double exact_val = exact_y(sol.x[i]);
        double error = fabs(sol.y[i] - exact_val);
        cout << format("{:^10.4f} {:^15.8f} {:^15.8f} {:^15.8f}\n",
                      sol.x[i], sol.y[i], exact_val, error);
    }
}

void save_to_file(const Solution& sol, const std::string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << "x,y,y_exact\n";
    for (size_t i = 0; i < sol.x.size(); i++) {
        file << sol.x[i] << "," << sol.y[i] << "," << exact_y(sol.x[i]) << "\n";
    }

    file.close();
}

double calculate_error(const Solution& num, const Solution& exact) {
    if (num.y.size() != exact.y.size()) return -1;

    double max_error = 0;
    for (size_t i = 0; i < num.y.size(); i++) {
        double error = fabs(num.y[i] - exact.y[i]);
        if (error > max_error) {
            max_error = error;
        }
    }

    return max_error;
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
