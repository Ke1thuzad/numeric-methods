#include "main.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <windows.h>

int main() {
    SetConsoleOutputCP(CP_UTF8);

    double a = 0.0;
    double b = 2.0;
    double h = 0.01;
    double bc_left = 0.0;
    double bc_right = 1.0;
    int steps = static_cast<int>((b - a) / h);

    double eps = 1e-8;
    double eta = shooting_method(a, b, h, bc_left, bc_right, eps);

    std::cout << "y(0) = " << eta << std::endl;

    Solution shooting_sol = runge_kutta_4_system(a, eta, bc_left, h, steps);

    int N = 200;
    BoundarySolution fd_sol = finite_difference_method(a, b, N);

    double abs_error_shooting = compute_absolute_error(shooting_sol.x, shooting_sol.y);
    double abs_error_fd = compute_absolute_error(fd_sol.x, fd_sol.y);

    BoundarySolution fd_sol_h = finite_difference_method(a, b, N);
    BoundarySolution fd_sol_2h = finite_difference_method(a, b, N/2);
    double rr_error_fd = runge_romberg_error(fd_sol_h, fd_sol_2h, 2);

    const Solution& shooting_sol_h = shooting_sol;
    Solution shooting_sol_2h = runge_kutta_4_system(a, eta, bc_left, h*2, steps/2);
    double rr_error_shooting = runge_romberg_error(shooting_sol_h, shooting_sol_2h, 4);

    std::cout << "Метод             | Абс погрешность | Рунге-Ромберг\n";
    std::cout << "---------------------------------------------------\n";
    std::cout << std::left << std::setw(31) << "Метод стрельбы"
              << "| " << std::setw(16) << abs_error_shooting
              << "| " << rr_error_shooting << "\n";
    std::cout << std::left << std::setw(30) << "Конечно-разностный"
              << "| " << std::setw(16) << abs_error_fd
              << "| " << rr_error_fd << "\n\n";

    save_solution_to_file(fd_sol, "finite_difference.csv");

    std::ofstream shooting_file("shooting.csv");
    if (shooting_file.is_open()) {
        shooting_file << "x,y,y_exact,error\n";
        for (size_t i = 0; i < shooting_sol.x.size(); i++) {
            double exact_val = exact_y(shooting_sol.x[i]);
            double error_val = fabs(shooting_sol.y[i] - exact_val);
            shooting_file << shooting_sol.x[i] << ","
                         << shooting_sol.y[i] << ","
                         << exact_val << ","
                         << error_val << "\n";
        }
        shooting_file.close();
    }

    return 0;
}