#include "main.h"
#include <iostream>
#include <cmath>
#include <windows.h>

using namespace std;

int main() {
    SetConsoleOutputCP(CP_UTF8);

    double x0 = 0.0;
    double y0 = 1.0;
    double z0 = 2.0;
    double h = 0.1;
    double x_end = 1.0;

    int steps = static_cast<int>((x_end - x0) / h);

    Solution sol_euler = euler_explicit(x0, y0, z0, h, steps);
    Solution sol_euler_cauchy = euler_cauchy(x0, y0, z0, h, steps);
    Solution sol_euler_improved = euler_improved(x0, y0, z0, h, steps);
    Solution sol_rk4 = runge_kutta_4(x0, y0, z0, h, steps);
    Solution sol_adams = adams_4(x0, y0, z0, h, steps);

    Solution exact;
    exact.x.resize(steps + 1);
    exact.y.resize(steps + 1);
    exact.z.resize(steps + 1);

    for (int i = 0; i <= steps; i++) {
        exact.x[i] = x0 + i * h;
        exact.y[i] = exact_y(exact.x[i]);
        exact.z[i] = exact_z(exact.x[i]);
    }

    print_solution(sol_euler, "Явный метод Эйлера");
    print_solution(sol_euler_cauchy, "Метод Эйлера-Коши");
    print_solution(sol_euler_improved, "Первый улучшенный метод Эйлера");
    print_solution(sol_rk4, "Метод Рунге-Кутты 4-го порядка");
    print_solution(sol_adams, "Метод Адамса 4-го порядка");

    cout << "\n" << string(60, '=') << "\n";
    cout << "Сравнение методов:\n";
    cout << format("{:<40} {:>20}\n", "Метод", "Макс. погрешность");
    cout << string(60, '-') << "\n";

    double err_euler = calculate_error(sol_euler, exact);
    double err_euler_cauchy = calculate_error(sol_euler_cauchy, exact);
    double err_euler_improved = calculate_error(sol_euler_improved, exact);
    double err_rk4 = calculate_error(sol_rk4, exact);
    double err_adams = calculate_error(sol_adams, exact);

    cout << format("{:<40} {:>20.10f}\n", "Явный метод Эйлера", err_euler);
    cout << format("{:<40} {:>20.10f}\n", "Метод Эйлера-Коши", err_euler_cauchy);
    cout << format("{:<40} {:>20.10f}\n", "Улучшенный метод Эйлера", err_euler_improved);
    cout << format("{:<40} {:>20.10f}\n", "Метод Рунге-Кутты 4-го порядка", err_rk4);
    cout << format("{:<40} {:>20.10f}\n", "Метод Адамса 4-го порядка", err_adams);


    cout << "\n" << string(60, '=') << "\n";
    cout << "Оценка погрешности методом Рунге-Ромберга:\n";

    int steps_2h = steps / 2;
    Solution sol_euler_2h = euler_explicit(x0, y0, z0, 2*h, steps_2h);
    Solution sol_euler_cauchy_2h = euler_cauchy(x0, y0, z0, 2*h, steps_2h);
    Solution sol_euler_improved_2h = euler_improved(x0, y0, z0, 2*h, steps_2h);
    Solution sol_rk4_2h = runge_kutta_4(x0, y0, z0, 2*h, steps_2h);
    Solution sol_adams_2h = adams_4(x0, y0, z0, 2*h, steps_2h);

    double rr_euler = runge_romberg_error(sol_euler, sol_euler_2h, 1);
    double rr_euler_cauchy = runge_romberg_error(sol_euler_cauchy, sol_euler_cauchy_2h, 2);
    double rr_euler_improved = runge_romberg_error(sol_euler_improved, sol_euler_improved_2h, 2);
    double rr_rk4 = runge_romberg_error(sol_rk4, sol_rk4_2h, 4);
    double rr_adams = runge_romberg_error(sol_adams, sol_adams_2h, 4);

    cout << format("{:<40} {:>25}\n", "Метод", "Оценка Рунге-Ромберга");
    cout << string(65, '-') << "\n";
    cout << format("{:<40} {:>25.10f}\n", "Метод Эйлера (p=1)", rr_euler);
    cout << format("{:<40} {:>25.10f}\n", "Метод Эйлера-Коши (p=2)", rr_euler_cauchy);
    cout << format("{:<40} {:>25.10f}\n", "Улучшенный метод Эйлера (p=2)", rr_euler_improved);
    cout << format("{:<40} {:>25.10f}\n", "Метод Рунге-Кутты (p=4)", rr_rk4);
    cout << format("{:<40} {:>25.10f}\n", "Метод Адамса (p=4)", rr_adams);

    save_to_file(sol_euler, "euler.csv");
    save_to_file(sol_euler_cauchy, "euler_cauchy.csv");
    save_to_file(sol_euler_improved, "euler_improved.csv");
    save_to_file(sol_rk4, "runge_kutta.csv");
    save_to_file(sol_adams, "adams.csv");
    save_to_file(exact, "exact.csv");

    return 0;
}