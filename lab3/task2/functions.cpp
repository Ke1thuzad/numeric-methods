#include "main.h"
#include "../../lab1/task2/main.h"

double find_h(std::vector<double> x, int index) {
    return x[index] - x[index - 1];
}


double* cubic_spline_equation(const std::vector<double> &x, const std::vector<double> &f) {
    const int n = x.size() - 1;

    Matrix<double> tridiag(n - 1, n - 1);
    std::vector<double> d_eq(n - 1);

    {
        double h = find_h(x, 2);
        double prev_h = find_h(x, 1);

        tridiag.coefficients[0][0] = 2 * (prev_h + h);
        tridiag.coefficients[0][1] = h;

        d_eq[0] = 3 * ((f[2] - f[1]) / h - (f[1] - f[0]) / prev_h);
    }

    for (int i = 3; i < n; ++i) {
        double h = find_h(x, i);
        double prev_h = find_h(x, i - 1);

        int index = i - 2;

        tridiag.coefficients[index][index - 1] = prev_h;
        tridiag.coefficients[index][index] = 2 * (prev_h + h);
        tridiag.coefficients[index][index + 1] = h;

        d_eq[index] = 3 * ((f[i] - f[i - 1]) / h - (f[i - 1] - f[i - 2]) / prev_h);
    }

    {
        double h = find_h(x, n);
        double prev_h = find_h(x, n - 1);

        int index = n - 2;

        tridiag.coefficients[index][index - 1] = prev_h;
        tridiag.coefficients[index][index] = 2 * (prev_h + h);

        d_eq[index] = 3 * ((f[n] - f[n - 1]) / h - (f[n - 1] - f[n - 2]) / prev_h);
    }

    return thomas_algorithm(tridiag, d_eq.data());
}

double cubic_spline(const std::vector<double>& x, const std::vector<double>& f, const double point_x) {
    const int n = x.size() - 1;

    double* solution = cubic_spline_equation(x, f);

    std::vector<double> a(n), b(n), c(n), d(n);

    c[0] = 0;

    for (int i = 0; i < n; ++i) {
        a[i] = f[i];

        double h = find_h(x, i + 1);

        c[i + 1] = solution[i];

        b[i] = (f[i + 1] - f[i]) / h - h / 3 * (c[i + 1] + 2 * c[i]);

        d[i] = (c[i + 1] - c[i]) / (3 * h);
    }

    delete[] solution;

    double result;

    std::string out;

    for (int i = 0; i < n; ++i) {
        if (point_x > x[i] && point_x < x[i + 1]) {
            const double line = point_x - x[i];
            result = a[i] + b[i] * line + c[i] * line * line + d[i] * line * line * line;
            out = std::format("{0:.6f} + {1:.6f}(x - {4}) + {2:.6f}(x - {4})^2 + {3:.6f}(x - {4})^3", a[i], b[i], c[i], d[i], x[i]);
            break;
        }
    }

    std::cout << out << std::endl;

    return result;
}
