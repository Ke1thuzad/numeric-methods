#include "main.h"
#include <format>

int main() {
    const std::vector<double> x = {0.1, 0.5, 0.9, 1.3, 1.7};
    const std::vector<double> f = {-2.3026, -0.69315, -0.10536, 0.26236, 0.53063};
    constexpr double x_point = 1.5;

    std::cout << std::format("f({:.1f}) = {:.6f}", x_point, cubic_spline(x, f, x_point)) << std::endl;

    return 0;
}