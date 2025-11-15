#include "main.h"

int main() {
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0}, y = {0.0, 1.0, 1.4142, 1.7321, 2.0};
    double x_point = 2.0;

    std::cout << "First derivative in point X*: " << compute_first_derivative(x, y, x_point) << std::endl;
    std::cout << "Second derivative in point X*: " << compute_second_derivative(x, y, x_point) << std::endl;

    return 0;
}