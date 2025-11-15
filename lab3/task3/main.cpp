#include "main.h"

int main() {
    const std::vector<double> x{0.1, 0.5, 0.9, 1.3, 1.7, 2.1};
    const std::vector<double> y{-2.3026, -0.69315, -0.10536, 0.26236, 0.53063, 0.74194};

    std::cout << "First degree approximation polynomial:" << std::endl;
    approximation_polynomial_first_degree(x, y);

    std::cout << std::endl;

    std::cout  << "Second degree approximation polynomial:" << std::endl;
    approximation_polynomial_second_degree(x, y);

    return 0;
}