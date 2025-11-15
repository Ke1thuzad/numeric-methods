#include "main.h"

int main() {
    constexpr double X_prec = 0.8;
    std::vector<double> X{0.2, 0.6, 1.0, 1.4};


    std::cout << "Lagrange Polynomial:" << std::endl;
    Lagrange_polynomial(function, X, X_prec);

    std::cout << std::endl;

    std::cout << "Newton Polynomial" << std::endl;
    Newton_polynomial(function, X, X_prec);

    return 0;
}

