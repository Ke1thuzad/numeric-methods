#include "main.h"


int main(const int argc, char **argv) {
    if (argc == 2) {
        const LinearEquation<double> equation = read_equation_from_file(argv[1]);

        std::cout << "Original matrix A:" << std::endl;
        std::cout << equation.matrix;

        std::cout << "Vector b:" << std::endl;
        for (int i = 0; i < equation.n; i++) {
            std::cout << equation.b[i] << " ";
        }
        std::cout << std::endl << std::endl;

        const LUPResult<double> lup = LUPDecomposition(equation.matrix);

        const double* solution = solveLUP(lup, equation.b);

        std::cout << "Solution x:" << std::endl;
        for (int i = 0; i < equation.n; i++) {
            std::cout << "x[" << i << "] = " << solution[i] << std::endl;
        }
        std::cout << std::endl;

        const double det = determinantLUP(lup);
        std::cout << "Determinant of A: " << det << std::endl << std::endl;

        const Matrix<double> inv = inverseMatrix(equation.matrix);
        std::cout << "Inverse matrix A^(-1):" << std::endl;
        std::cout << inv;

        std::cout << "Inverse check: " << std::endl;
        std::cout << equation.matrix * inv;
    } else {
        throw std::runtime_error("Incorrect arguments. Add path to your test file");
    }
}

LinearEquation<double> read_equation_from_file(const char *path) {
    std::ifstream stream(path);
    if (!stream.is_open())
        throw std::runtime_error("Cannot open file");

    int n;
    stream >> n;

    if (n <= 0)
        throw std::domain_error("Matrix size cannot be less than 1");

    LinearEquation<double> equation(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (!(stream >> equation.matrix.coefficients[i][j])) {
                throw std::runtime_error("File format is not correct");
            }
        }
    }

    for (int i = 0; i < n; i++) {
        if (!(stream >> equation.b[i])) {
            throw std::runtime_error("File format is not correct");
        }
    }

    return equation;
}