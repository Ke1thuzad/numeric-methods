#include "main.h"

int main(const int argc, char **argv) {
    if (argc == 2) {
        auto read = read_tridiagonal_matrix_from_file(argv[1]);
        const Matrix<double> matrix = read.first;
        double *d = read.second;

        const double* res = thomas_algorithm(matrix, d);



        for (int i = 0; i < matrix.n; ++i) {
            std::cout << res[i] << ", ";
        }

        std::cout << std::endl;

        delete[] res;
        delete[] d;

    } else if (argc > 2) {
        std::cout << "Invalid arguments. Usage: ./main <matrix filepath>" << std::endl;
    }

    return 0;
}

