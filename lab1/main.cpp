#include <iostream>
#include "Matrix.h"
#include <armadillo>
#include "Func_image.h"

double f_s(double a) {
    return a;
}

int main() {
    //import_array();
    Matrix m1{{{2.0, 1.0}, {1.0, 2.0}}};
    Matrix m2{{{1.0, 2.0}, {4.0, 5.0}}};
    Matrix::Matrix_mul(m1, m2);
    auto m = Matrix::Matrix_mul({{1.0, 2.0}, {4.0, 5.0}}, {{1.0, 2.0}, {4.0, 5.0}});
    m.print();
    auto m_t = Matrix::transpose(m);
    m_t.print();
    Func_image::plot_image(&f_s);


    std::vector<std::vector<double>> v;
    int M = 10;
    int N = M;
    for (int i = 0; i < M; i++) {
        
        for (size_t j = 0; j < N; j++)
        {
            /* code */
        }
        
    }

    return 0;
}