#include "Matrix.h"


arma::Mat<double> convertNestedStdVectorToArmadilloMatrix(const std::vector<std::vector<double>> &V) {
    arma::Mat<double> A(V.size(), V[0].size());

    for (size_t i = 0; i < V.size(); i++) {
        A.row(i) = arma::conv_to< arma::Row<double> >::from(V[i]);
    }
    return A;
}

Matrix::Matrix(std::vector<std::vector<double>> m) {
    matrix = convertNestedStdVectorToArmadilloMatrix(m);
    matrix_vec = m;
}

Matrix Matrix::Matrix_mul(Matrix& A, Matrix& B) {
    arma::Mat<double> a = A.matrix; //.matrix;//{{1.0, 2}, {3, 4}}; 
    arma::Mat<double> b = B.matrix; //.matrix;//{{1.0, 2}, {3, 4}};

    //std::cout << a*b << std::endl;
    //arma::Mat<double> c = a*b;
    std::cout << a.n_rows << " " << a.n_cols << std::endl;
    
    std::cout << (a*b) << std::endl;
    return Matrix();
}

Matrix Matrix::Matrix_mul(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B) {
    std::vector<std::vector<double>> C;
    for (size_t i = 0; i < A.size(); i++)
    {   
        std::vector<double> line;
        for (size_t j = 0; j < B[0].size(); j++) {
            double sum = 0;
            for (size_t k = 0; k < B.size(); k++) {
                sum += A[i][k] * B[k][j];
            }
            line.push_back(sum);
        }
        C.push_back(line);
    }
    return Matrix(C);
}

void Matrix::print() {
    matrix.print();
}

std::vector<double> Matrix::Vector_scalar(std::vector<double> A, std::vector<double> B) {
    std::vector<double> C;
    for (size_t i = 0; i < A.size(); i++) {
        C.push_back(A[i] * B[i]);
    }
    return C;
}

double Matrix::Euclidean_norm(std::vector<double> A) {
    double sum = 0;
    for (size_t i = 0; i < A.size(); i++) {
        sum += A[i] * A[i];
    }
    return sqrt(sum);
}

Matrix Matrix::transpose(Matrix& A) {
    std::vector<std::vector<double>> c = A.matrix_vec;
    for (size_t i = 0; i < c.size(); i++) {
        for (size_t j = i; j < c[0].size(); j++) 
        {
            c[i][j] = A.matrix_vec[j][i];
        }
        
    }
    return Matrix(c);
}