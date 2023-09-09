#pragma once
#include <iostream>
#include <string>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "matplotlibcpp.h"
#include <numpy/arrayobject.h>
#include <armadillo>



class Matrix {
    private:
        arma::Mat<double> matrix;
        std::vector<std::vector<double>> matrix_vec; 
    public:
        Matrix(){};
        Matrix(std::vector<std::vector<double>> m);
        static Matrix Matrix_mul(Matrix& A, Matrix& B);
        static Matrix Matrix_mul(std::vector<std::vector<double>>, std::vector<std::vector<double>>);
        static std::vector<double> Vector_scalar(std::vector<double> A, std::vector<double> B);
        static double Euclidean_norm(std::vector<double> A);
        static Matrix transpose(Matrix& A);
        void print();
};