#include <iostream>
#include <unistd.h>
#include <vector>
#include <future>
#include <matplot/matplot.h>
namespace plt = matplot;


double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

std::vector<std::vector<double>> init_matrix(int N) {
    std::vector<std::vector<double>> out;

    for (int i = 0; i < N; i++) {
        std::vector<double> row;
        for(int j = 0; j < N; j++) {
            row.push_back(0.0);
        }
        out.push_back(row);
    }

    return out;
}

std::vector<std::vector<double>> gener_matrix(int N, double fMin, double fMax) {
    std::vector<std::vector<double>> out = init_matrix(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            out[i][j] = fRand(fMin, fMax);
        }
    }
    return out;
}

std::vector<std::vector<double>> sum_matrix(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    int n = A.size();
    auto C = init_matrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j]+B[i][j];
        }
    }
    return C;
}

std::vector<std::vector<double>> sub_matrix(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    int n = A.size();
    auto C = init_matrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j]-B[i][j];
        }
    }

    return C;
}

std::vector<std::vector<double>> mult_matrix(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    std::vector<std::vector<double>> out;
    int N = A.size();

    out = init_matrix(N);

    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            out[i][j] = 0;
            for(int k = 0; k < N; k++){
                out[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return out;
}

std::vector<std::vector<double>> strass(const std::vector<std::vector<double>>& A,const std::vector<std::vector<double>>& B, int n_min) {
    int n = A.size();
    if (n <= n_min) {
        return mult_matrix(A, B); 
    } else {
        int m = n / 2;
        auto a11 = init_matrix(m);
        auto a12 = init_matrix(m);
        auto a21 = init_matrix(m);
        auto a22 = init_matrix(m);

        auto b11 = init_matrix(m);
        auto b12 = init_matrix(m);
        auto b21 = init_matrix(m);
        auto b22 = init_matrix(m);

        for (int u = 0; u < m; u++) {
            for (int delta = 0; delta < m; delta++) {
                a11[u][delta] = A[u][delta];
                a12[u][delta] = A[u][delta + m];
                a21[u][delta] = A[u+m][delta];
                a22[u][delta] = A[u+m][delta+m];

                b11[u][delta] = B[u][delta];
                b12[u][delta] = B[u][delta + m];
                b21[u][delta] = B[u+m][delta];
                b22[u][delta] = B[u+m][delta+m];
            }
        }

        auto P1 = strass(sum_matrix(a11,a22), sum_matrix(b11, b22), n_min);
        auto P2 = strass(sum_matrix(a21, a22), b11, n_min); 
        auto P3 = strass(a11, sub_matrix(b12, b22), n_min);
        auto P4 = strass(a22, sub_matrix(b21, b11), n_min);
        auto P5 = strass(sum_matrix(a11, a12), b22, n_min);
        auto P6 = strass(sub_matrix(a21, a11), sum_matrix(b11, b12), n_min);
        auto P7 = strass(sub_matrix(a12, a22), sum_matrix(b21, b22), n_min);

        auto C11 = sum_matrix(sub_matrix(sum_matrix(P1, P4), P5), P7);
        auto C12 = sum_matrix(P3, P5);
        auto C21 = sum_matrix(P2, P4);
        auto C22 = sum_matrix(sub_matrix(sum_matrix(P1, P3), P2), P6);

        auto C = init_matrix(n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                C[i][j] = C11[i][j];
                C[i][j+m] = C12[i][j];
                C[i+m][j] = C21[i][j];
                C[i+m][j+m] = C22[i][j];
            }
        }
        return C;
    }
}


std::vector<std::vector<double>> strass_parallel(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B, int n_min) {
    int n = A.size();
    if (n <= n_min) {
        return mult_matrix(A, B); 
    } else {
        int m = n / 2;
        auto a11 = init_matrix(m);
        auto a12 = init_matrix(m);
        auto a21 = init_matrix(m);
        auto a22 = init_matrix(m);

        auto b11 = init_matrix(m);
        auto b12 = init_matrix(m);
        auto b21 = init_matrix(m);
        auto b22 = init_matrix(m);

        for (int u = 0; u < m; u++) {
            for (int delta = 0; delta < m; delta++) {
                a11[u][delta] = A[u][delta];
                a12[u][delta] = A[u][delta + m];
                a21[u][delta] = A[u+m][delta];
                a22[u][delta] = A[u+m][delta+m];

                b11[u][delta] = B[u][delta];
                b12[u][delta] = B[u][delta + m];
                b21[u][delta] = B[u+m][delta];
                b22[u][delta] = B[u+m][delta+m];
            }
        }

        auto task1 = std::async(std::launch::async, strass, sum_matrix(a11,a22), sum_matrix(b11, b22), n_min);
        auto task2 = std::async(std::launch::async, strass, sum_matrix(a21, a22), b11, n_min); 
        auto task3 = std::async(std::launch::async, strass, a11, sub_matrix(b12, b22), n_min);
        auto task4 = std::async(std::launch::async, strass, a22, sub_matrix(b21, b11), n_min);
        auto task5 = std::async(std::launch::async, strass, sum_matrix(a11, a12), b22, n_min);
        auto task6 = std::async(std::launch::async, strass, sub_matrix(a21, a11), sum_matrix(b11, b12), n_min);
        auto task7 = std::async(std::launch::async, strass, sub_matrix(a12, a22), sum_matrix(b21, b22), n_min);

        auto P1 = task1.get(); 
        auto P2 = task2.get(); 
        auto P3 = task3.get(); 
        auto P4 = task4.get(); 
        auto P5 = task5.get(); 
        auto P6 = task6.get(); 
        auto P7 = task7.get(); 

        auto C11 = sum_matrix(sub_matrix(sum_matrix(P1, P4), P5), P7);
        auto C12 = sum_matrix(P3, P5);
        auto C21 = sum_matrix(P2, P4);
        auto C22 = sum_matrix(sub_matrix(sum_matrix(P1, P3), P2), P6);

        auto C = init_matrix(n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                C[i][j] = C11[i][j];
                C[i][j+m] = C12[i][j];
                C[i+m][j] = C21[i][j];
                C[i+m][j+m] = C22[i][j];
            }
        }
        return C;
    }
}



void print_matrix(std::vector<std::vector<double>> A) {
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A.size(); j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
}

void print_vector(std::vector<double> v) {
    for(auto elem: v) {
        std:: cout << elem << " ";
    }
    std::cout << std::endl;
}

int main() {
    int N = 16;
    auto A = gener_matrix(N, -10, 10);
    auto B = gener_matrix(N, -10, 10);

    auto C = strass(A, B, 4);

    print_matrix(C);

    print_matrix(mult_matrix(A, B));


    // Задаем размер матрицы
    int maxN = 1 << 11;

    // Векторы для хранения размеров и времен
    std::vector<int> sizes;
    std::vector<double> standardTimes;
    std::vector<double> strassenTimes;
    std::vector<double> parallelStrassenTimes;

    int n_min = 64;
    for (int N = 2; N <= maxN; N *= 2) {
        sizes.push_back(N);
        std::cout << N << std::endl;
        // Генерируем матрицы A и B

        auto A = gener_matrix(N, -10, 10);
        auto B = gener_matrix(N, -10, 10);

        // Замеряем время для стандартного алгоритма
        auto start = std::chrono::high_resolution_clock::now();
        auto m1 = mult_matrix(A, B);
        auto end = std::chrono::high_resolution_clock::now();
        double standardTime = std::chrono::duration<double>(end - start).count();
        standardTimes.push_back(standardTime);

        // Замеряем время для метода Штрассена
        start = std::chrono::high_resolution_clock::now();
        auto m2 = strass(A, B, n_min);
        end = std::chrono::high_resolution_clock::now();
        double strassenTime = std::chrono::duration<double>(end - start).count();
        strassenTimes.push_back(strassenTime);

        // Замеряем время для метода Штрассена с многопоточностью
        start = std::chrono::high_resolution_clock::now();
        auto m3 = strass_parallel(A, B, n_min);
        end = std::chrono::high_resolution_clock::now();
        double parallelStrassenTime = std::chrono::duration<double>(end - start).count();
        parallelStrassenTimes.push_back(parallelStrassenTime);

        
    }

    // Строим график зависимости времени от размера матрицы
    plt::plot(sizes, standardTimes, sizes, strassenTimes, sizes, parallelStrassenTimes);
    
    plt::xlabel("Matrix Size (N)");
    plt::ylabel("Time (s)");
    auto ax1 = plt::nexttile();
    plt::title("Matrix Multiplication Time Comparison");

    ::matplot::legend(ax1, {"Classic mult", "Strassen", "Parallel Strassen"});


    // // Отображаем график
    plt::save("info.png");
    plt::show();

    print_vector(standardTimes);
    print_vector(strassenTimes);
    print_vector(parallelStrassenTimes);

    return 0;
}

