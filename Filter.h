#include "Kim_Filter.h" // Include the header file where kim_filter_PQ is defined
#include <iostream>
void transposeMatrix(double *matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            std::swap(matrix[i * n + j], matrix[j * n + i]);
        }
    }
}
void Filter(double u[], double ufilter[], double kc, double eps, int n) {
    // Define the matrices P and Q using the Kim method
    double *P = new double[n * n];
    double *Q = new double[n * n];
    kim_filter_PQ(P, Q, kc, eps, n); // Assuming kim_filter_PQ is your function for setting up P and Q

    transposeMatrix(P, n);
    transposeMatrix(Q, n);

    // Compute the LU decomposition of the matrix P
    int *ipiv = new int[n];
    int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, P, n, ipiv);

    if (info != 0) {
        std::cerr << "LU factorization failed: " << info << std::endl;
        delete[] P;
        delete[] Q;
        delete[] ipiv;
        return;
    }

    // Compute the inverse of the matrix P
    double *Pinv = new double[n * n];
    std::memcpy(Pinv, P, n * n * sizeof(double));
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, Pinv, n, ipiv);

    if (info != 0) {
        std::cerr << "Matrix inversion failed: " << info << std::endl;
        delete[] P;
        delete[] Q;
        delete[] ipiv;
        delete[] Pinv;
        return;
    }

    // Compute the product of the inverted matrix Pinv and matrix Q
    double *RF = new double[n * n];
    double alpha = 1.0;
    double beta = 0.0;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, Pinv, n, Q, n, beta, RF, n);


    // Compute the product of matrix R and u
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, alpha, RF, n, u, 1, beta, ufilter, 1);

    // Clean up memory
    delete[] P;
    delete[] Q;
    delete[] ipiv;
    delete[] Pinv;
    delete[] RF;
}
