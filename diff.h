#include <lapacke.h>
#include <cblas.h>
#include "kim.h"
#include <iostream>
#include <cstring>
#include <fstream>
#include <cmath>

#define IDX(i,j) ( (i) + n * (j) )

enum DerivativeType {
    KIM4,
    FINITE_DIFFERENCE,
};

DerivativeType derivativeType = KIM4;

void derivatives(double t, double u[], double uvals[], double c, double dx, int n) {
    switch(derivativeType) {
        case FINITE_DIFFERENCE:
            for (int j = 0; j < n; j++) {
                if (j == 0) {
                    uvals[j] = (c / (2 * dx)) * (u[j + 1] - u[n - 1]);
                } else if (j == n - 1) {
                    uvals[j] = (c / (2 * dx)) * (u[0] - u[j - 1]);
                } else {
                    uvals[j] = (c / (2 * dx)) * (u[j + 1] - u[j - 1]);
                }
            }
            break;


        case KIM4: {
            // Define the matrices P and Q using the Kim method
            double *P = new double[n * n];
            double *Q = new double[n * n];
            derivKim4_PQ(P, Q, n); // Assuming derivKim4_PQ is your function for setting up P and Q

                    // Take the transpose of matrix P
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    std::swap(P[i * n + j], P[j * n + i]);
                }
            }

            // Take the transpose of matrix Q
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    std::swap(Q[i * n + j], Q[j * n + i]);
                }
            }

            // Compute the LU decomposition of the matrix P
            int *ipiv = new int[n];
            int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, P, n, ipiv);

            if (info != 0) {
                std::cerr << "LU factorization failed: " << info << std::endl;
                return;
            }

            // Compute the inverse of the matrix P
            double *Pinv = new double[n * n];
            std::memcpy(Pinv, P, n * n * sizeof(double));
            info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, Pinv, n, ipiv);

            if (info != 0) {
                std::cerr << "Matrix inversion failed: " << info << std::endl;
                return;
            }

            // Compute the product of the inverted matrix Pinv and matrix Q
            double *R = new double[n * n];
            double alpha = 1.0;
            double beta = 0.0;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, Pinv, n, Q, n, beta, R, n);

                            // Multiply R by dx
                for (int i = 0; i < n * n; i++) {
                    R[i] *= 1/dx;
                }
            // Compute the product of matrix R and u
            cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, alpha, R, n, u, 1, beta, uvals, 1);

            // Clean up memory
            delete[] P;
            delete[] Q;
            delete[] ipiv;
            delete[] Pinv;
            delete[] R;
            break;
        }
    }
}
