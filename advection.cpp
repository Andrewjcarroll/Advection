#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include "initial_data.h"
#include "diff.h"
#include "Filter.h"

using namespace std;

void printMatrix(double* matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cout << matrix[i * cols + j] << " ";
        }
        cout << endl;
    }
}

int main() {
    // wave speed
    double c = 1;

    // spatial domain
    double xmin = 0;
    double xmax = 1;

    // time domain
    int m = 1000;  // num of time steps
    double tmin = 0;
    double tmax = 100;

    int n = 200;  // num of grid points

    double kc = 0.88;
    double eps = 0.25;

    // Create a directory called "plots" if it doesn't exist
    system("mkdir -p plots");

    // x grid of n points
    double* X = new double[n]; // Dynamically allocate memory for X
    double dx = (xmax - xmin) / n;
    for (int i = 0; i < n; i++) {
        X[i] = xmin + i * dx;
    }

    // for CFL of 0.1
    double CFL = 0.1;
    double dt = CFL * dx / c;

    // each value of the U array contains the solution for all x values at each timestep
    double** U = new double*[m+1]; // Dynamically allocate memory for U
    for (int i = 0; i <= m; i++) {
        U[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        U[0][i] = initial_u(X[i]);
    }

    // RK stepper
    for (int k = 0; k < m; k++) {
        double t = tmin + k * dt;
        double* k1 = new double[n];
        double* k2 = new double[n];
        double* k3 = new double[n];
        double* k4 = new double[n];
        double* u = new double[n];
        for (int i = 0; i < n; i++) {
            u[i] = U[k][i];
        }
        derivatives(t, u, k1, c, dx, n);
        derivatives(t + 0.5 * dt, u, k2, c, dx, n);
        derivatives(t + 0.5 * dt, u, k3, c, dx, n);
        derivatives(t + dt, u, k4, c, dx, n);

            // Apply filtering to each step
        // Filter(u, k1, kc, eps, n);
        // Filter(u, k2, kc, eps, n);
        // Filter(u, k3, kc, eps, n);
        // Filter(u, k4, kc, eps, n);


        for (int i = 0; i < n; i++) {
            U[k + 1][i] = u[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) * dt / 6;
        }
        
        delete[] k1;
        delete[] k2;
        delete[] k3;
        delete[] k4;
        delete[] u;

        // Plot using gnuplot every 25 time steps
        if (k % 25 == 0) {
            ofstream plotData("plot_data.dat");
            for (int i = 0; i < n; i++) {
                plotData << X[i] << " " << U[k + 1][i] << endl;
            }
            plotData.close();

            // Generate a unique filename based on the step
            ostringstream command;
            command << "gnuplot -e \"set terminal png; set output 'plots/plot_step_" << k << ".png'; plot 'plot_data.dat' using 1:2 with lines title 'Numerical'\"";
            system(command.str().c_str());
        }
    }

    // Free dynamically allocated memory
    delete[] X;
    for (int i = 0; i <= m; i++) {
        delete[] U[i];
    }
    delete[] U;

    return 0;
}
