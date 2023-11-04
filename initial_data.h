#ifndef INITIAL_DATA_H
#define INITIAL_DATA_H

#include <cmath>

enum InitialConditionType {
    KIM,
    GAUSSIAN,
    SIN
};

InitialConditionType initialConditionType = GAUSSIAN;

double initial_u(double x) {
    double mean, variance;
    switch(initialConditionType) {
        case KIM:
            return sin(17 * M_PI * x) + cos(29 * M_PI * x);
        case GAUSSIAN:
            mean = 0.5;
            variance = 0.1;
            return exp(-(x - mean)*(x - mean) / (2 * variance * variance)) / (sqrt(2 * M_PI) * variance);
        case SIN:
            return sin(M_PI*x);
    }
    return 0.0; // Default return value if the initial condition type is not recognized
}

#endif // INITIAL_DATA_H
