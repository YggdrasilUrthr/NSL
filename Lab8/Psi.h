#pragma once

#include <cmath>

// This file contains all the QM inline functions

inline double psi(double x, double mu, double sigma) {

    return exp(-0.5 * pow((x - mu) / sigma, 2)) + exp(-0.5 * pow((x + mu) / sigma, 2));

}

inline double psi2(double x, double mu, double sigma) {

    return pow(psi(x, mu, sigma), 2);

}

inline double psi_sec(double x, double mu, double sigma) {

    return -psi(x, mu, sigma) / pow(sigma, 2) + 
           exp(-pow((x - mu), 2) / (2.0 * pow(sigma, 2))) * pow((x - mu), 2) / pow(sigma, 4) + 
           exp(-pow((x + mu), 2) / (2.0 * pow(sigma, 2))) * pow((x + mu), 2) / pow(sigma, 4);

}

inline double v_pot(double x) {

    return pow(x, 4) - 5.0 / 2.0 * pow(x, 2);

}

inline double H_psi(double x, double mu, double sigma) {

    return -0.5 * psi_sec(x, mu, sigma) / psi(x, mu, sigma) + v_pot(x);

}