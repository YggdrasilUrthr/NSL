#pragma once

#include <cmath>
#include <inttypes.h>
#include <vector>
#include <algorithm>

// Header for mutation functions and distance estimators

const static uint32_t MutationsNumber = 4;

// City represented as a point in 2D carthesian space

struct City {

    double x;
    double y;

};

// Distance estimators

inline double L2_Distance(City A, City B);
inline double L1_Distance(City A, City B);

void Permute(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m);           // Permute alleles at position n, m
void N_Shift(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m);           // Shift m contiguous alleles of n positions
void M_Permute(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m);         // Swap first n allels with n alleles starting at position m
void M_Invert(std::vector<uint32_t> &In_Vect, uint32_t m);                      // Invert the first m cities

inline double L2_Distance(City A, City B) {

    return sqrt(pow((A.x - B.x), 2) + pow((A.y - B.y), 2));

}

inline double L1_Distance(City A, City B) {

    return fabs(A.x - B.x) + fabs(A.y - B.y);

}