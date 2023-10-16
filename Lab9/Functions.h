#pragma once

#include <cmath>
#include <inttypes.h>
#include <vector>
#include <algorithm>

const static uint32_t MutationsNumber = 4;

struct City {

    double x;
    double y;

};

inline double L2_Distance(City A, City B);
inline double L1_Distance(City A, City B);

void Permute(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m);
void N_Shift(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m);
void M_Permute(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m);
void M_Invert(std::vector<uint32_t> &In_Vect, uint32_t m);

inline double L2_Distance(City A, City B) {

    return sqrt(pow((A.x - B.x), 2) + pow((A.y - B.y), 2));

}

inline double L1_Distance(City A, City B) {

    return fabs(A.x - B.x) + fabs(A.y - B.y);

}