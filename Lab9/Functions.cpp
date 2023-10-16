#include "Functions.h"

void Permute(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m) {

    std::swap(In_Vect[n], In_Vect[m]);

}

void N_Shift(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m) {

    std::vector<uint32_t> temp_vect(In_Vect.size() + n - 1, 0);
    std::copy(In_Vect.begin() + 1, In_Vect.end(), temp_vect.begin() + n);
    std::copy(temp_vect.begin() + In_Vect.size() - 1, temp_vect.end(), temp_vect.begin());
    std::copy(temp_vect.begin(), temp_vect.begin() + In_Vect.size() - 1, In_Vect.begin() + 1);

}

void M_Permute(std::vector<uint32_t> &In_Vect, uint32_t n, uint32_t m) {

    std::swap_ranges(In_Vect.begin() + 1, In_Vect.begin() + n + 1, In_Vect.begin() + 1 + m + n);

}

void M_Invert(std::vector<uint32_t> &In_Vect, uint32_t m) {

    std::reverse(In_Vect.begin() + 1, In_Vect.begin() + m);

}