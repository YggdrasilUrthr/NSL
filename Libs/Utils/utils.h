#include <cmath>
#include <vector>
#include <algorithm>

double error(double avg_unif, double avg2_unif, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt(fabs(avg2_unif - avg_unif * avg_unif) / n);

    }

}

inline double sgn(double x) {

    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);

}

inline double cart_dist(std::vector<double> cart_coords) {

    return sqrt(pow(cart_coords[0], 2) + pow(cart_coords[1], 2) + pow(cart_coords[2], 2));

}

std::vector<double> cart_to_sph(std::vector<double> cart_coords) {

    std::vector<double> sph_coord(3,0);

    double r_sqr = 0.0;

    for (size_t i = 0; i < 3; ++i) {

        r_sqr += pow(cart_coords[i], 2);

    }

    sph_coord[0] = sqrt(r_sqr);
    sph_coord[1] = acos(cart_coords[2] / sph_coord[0]);
    sph_coord[2] = sgn(cart_coords[1]) * acos(cart_coords[0] / sqrt(pow(cart_coords[0], 2) + pow(cart_coords[1], 2)));

    return sph_coord;

}