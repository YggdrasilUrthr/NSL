#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <cmath>
#include <functional>
#include <algorithm>

const static std::string csv_path = "./CSV/";

inline double sgn(double x) {

    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);

}

std::vector<double> cart_to_sph(std::vector<double> cart_coords) {

    std::vector<double> sph_coord(3,0);

    double r_sqr = 0.0;

    for (size_t i = 0; i < 3; ++i) {

        r_sqr += pow(cart_coords[i], 2);

    }
    
    // TODO check divergences??

    sph_coord[0] = sqrt(r_sqr);
    sph_coord[1] = acos(cart_coords[2] / sph_coord[0]);
    sph_coord[2] = sgn(cart_coords[1]) * acos(cart_coords[0] / sqrt(pow(cart_coords[0], 2) + pow(cart_coords[1], 2)));

    return sph_coord;

}

double H_GS(std::vector<double> sph_coords) {

    double psi = exp(-sph_coords[0]) / sqrt(M_PI);
    return pow(psi, 2);

}

double H_2P(std::vector<double> sph_coords) {

    double psi = 1.0 / 8.0 * sqrt(2.0 / M_PI) * sph_coords[0] * exp(-sph_coords[0] / 2.0) * cos(sph_coords[1]);
    return pow(psi, 2);

}

bool Move(Random &rnd, double metr_step, std::vector<double> &x_n, std::function<double(std::vector<double>)> psi) {

    std::vector<double> x_n_m(3,0);

    for (size_t j = 0; j < 3; j++) {

        x_n_m[j] = x_n[j] + (rnd.Rannyu() * 2.0 - 1.0) * metr_step;

    }

    std::vector<double> x_n_sph = cart_to_sph(x_n);
    std::vector<double> x_n_m_sph = cart_to_sph(x_n_m);

    double p = std::min(1.0, psi(x_n_m_sph) / psi(x_n_sph));
    double r = rnd.Rannyu();

    if (r < p) {
    
        x_n = x_n_m;
        return true;
        
    }

    return false;

}

int main() {

    const uint32_t M = 10e4;            // Final goal 10e6
    double metr_step = 1.1;

    Random rnd;
    rnd.Init_Random_Gen();

    csvwriter writer(csv_path + "Ex5_1_GS.csv");

    std::vector<double> x_n(3, 0);
    uint32_t attempted = 0;
    uint32_t accepted = 0;

    for (size_t i = 0; i < M; ++i) {

        if (Move(rnd, metr_step, x_n, H_GS)) {

            accepted++;
            writer.write_csv_line(x_n);

        }

        attempted++;

    }

    std::cout << "Acceptance ratio GS: " << (double)accepted / attempted << std::endl;

    writer.change_file(csv_path + "Ex5_1_2P.csv");

    x_n = {0.5, 0, 0};
    metr_step = 2.2;
    attempted = 0;
    accepted = 0;

    for (size_t i = 0; i < M; ++i) {

        if (Move(rnd, metr_step, x_n, H_2P)) {

            accepted++;
            writer.write_csv_line(x_n);

        }

        attempted++;

    }

    std::cout << "Acceptance ratio 2P: " << (double)accepted / attempted << std::endl;

    return 0;

}