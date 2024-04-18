#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <array>
#include <math.h>

//#define DUMP_POLAR_POS

const static std::string csv_path = "./CSV/";

#ifdef DUMP_POLAR_POS

csvwriter pos_writer_polar(csv_path + "Ex2_2_continue_polar.csv");

#endif

bool in_interval(double lower, double upper, double num) {

    return (num < upper && num >= lower);

}

void move_discr(Random &rnd, std::array<double, 3> &position, double a) {

    double r = rnd.Rannyu();

    // X, -X, Y, -Y, Z, -Z

    double step = 1.0 / 6.0;

    for (size_t i = 0; i < 6; ++i) {

        if(in_interval(i * step, (i + 1) * step, r)) {

            int fw_bw = (i % 2) ? -1.0 : 1.0; 
            size_t idx = ceil((i + 1) / 2.0) - 1.0;

            position[idx] += fw_bw * a;

            break;

        }

    }
    
}

std::vector<double> sample_polar(Random &rnd) {

    // See https://mathworld.wolfram.com/SpherePointPicking.html

    std::vector<double> polar_vect;

    double u = rnd.Rannyu();
    double v = rnd.Rannyu();

    double theta = acos(2.0 * u - 1.0);
    polar_vect.push_back(theta);

    double phi = 2.0 * M_PI * v;
    polar_vect.push_back(phi);

#ifdef DUMP_POLAR_POS

    pos_writer_polar.write_csv_line(polar_vect);

#endif

    return polar_vect;

}

void move_cont(Random &rnd, std::array<double, 3> &position, double a) {

    std::vector<double> polar_vect = sample_polar(rnd);

    position[0] += a * sin(polar_vect[0]) * cos(polar_vect[1]);
    position[1] += a * sin(polar_vect[0]) * sin(polar_vect[1]);
    position[2] += a * cos(polar_vect[0]);

}

double norm2(std::array<double, 3> position) {

    double norm2 = 0;

    for (auto coord : position) {

        norm2 += pow(coord, 2);

    }

    return norm2;

}

double error(double avg, double avg2, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt((avg2 - avg * avg) / n);

    }

}

int main() {

    Random rnd;
    rnd.Init_Random_Gen();

    const uint32_t M = 10000;       // Number of random walks
    const uint32_t N = 100;         // Number of steps per walker
    const uint32_t B = 100;         // Number of blocks
    const uint32_t L = M/N;         // Number of walkers per block
    
    const uint32_t a = 1;

    std::vector<std::vector<double>> r_n2_discr(B, std::vector<double>(N, 0));
    std::vector<std::vector<double>> r_n2_cont(B, std::vector<double>(N, 0));

    for (size_t i = 0; i < B; ++i) {

        for (size_t j = 0; j < L; ++j) {

            std::array<double, 3> position_discr = {0.0, 0.0, 0.0};
            std::array<double, 3> position_cont = {0.0, 0.0, 0.0};

            for (size_t k = 0; k < N; ++k) {

                move_discr(rnd, position_discr, a);
                r_n2_discr[i][k] += norm2(position_discr) / double(L); // Divide by L now to save another for loop

                move_cont(rnd, position_cont, a);
                r_n2_cont[i][k] += norm2(position_cont) / double(L); // Divide by L now to save another for loop

            }

        }

    }

    // Calculate errors and means over blocks at fixed step

    std::vector<std::vector<double>> r_n_discr_data(N, std::vector<double>(2, 0));
    std::vector<double> r_n_discr_avg2(N, 0);

    std::vector<std::vector<double>> r_n_cont_data(N, std::vector<double>(2, 0));
    std::vector<double> r_n_cont_avg2(N, 0);

    for (size_t i = 0; i < N; ++i) {

        for (size_t j = 0; j < B; ++j) {

            r_n_discr_data[i][0] += r_n2_discr[j][i] / double(B);
            r_n_discr_avg2[i] += pow(r_n2_discr[j][i],2) / double(B);

            r_n_cont_data[i][0] += r_n2_cont[j][i] / double(B);
            r_n_cont_avg2[i] += pow(r_n2_cont[j][i],2) /double(B);

        }

        r_n_discr_data[i][1] = i != 0 ? error(r_n_discr_data[i][0], r_n_discr_avg2[i], B) : 0.0;   // This is the error on r_n^2, skip i = 0
        r_n_discr_data[i][0] = sqrt(r_n_discr_data[i][0]);
        r_n_discr_data[i][1] /= (2 * r_n_discr_data[i][0]);                         // Propagate errors to sqrt

        r_n_cont_data[i][1] = i != 0 ? error(r_n_cont_data[i][0], r_n_cont_avg2[i], B) : 0.0;   // This is the error on r_n^2, skip i = 0
        r_n_cont_data[i][0] = sqrt(r_n_cont_data[i][0]);
        r_n_cont_data[i][1] /= (2 * r_n_cont_data[i][0]);                         // Propagate errors to sqrt

    }
    

    csvwriter writer(csv_path + "Ex2_2_discrete.csv");

    for (auto item : r_n_discr_data) {

        writer.write_csv_line(std::vector<double>{item});

    }

    writer.change_file(csv_path + "Ex2_2_continue.csv");

    for (auto item : r_n_cont_data) {

        writer.write_csv_line(std::vector<double>{item});

    }
    

    return 0;

}