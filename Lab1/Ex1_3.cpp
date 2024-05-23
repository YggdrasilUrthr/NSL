#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <math.h>
#include <array>

const static std::string csv_path = "./CSV/";

double error(double avg, double avg2, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt((avg2 - avg * avg) / n);

    }

}

std::vector<double> sample_vertex(Random &rnd) {

    double x= 0.0;
    double y = 0.0;
    double r2 = 0.0;

    // Sample point on unitary circumference

    do {

        x = rnd.Rannyu() * 2.0 - 1.0;
        y = rnd.Rannyu() * 2.0 - 1.0;
        r2 = pow(x, 2) + pow(y, 2);

    } while (r2 >= 1);

    return std::vector<double>({x / sqrt(r2), y / sqrt(r2)});

}

int main() {

    Random rnd;
    rnd.Init_Random_Gen();

    const double d = 0.5;                   // Line distance
    const double L = 0.1;                   // Needle length

    // The sistem is periodic, is then possible to check hits only on a single cell (0, d)

    const uint32_t M = 100000;              // Number of throws
    const uint32_t N = 100;                 // Number of blocks
    uint32_t K = M / N;

    std::vector<double> pi_avgs;
    std::vector<double> pi_avgs2;

    csvwriter writer(csv_path + "Ex1_3_pos.csv");

    for (size_t i = 0; i < N; ++i) {

        uint32_t N_hit = 0;

        for (size_t j = 0; j < K; ++j) {

            double center_x = rnd.Rannyu() * d; // Lines are infinite => only one center coordinate (always inside cell) is needed
            std::vector<double> vertex1 = sample_vertex(rnd);

            writer.write_csv_line<double>(vertex1); 

            std::vector<double> vertex2 = {center_x + L * vertex1[0], L * vertex1[1]};      // This can be inside or outside cell

            // Lines line in x = 0 and x = d => needle cross line if lower_x < 0 v upper_x > d

            if(vertex2[0] > d || vertex2[0] < 0) {

                N_hit++;

            }

        }

        double pi_avg = 2.0 * L * (double)K / (N_hit * d);
        double pi_avg2 = pow(pi_avg, 2);

        pi_avgs.push_back(pi_avg);
        pi_avgs2.push_back(pi_avg2);

    }
    
    std::vector<std::array<double, 2>> sum_prog_vect;
    std::vector<std::array<double, 2>> data_out;

    for (size_t i = 0; i < N; ++i) {

        sum_prog_vect.push_back({0.0, 0.0});
        data_out.push_back({0.0, 0.0});

        for (size_t j = 0; j < i + 1; ++j) {

            sum_prog_vect[i][0] += pi_avgs[j];
            sum_prog_vect[i][1] += pi_avgs2[j];

        }

        sum_prog_vect[i][0] /= (i + 1);
        sum_prog_vect[i][1] /= (i + 1);
        
        data_out[i][0] = sum_prog_vect[i][0];
        data_out[i][1] = error(sum_prog_vect[i][0], sum_prog_vect[i][1], i);

    }

    writer.change_file(csv_path + "Ex1_3.csv");

    for(auto item : data_out) {

        writer.write_csv_line<double>(std::vector<double>(item.begin(), item.end()));

    }    

    return 0;

}