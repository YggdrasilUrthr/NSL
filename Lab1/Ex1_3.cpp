#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <math.h>
#include <array>

const static std::string rangen_lib_path = "../Libs/RanGen/";
const static std::string csv_path = "./CSV/";

double error(double avg, double avg2, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt((avg2 - avg * avg) / n);

    }

}

void init_random_gen(Random &rnd) {

    int seed[4];
    int p1, p2;
    std::ifstream Primes(rangen_lib_path + "Primes");
    
    if (Primes.is_open()){
        
        Primes >> p1 >> p2 ;
   
    } else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
   
    Primes.close();

    std::ifstream input(rangen_lib_path + "seed.in");
    std::string property;
    
    if(input.is_open()){
        
        while (!input.eof() ){
            
            input >> property;
            
            if( property == "RANDOMSEED" ){
            
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            
            }

        }
        
        input.close();
   
    } else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;

    rnd.SaveSeed();

}

double sample_angle(Random &rnd) {

    double x= 0.0;
    double y = 0.0;
    double r = 0.0;

    // Sample angle between zero and pi / 2

    do {

        x = rnd.Rannyu();
        y = rnd.Rannyu();
        r = pow(x, 2) + pow(y, 2);

    } while (r >= 1);
    
    double phi = (y >= 0) ? acos(x / sqrt(r)) : 2 * M_PI -  acos(x / sqrt(r));

    return phi;

}

bool in_interval(double lower, double upper, double num) {

    return (num < upper && num >= lower);

}

int main() {

    Random rnd;
    init_random_gen(rnd);

    const double d = 0.5;
    const double L = 0.1;

    const uint32_t W = 10;                  // Plane width
    const uint32_t H = 10;                  // Plane height

    uint32_t N_lines = W / d;

    const uint32_t M = 10000;
    const uint32_t N = 100;
    uint32_t K = M / N;

    std::vector<uint32_t> pi_avgs;
    std::vector<uint32_t> pi_avgs2;

    for (size_t i = 0; i < N; ++i) {

        uint32_t N_hit = 0;

        for (size_t j = 0; j < K; ++j) {

            double center_x = rnd.Rannyu() * W;
            double center_y = rnd.Rannyu() * H;
            double angle = sample_angle(rnd);



            std::vector<double> upper = {center_x + L / 2.0 * cos(angle), center_y + L / 2.0 * sin(angle)};
            std::vector<double> lower = {center_x - L / 2.0 * cos(angle), center_y - L / 2.0 * sin(angle)};

            for (size_t k = 0; k < N_lines; ++k) {

                if(in_interval(lower[0], upper[0], k * d)) {

                    N_hit++;

                }

            }

        }

        double pi_avg = 2 * L * K / (N_hit * d);
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

    csvwriter writer(csv_path + "Ex1_3.csv");

    for(auto item : data_out) {

        writer.write_csv_line<double>(std::vector<double>(item.begin(), item.end()));

    }    

    return 0;

}