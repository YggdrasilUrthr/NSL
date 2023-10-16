#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <array>
#include <math.h>

const static std::string rangen_lib_path = "../Libs/RanGen/";
const static std::string csv_path = "./CSV/";

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

    std::vector<double> polar_vect;

    double x_theta = 0.0;
    double y_theta = 0.0;
    double r_theta = 0.0;
    
    do {

        x_theta = -1.0 + 2.0 * rnd.Rannyu();
        y_theta = rnd.Rannyu();
        r_theta = pow(x_theta, 2) + pow(y_theta, 2);

    } while (r_theta >= 1);
    
    double theta = acos(x_theta / sqrt(r_theta));
    polar_vect.push_back(theta);

    double x_phi = 0.0;
    double y_phi = 0.0;
    double r_phi = 0.0;

    do {

        x_phi = -1.0 + 2.0 * rnd.Rannyu();
        y_phi = -1.0 + 2.0 * rnd.Rannyu();
        r_phi = pow(x_phi, 2) + pow(y_phi, 2);

    } while (r_phi >= 1);
    
    double phi = (y_phi >= 0) ? acos(x_phi / sqrt(r_phi)) : 2 * M_PI -  acos(x_phi / sqrt(r_phi));
    polar_vect.push_back(phi);

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

int main() {

    Random rnd;
    init_random_gen(rnd);

    const uint32_t M = 10000;
    const uint32_t N = 100;
    const uint32_t a = 1;

    std::vector<double> r_n_discr(100, 0);
    std::vector<double> r_n_cont(100, 0);

    for (size_t i = 0; i < M; ++i) {

        std::array<double, 3> position_discr = {0.0, 0.0, 0.0};
        std::array<double, 3> position_cont = {0.0, 0.0, 0.0};

        for (size_t j = 0; j < N; ++j) {

            move_discr(rnd, position_discr, a);
            r_n_discr[j] += norm2(position_discr);

            move_cont(rnd, position_cont, a);
            r_n_cont[j] += norm2(position_cont);

        }

    }

    csvwriter writer(csv_path + "Ex2_2_discrete.csv");

    for (auto item : r_n_discr) {

        writer.write_csv_line(std::vector<double>{item / M});

    }

    writer.change_file(csv_path + "Ex2_2_continue.csv");

    for (auto item : r_n_cont) {

        writer.write_csv_line(std::vector<double>{item / M});

    }
    

    return 0;

}