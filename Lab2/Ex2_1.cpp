#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <array>
#include <functional>
#include <math.h>

const static std::string rangen_lib_path = "../Libs/RanGen/";
const static std::string csv_path = "./CSV/";

//TODO insert tis function in library
//Check non uniform integral result

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

double error(double avg_unif, double avg2_unif, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt((avg2_unif - avg_unif * avg_unif) / n);

    }

}

// generalized accept-reject

double ran_ar(Random &rnd, std::function<double(double)> distrib, double xmin, double xmax, double pmax) {

    double r;
    double x;

    do {

        x = xmin + (xmax - xmin) * rnd.Rannyu();
        r = rnd.Rannyu();

    } while (r >= (distrib(x) / pmax));
    
    return x;

}

int main() {

    Random rnd;
    init_random_gen(rnd);

    const uint32_t M = 10000;
    const uint32_t N = 100;

    uint32_t L = M / N;

    std::vector<double> avg_unif;
    std::vector<double> avg2_unif;
    std::vector<double> avg_quad;
    std::vector<double> avg2_quad;

    std::function<double(double)> quadratic = [](double x) {return (-pow(x, 2) + 1.0) * (3.0 / 2.0);}; //integ 2 / 3

    csvwriter writer(csv_path + "Ex2_1_quadratic_distr.csv");

    for (size_t i = 0; i < N; ++i) {

        double sum_unif = 0;
        double sum_quad = 0;

        for (size_t j = 0; j < L; ++j) {
        
            double x_unif = rnd.Rannyu();
            double x_quad = ran_ar(rnd, quadratic, 0, 1, (3.0 / 2.0));

            writer.write_csv_line(std::vector<double>({x_quad}));

            sum_unif += cos(M_PI * x_unif / 2);
            sum_quad += (cos(M_PI * x_quad / 2) / quadratic(x_quad)); 

        }

        avg_unif.push_back(sum_unif * (M_PI / 2) / L);
        avg2_unif.push_back(pow(avg_unif[i], 2));
        avg_quad.push_back(sum_quad * (M_PI / 2) / L);
        avg2_quad.push_back(pow(avg_quad[i], 2));

    }
    
    std::vector<std::array<double, 2>> G_N_unif;
    std::vector<std::array<double, 2>> G_N_quad;
    double sum_prog_unif = 0;
    double sum_prog2_unif = 0;
    double sum_prog_quad = 0;
    double sum_prog2_quad = 0;

    for (size_t i = 0; i < N; ++i) {
        
        sum_prog_unif += avg_unif[i];
        sum_prog2_unif += avg2_unif[i];
        sum_prog_quad += avg_quad[i];
        sum_prog2_quad += avg2_quad[i];

        G_N_unif.push_back({sum_prog_unif / (i + 1), error(sum_prog_unif / (i + 1), sum_prog2_unif / (i + 1), i)});
        G_N_quad.push_back({sum_prog_quad / (i + 1), error(sum_prog_quad / (i + 1), sum_prog2_quad / (i + 1), i)});

    }
    

    writer.change_file(csv_path + "Ex2_1_uniform.csv");

    for (auto item : G_N_unif) {

        writer.write_csv_line(std::vector<double>(item.begin(), item.end()));

    }

    writer.change_file(csv_path + "Ex2_1_quad.csv");

    for (auto item : G_N_quad) {

        writer.write_csv_line(std::vector<double>(item.begin(), item.end()));

    }

    return 0;

}