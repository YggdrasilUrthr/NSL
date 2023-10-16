#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <cmath>
#include <functional>

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

double H_GS(double r, double theta, double phi) {

    double psi = exp(-r) / sqrt(M_PI);
    return pow(psi, 2);

}

double H_2P(double r, double theta, double phi) {

    double psi = 1.0 / 8.0 * sqrt(2.0 / M_PI) * r * exp(-r / 2.0) * cos(theta);
    return pow(psi, 2);

}

int main() {

    const uint32_t M = 10e6;

    std::function<double(double, double, double)> H_prob_density;

    double x_n = 0.0;

    for (size_t i = 0; i < M; ++i) {

        

    }
    
    return 0;

}