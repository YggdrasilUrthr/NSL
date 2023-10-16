#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <cmath>
#include <functional>
#include <iostream>

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

inline double psi(double x, double mu, double sigma) {

    return exp(-pow((x - mu), 2) / (2.0 * pow(sigma, 2))) + exp(-pow((x + mu), 2) / (2.0 * pow(sigma, 2)));

}

inline double psi_sec(double x, double mu, double sigma) {

    double b = psi(x, mu, sigma);
    double c = -psi(x, mu, sigma) / pow(sigma, 2);

    return -psi(x, mu, sigma) / pow(sigma, 2) + 
           pow((x - mu), 2) / pow(sigma, 4) * exp(-pow((x - mu), 2) / (2.0 * pow(sigma, 2))) + 
           pow((x + mu), 2) / pow(sigma, 4) * exp(-pow((x + mu), 2) / (2.0 * pow(sigma, 2)));

}

inline double v_pot(double x) {

    return pow(x, 4) - 5.0 / 2.0 * pow(x, 2);

}

inline double H_psi(double x, double mu, double sigma) {

    return -1.0 / 2.0 * psi_sec(x, mu, sigma) / psi(x, mu, sigma) + v_pot(x);

}

double calc_energy(double mu, double sigma, Random &rnd, const double step) {

    const uint32_t M = 10000;
    const uint32_t N = 100; // Change this
    const uint32_t L = M / N;

    std::vector<double> energy_blk(N, 0);

    for (size_t i = 0; i < N; ++i) {

        uint32_t accepted = 0;
        uint32_t attempted = 0;

        double x = mu;

        for (size_t j = 0; j < L; ++j) {

            double x_mvd = x + step * (rnd.Rannyu() - 0.5) * 2;
            double alpha = std::min(1.0, pow(psi(x_mvd, mu, sigma) / psi(x, mu, sigma), 2));

            double r = rnd.Rannyu();

            if (r <= alpha) {
                
                x = x_mvd;
                ++accepted;

            }

            energy_blk[i] += H_psi(x, mu, sigma);
            ++attempted;

        }

        //std::cout << "Acceptance rate: " << (static_cast<double>(accepted) / static_cast<double>(attempted)) << std::endl;

        energy_blk[i] /= L;

    }
    
    //csvwriter writer(csv_path + "Ex8_E.csv");
    double energy_acc = 0;

    for (size_t i = 0; i < N; ++i) {

        energy_acc += energy_blk[i];
        //writer.write_csv_line(std::vector<double>({i + 1.0, energy_acc / (i + 1)}));

    }

    return energy_acc / N;

}

inline double Boltzmann(double x, double beta) {

    return exp(-beta * x);

}

int main() {

    Random rnd;
    init_random_gen(rnd);

    const double step = 2.0;

    double mu = 1.0;
    double sigma = 1.0;

    //double step_mu = 0.1;
    //double step_sigma = 0.1;

    const double step_sa = 1.5;
    double T_i = 1e-2;
    double T_f = 1e-6;

    const uint32_t N = 20;
    const uint32_t M = 100;

    double energy = calc_energy(mu, sigma, rnd, step);
    uint32_t accepted = 0;
    uint32_t attempted = 0;

    double T = T_i;
    double T_step = (T_i - T_f) / N; 
    double beta = 1.0 / T;

    for (size_t i = 0; i < N; ++i) {

        for (size_t j = 0; j < M; ++j) {

            double mu_mvd;
            double sigma_mvd;

            do {

                mu_mvd = mu + step_sa * (rnd.Rannyu() - 0.5) * 2;
                sigma_mvd = sigma + step_sa * (rnd.Rannyu() - 0.5) * 2;

            } while(mu_mvd < 0 || sigma_mvd <= 0);

            double energy_mvd = calc_energy(mu_mvd, sigma_mvd, rnd, step);

            double x = energy_mvd - energy;
            double alpha = std::min(1.0, Boltzmann(x, beta));

            //std::cout << alpha << "\t";

            double r = rnd.Rannyu();

            if (r <= alpha) {
                
                mu = mu_mvd;
                sigma = sigma_mvd;
                energy = energy_mvd;

                ++accepted;

            }

            ++attempted;

        }

        T -= T_step;
        beta = 1.0 / T;

        std::cout << "Acceptance rate: " << (static_cast<double>(accepted) / static_cast<double>(attempted)) << std::endl;
        std::cout << "Temperature: " << T << std::endl;
        std::cout << "Beta: " << beta << std::endl;
        std::cout << "Mu: " << mu << std::endl;
        std::cout << "Sigma: " << sigma << std::endl;
        std::cout << "Energy: " << energy << std::endl;
        std::cout << "--------------------------------" << std::endl;

    }

    return 0;

}