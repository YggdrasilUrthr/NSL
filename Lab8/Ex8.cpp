#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <cmath>
#include <functional>
#include <iostream>

const static std::string csv_path = "./CSV/";

inline double psi(double x, double mu, double sigma) {

    return exp(-0.5 * pow((x - mu) / sigma, 2)) + exp(-0.5 * pow((x + mu) / sigma, 2));

}

inline double psi2(double x, double mu, double sigma) {

    return pow(psi(x, mu, sigma), 2);

}

inline double psi_sec(double x, double mu, double sigma) {

    //double b = psi(x, mu, sigma);
    //double c = -psi(x, mu, sigma) / pow(sigma, 2);

    return -psi(x, mu, sigma) / pow(sigma, 2) + 
           exp(-pow((x - mu), 2) / (2.0 * pow(sigma, 2))) * pow((x - mu), 2) / pow(sigma, 4) + 
           exp(-pow((x + mu), 2) / (2.0 * pow(sigma, 2))) * pow((x + mu), 2) / pow(sigma, 4);

}

inline double v_pot(double x) {

    return pow(x, 4) - 5.0 / 2.0 * pow(x, 2);

}

inline double H_psi(double x, double mu, double sigma) {

    return -0.5 * psi_sec(x, mu, sigma) / psi(x, mu, sigma) + v_pot(x);

}

double calc_energy(
    
    double mu, double sigma, 
    Random &rnd, const double step, double &x, 
    const uint32_t M = 10000, const uint32_t N = 100,
    bool save = false, std::string file_name = csv_path + "E.csv"
    
) {

    // TODO add M,N as params

    //const uint32_t M = 10000;
    //const uint32_t N = 100; // Change this
    const uint32_t L = M / N;

    std::vector<double> energy_blk(N, 0.0);

    //double x = mu;

    for (size_t i = 0; i < N; ++i) {

        uint32_t accepted = 0;
        uint32_t attempted = 0;

        for (size_t j = 0; j < L; ++j) {

            double x_mvd = x + step * (rnd.Rannyu() - 0.5) * 2;
            double alpha = std::min(1.0, psi2(x_mvd, mu, sigma) / psi2(x, mu, sigma)); // Sample psi^2

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
    
    csvwriter writer(file_name);

    double energy_acc = 0.0;

    for (size_t i = 0; i < N; ++i) {

        energy_acc += energy_blk[i];

        if(save) {

            writer.write_csv_line(std::vector<double>({energy_acc / (i + 1)}));

        }

    }

    return energy_acc / N;

}

inline double Boltzmann(double x, double beta) {

    return exp(-beta * x);

}

int main() {

    Random rnd;
    rnd.Init_Random_Gen();

    const double step = 2.0;
    const double x_0 = 2.0;

    double x = x_0;
    double mu = 1.0; // was 1.0
    double sigma = 1.0; // was 1.0

    //double step_mu = 0.1;
    //double step_sigma = 0.1;

    double step_sa = 0.5;         // Lower steps produced exploding results (rounding errors??)
    double T_i = 1; // was 1e-2
    double T_f = 1e-6; // was 1e-6

    // If sigma gets too small computations wail and the results become meaningless (Dirac delta)

    const uint32_t N = 20;              // N SA steps // was 100
    const uint32_t M = 100;

    std::vector<double> energies(N, 0.0);

    double energy = calc_energy(mu, sigma, rnd, step, x);
    uint32_t accepted = 0;
    uint32_t attempted = 0;

    double T = T_i;
    double T_step_lin = (T_i - T_f) / N;
    double T_step_pow = pow(T_f / T_i, 1.0 / N);                                // T_f = T_i * (x ^ N) => x = (T_f / T_i) ^ (1 / N)
    double beta = 1.0 / T;

    for (size_t i = 0; i < N; ++i) {

        x = x_0;
        calc_energy(mu, sigma, rnd, step, x, 50, 1);                                   // Equilibrate // TODO do this with less steps

        for (size_t j = 0; j < M; ++j) {

            double mu_mvd;
            double sigma_mvd;

            mu_mvd = mu + step_sa * (rnd.Rannyu() - 0.5) * 2;                   // Psi symmetric in mu ?? / beta
            sigma_mvd = sigma + step_sa * (rnd.Rannyu() - 0.5) * 2;             // Psi symmetric in nu ?? / beta
            
            double energy_mvd = calc_energy(mu_mvd, sigma_mvd, rnd, step, x);      

            double delta_e = energy_mvd - energy;
            //double alpha = std::min(1.0, Boltzmann(delta_e, beta));
            double alpha = delta_e > 0 ? Boltzmann(delta_e, beta) : 1.0;

            double r = rnd.Rannyu();

            if (r <= alpha) {
                
                mu = mu_mvd;
                sigma = sigma_mvd;
                energy = energy_mvd;

                ++accepted;

            }

            ++attempted;

        }

        energies[i] = energy;
        // T -= T_step_lin;                     // Linear update law, same weight to high and low temps
        T *= T_step_pow;                        // Power update law, more time spent at low temps
        beta = 1.0 / T;

        std::cout << "Acceptance rate: " << (static_cast<double>(accepted) / static_cast<double>(attempted)) << std::endl;
        std::cout << "Temperature: " << T << std::endl;
        std::cout << "Beta: " << beta << std::endl;
        std::cout << "Mu: " << mu << std::endl;
        std::cout << "Sigma: " << sigma << std::endl;
        std::cout << "Energy: " << energy << std::endl;
        std::cout << "--------------------------------" << std::endl;

    }

    //calc_energy(mu, sigma, rnd, step, x, 10000, 100, true, csv_path + "E_final.csv");

    csvwriter writer(csv_path + "Ex8_E.csv");
    
    double idx = 1; // TODO fix this in a better way

    for(auto val : energies) {

        writer.write_csv_line(std::vector<double>({idx, val}));
        ++idx;

    }

    return 0;

}