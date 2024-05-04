#include "./Psi.h"
#include "./Ex8.h"

inline double Boltzmann(double x, double beta) {

    return exp(-beta * x);

}

int main() {

    Init();

    csvwriter writer_SA(file_name_SA);
    csvwriter writer_SA_evol(file_name_SA_evol);

    for (size_t i = 0; i < N_SA; ++i) {

        for (size_t j = 0; j < L_SA; ++j) {

            Move_SA();

        }

        if(save_SA) {

            writer_SA.write_csv_line(std::vector<double>({T, energy, energy_err}));
            writer_SA_evol.write_csv_line(std::vector<double>({energy, mu, sigma}));

        }

        //T -= T_step_lin;                      // Linear update law, same weight to high and low temps
        //T *= T_step_pow;                      // Power update law, more time spent at low temps
        //beta = 1.0 / T;                       // Used when a T update law is selected 

        beta += beta_step;                      // Used when a beta update law is selected
        T = 1.0 / beta;                     

        std::cout << "Block number: " << i + 1 << "/" << N_SA << std::endl;
        std::cout << "Acceptance rate: " << (static_cast<double>(accepted_SA) / static_cast<double>(attempted_SA)) << std::endl;
        std::cout << "Temperature: " << T << std::endl;
        std::cout << "Beta: " << beta << std::endl;
        std::cout << "Mu: " << mu << std::endl;
        std::cout << "Sigma: " << sigma << std::endl;
        std::cout << "Energy: " << energy << std::endl;
        std::cout << "--------------------------------" << std::endl;

    }

    csvwriter SA_param_writer(file_name_SA_params);
    SA_param_writer.write_csv_line(std::vector<double>({mu, sigma}));

    if(save_final) {
        
        std::cout << "Measuring energy and filling histogram for final configuration..." << std::endl;
        save = save_final;
        CalcEnergy(N_f, L_f, file_name_E_final);
        Hist();                                 // Create an histogram of psi^2 with final parameters

    }

    return 0;

}

void Input() {

    // Read input informations

    std::ifstream ReadInput;
    ReadInput.open("input.dat");

    // Read SA parameters
    ReadInput >> step_sa;
    ReadInput >> N_SA;
    ReadInput >> L_SA;

    // Read temperature parameters
    ReadInput >> T_i;
    ReadInput >> T_f;

    // Read SA starting parameters
    ReadInput >> mu_0;
    ReadInput >> sigma_0;

    // Read Metropolis paramters
    ReadInput >> step;
    ReadInput >> N;
    ReadInput >> L;

    // Read Metropolis starting parameters
    ReadInput >> x_0;

    // Save files?
    ReadInput >> save;
    ReadInput >> save_SA;
    ReadInput >> save_final;
  
    // Final configuration parameter
    ReadInput >> N_f;
    ReadInput >> L_f;
    ReadInput >> nbins;

    ReadInput.close();

}

void Init() {

    rnd.Init_Random_Gen();

    Input();

    x = x_0;
    mu = mu_0;
    sigma = sigma_0;                                                 

    CalcEnergy(50, 1);                                                      // Equilibrate
    energy = 0;
    energy_err = 0;
    uint32_t accepted = 0;
    uint32_t attempted = 0;

    T = T_i;
    T_step_lin = (T_i - T_f) / N_SA;
    T_step_pow = pow(T_f / T_i, 1.0 / N_SA);                                // T_f = T_i * (x ^ N) => x = (T_f / T_i) ^ (1 / N)
    beta_step = (1.0 / T_f - 1.0 / T_i) / N_SA;                             // Linear update law in beta
    beta = 1.0 / T;

}

void CalcEnergy(uint32_t N, uint32_t L, std::string filename) {

    csvwriter writer(filename);

    double avg = 0.0;
    double avg2 = 0.0;

    for (size_t i = 0; i < N; ++i) {

        uint32_t accepted = 0;
        uint32_t attempted = 0;

        double energy_blk = 0.0;

        for (size_t j = 0; j < L; ++j) {

            double x_mvd = x + step * (rnd.Rannyu() - 0.5) * 2;
            double alpha = std::min(1.0, psi2(x_mvd, mu, sigma) / psi2(x, mu, sigma)); // Sample psi^2

            double r = rnd.Rannyu();

            if (r < alpha) {
                
                x = x_mvd;
                ++accepted;

            }

            energy_blk += (H_psi(x, mu, sigma) / double(L));

            ++attempted;

        }

        avg += energy_blk;
        avg2 += pow(energy_blk, 2.0);

        double iblk = double(i + 1);

        if(save) {

            writer.write_csv_line(std::vector<double>({iblk, energy_blk, avg / iblk, Error(avg, avg2, iblk)}));

        }

        // Only the final value of energy is necessary for SA

        energy = avg / iblk;
        energy_err = Error(avg, avg2, iblk);

        //std::cout << "Acceptance rate: " << (static_cast<double>(accepted) / static_cast<double>(attempted)) << std::endl;

    }

}

void Hist() {

    hist = std::vector<double>(nbins, 0.0);
    double bin_size = (hist_upper - hist_lower) / nbins;

    for (size_t i = 0; i < N_f; ++i) {

        uint32_t accepted = 0;
        uint32_t attempted = 0;

        for (size_t j = 0; j < L_f; ++j) {

            double x_mvd = x + step * (rnd.Rannyu() - 0.5) * 2;
            double alpha = std::min(1.0, psi2(x_mvd, mu, sigma) / psi2(x, mu, sigma)); // Sample psi^2

            double r = rnd.Rannyu();

            if (r < alpha) {
                
                x = x_mvd;
                ++accepted;

            }

            uint32_t idx = floor((x - hist_lower) / bin_size);
            hist[idx]++;

            ++attempted;

        }

    }

    std::vector<double> r_bins(nbins, 0.0);

    for (size_t i = 0; i < nbins; i++) {

        r_bins[i] = hist_lower + bin_size * i;

    }
    
    csvwriter hist_writer(file_name_hist);
    hist_writer.write_csv_line(r_bins);
    hist_writer.write_csv_line(hist);

}

void Move_SA() {

    double prev_mu = mu;
    double prev_sigma = sigma;
    double prev_energy = energy;
    double prev_energy_err = energy_err;

    mu += step_sa * (rnd.Rannyu() - 0.5) * 2 / beta;       
    sigma += step_sa * (rnd.Rannyu() - 0.5) * 2 / beta;             
    
    CalcEnergy();   
    double delta_e = energy - prev_energy;
    double alpha = std::min(1.0, Boltzmann(delta_e, beta));

    double r = rnd.Rannyu();

    if (r < alpha) {

        ++accepted_SA;

    } else {

        mu = prev_mu;
        sigma = prev_sigma;
        energy = prev_energy;
        energy_err = prev_energy_err;

    }

    ++attempted_SA;

}

double Error(double sum, double sum2, double iblk) {

    return sqrt(fabs(sum2 / iblk - pow(sum / iblk, 2)) / iblk);

}