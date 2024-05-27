#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"
#include "../Libs/Utils/utils.h"

#include <cmath>
#include <functional>
#include <algorithm>

const static std::string csv_path = "./CSV/";

// Ground state Psi

double H_GS(std::vector<double> sph_coords) {

    double psi = exp(-sph_coords[0]) / sqrt(M_PI);
    return pow(psi, 2);

}

// Excited state Psi

double H_2P(std::vector<double> sph_coords) {

    double psi = 1.0 / 8.0 * sqrt(2.0 / M_PI) * sph_coords[0] * exp(-sph_coords[0] / 2.0) * cos(sph_coords[1]);
    return pow(psi, 2);

}

// Perform a Metropolis move

bool Move(
    
    Random &rnd, 
    double metr_step, 
    std::vector<double> &x_n, 
    std::function<double(std::vector<double>)> psi2, 
    bool gauss_sample = false
    
) {

    std::vector<double> x_n_m(3,0);

    for (size_t j = 0; j < 3; j++) {        

        if(gauss_sample) {

            // Use metr_step as std. dev. => implemented version in Random = my implementation, see Ex3_1_1
            x_n_m[j] = x_n[j] + rnd.Gauss(0, metr_step);     

        } else {

            x_n_m[j] = x_n[j] + rnd.Rannyu(-metr_step, metr_step);

        }

    }

    std::vector<double> x_n_sph = cart_to_sph(x_n);
    std::vector<double> x_n_m_sph = cart_to_sph(x_n_m);

    double p = std::min(1.0, psi2(x_n_m_sph) / psi2(x_n_sph));
    double r = rnd.Rannyu();

    if (r < p) {
    
        x_n = x_n_m;
        return true;
        
    }

    return false;

}

// Run the simulation

void Run_Sim(
    
    const uint32_t n_blocks, 
    const uint32_t step_block, 
    Random &rnd, double metr_step, 
    std::function<double(std::vector<double>)> psi2,
    csvwriter &writer_r,
    csvwriter *writer_pos = nullptr,
    std::vector<double> start_x = std::vector<double>({3, 0}),
    bool gauss_sample = false

) {

    std::vector<double> x_n = start_x;              // Metropolis starting point
    std::vector<double> prog_avg(3, 0);

    for (size_t i = 0; i < n_blocks; ++i) {

        uint32_t attempted = 0;
        uint32_t accepted = 0;

        double r = 0.0;

        for (size_t j = 0; j < step_block; ++j) {

            if (Move(rnd, metr_step, x_n, psi2, gauss_sample)) {

                accepted++;

                if (writer_pos != nullptr) {

                    writer_pos->write_csv_line(x_n);

                }

            }

            r += (cart_dist(x_n) / step_block);
            attempted++;
        
        }

        prog_avg[0] += r;
        prog_avg[1] += pow(r, 2);

        prog_avg[2] = error(prog_avg[0] / (i + 1), prog_avg[1] / (i + 1), i);       

        writer_r.write_csv_line(std::vector<double>({prog_avg[0] / (i + 1), prog_avg[2]})); 

        std::cout << "Acceptance ratio in block " + std::to_string(i + 1) + ": " << (double)accepted / attempted << std::endl;

    }

}

// Run Equilibration (same as Run_Sim but without data-blocking)

std::vector<double> Equilibrate(
    
    Random &rnd, double metr_step, 
    std::function<double(std::vector<double>)> psi2,
    double start_r = 50.0,                                   // Completely random large value
    csvwriter *writer_r = nullptr,
    bool gauss_sample = false

) {

    const uint32_t N_eq = 1e3;
    const uint32_t L_eq = 1;

    std::vector<double> x_n = std::vector<double>({start_r, 0.0, 0.0});    // Start far from the origin to see equilibration effect

    for (size_t i = 0; i < N_eq; ++i) {

        uint32_t attempted = 0;
        uint32_t accepted = 0;

        double r = 0.0;

        for (size_t j = 0; j < L_eq; ++j) {

            Move(rnd, metr_step, x_n, psi2, gauss_sample);
            r += (cart_dist(x_n) / L_eq);
        
        }    

        if (writer_r != nullptr) {

            writer_r->write_csv_line(std::vector<double>({r})); 

        }

    }

    return x_n;

}

int main(int argc, char ** argv) {

    const uint32_t M = 1e6;        // Number of throws
    const uint32_t N = 100;         // Number of blocks
    const uint32_t L = M / N;       // Throws per block => 10e4 throws per block
    double metr_step_gs = 1.2; 
    double metr_step_2P = 3.0;

    Random rnd;
    rnd.Init_Random_Gen();

    std::vector<double> eq_x_GS;
    std::vector<double> eq_x_2P;

    bool gauss_sample = false;
    std::string gauss_fn = "";

    // Second CL argument: Perform Gauss sampling??
    if (argc > 2 && std::stoi(std::string(argv[2]))) {
    
        gauss_sample = true;
        gauss_fn = "_gauss";

        // Change metr_step to have 50% prob.

        metr_step_gs = 0.8;
        metr_step_2P = 2.0;

    }
    
    csvwriter writer_pos(csv_path + "Ex5_1_GS_pos" + gauss_fn + ".csv");
    csvwriter writer_r(csv_path + "Ex5_1_GS_r" + gauss_fn + ".csv");

    // First CL argument: Dump eq. data??
    if (argc > 1 && std::stoi(std::string(argv[1]))) {

        // Dump equilibration data if asked by user

        writer_r.change_file(csv_path + "Ex5_1_GS_eq" + gauss_fn + ".csv");
        eq_x_GS = Equilibrate(rnd, metr_step_gs, H_GS, 15.0, &writer_r, gauss_sample);       // Start 10 <r0> away
        writer_r.change_file(csv_path + "Ex5_1_2P_eq" + gauss_fn + ".csv");
        eq_x_2P = Equilibrate(rnd, metr_step_2P, H_2P, 50.0, &writer_r, gauss_sample);       // 2P orbital extends more in space (see <r0>), start farther away
        writer_r.change_file(csv_path + "Ex5_1_GS_r" + gauss_fn + ".csv");

    } else {

        eq_x_GS = Equilibrate(rnd, metr_step_gs, H_GS, 12.0, nullptr, gauss_sample);
        eq_x_2P = Equilibrate(rnd, metr_step_2P, H_2P, 24.0, nullptr, gauss_sample);

    }

    // Run simulation after equilibration
    std::cout << "Simulating position for GS orbital..." << std::endl;
    Run_Sim(N, L, rnd, metr_step_gs, H_GS, writer_r, &writer_pos, eq_x_GS, gauss_sample);

    writer_r.change_file(csv_path + "Ex5_1_2P_r" + gauss_fn + ".csv");
    writer_pos.change_file(csv_path + "Ex5_1_2P_pos" + gauss_fn + ".csv");
    std::cout << "Simulating position for 2P orbital..." << std::endl;
    Run_Sim(N, L, rnd, metr_step_2P, H_2P, writer_r, &writer_pos, eq_x_2P, gauss_sample);

    return 0;

}