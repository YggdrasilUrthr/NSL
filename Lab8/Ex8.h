#pragma once

#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <functional>
#include <iostream>

// File paths and save variable
const static std::string csv_path = "./CSV/";
const static std::string file_name_E = csv_path + "E.csv";
const static std::string file_name_SA = csv_path + "SA_E.csv";
const static std::string file_name_SA_params = csv_path + "SA_params.csv";
const static std::string file_name_E_final = csv_path + "E_final.csv";
const static std::string file_name_hist = csv_path + "final_hist.csv";
const static std::string file_name_SA_evol = csv_path + "SA_evol.csv";
bool save;
bool save_SA;
bool save_final;

// Random generator
Random rnd;

// SA parameters

double step_sa;                         // SA step
double T_i;                             // Initial temperature
double T_f;                             // Final temperature
double T_step_lin, T_step_pow;          // Temperature steps for the 2 update laws
double beta_step;                       // Update step for beta
double T, beta;                         // Temperature and beta
uint32_t N_SA;                          // Number of SA blocks
uint32_t L_SA;                          // Number of steps per SA block
uint32_t attempted_SA, accepted_SA;
double mu_0;                            // Starting value for mu
double sigma_0;                         // Starting value for sigma
double mu, sigma;                       // Parameter to optimize
double energy, energy_err;              // Current energy value and error

// Psi sampling parameters

double step;                            // Psi Metropolis algorithm step
double x_0;                             // Starting point for psi sampling
double x;                               // Sampling position
uint32_t N;
uint32_t L;

// Final configuration parameter

uint32_t N_f;
uint32_t L_f;
uint32_t nbins;
std::vector<double> hist;
const double hist_upper = 3.0;
const double hist_lower = -3.0;

// Functions 

void Init();
void Input();
void CalcEnergy(uint32_t N = N, uint32_t L = L, std::string filename = file_name_E);
void Hist();
void Move_SA();
double Error(double, double, double);