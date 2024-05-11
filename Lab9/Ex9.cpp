#include <iostream>

#include "Population.h"
#include "../Libs/CsvWriter/csvwriter.h"

const static std::string csv_path = "./CSV/";
void Input(uint32_t &, uint32_t &, double &, double &, double &);

int main() {

    const static uint32_t N_Cities = 34;                                    // Number of cities (chromosome sequence length)
    uint32_t N_Gens = 300;                                                  // Number of generations
    uint32_t N_Chroms = 1e3;                                                // Number of chromosomes per generation
    double p = 2.0;                                                         // Selection op p parameter, default 2.0
    double Mutate_Prob = 0.1;                                               // Mutation probability, default 10%
    double Crossover_Prob = 0.5;                                            // Crossover probability, default 50%

    Input(N_Gens, N_Chroms, p, Mutate_Prob, Crossover_Prob);

    std::vector<City> Cities;

    // Generate cities as randomly placed on a circumference

    Random Rnd;
    Rnd.Init_Random_Gen();

    for (size_t i = 0; i < N_Cities; ++i) {

        //double theta = 2 * M_PI / N_Cities * i;                           // Use this if cities must be regularly distributed
        double theta = Rnd.Rannyu() * (2.0 * M_PI);                         // Use this if cities must be randomly distributed
        Cities.push_back(City());

        Cities.back().x = cos(theta);
        Cities.back().y = sin(theta);

    }
    

    Population Pop(N_Chroms, N_Gens, Cities);
    std::vector<double> Mutate_Probs(MutationsNumber, Mutate_Prob);          
    Pop.SetMutateProb(Mutate_Probs);
    Pop.SetCrossoverProb(Crossover_Prob);                                      
    Pop.ChangeP(p);                                                   

    // Evolve population
    std::cout << "Starting TSP with cities randomly placed on a circumference: " << std::endl;
    Pop.Evolve();

    // Write out results

    csvwriter writer(csv_path + "Ex9_1_Hist.csv");
    
    std::vector<std::vector<double>> Hist = Pop.GetHist();

    for (size_t i = 0; i < Hist[0].size(); ++i) {

        writer.write_csv_line(std::vector<double>({Hist[0][i], Hist[1][i]}));

    }

    Chromosome Best_Fit = Pop.GetBest();
    std::cout << "Evolution completed, best sequence and fit obtained: " << std::endl;
    Best_Fit.PrintData();

    writer.change_file(csv_path + "Ex9_1_Cities.csv");
    uint32_t idx = 1;

    for (auto city : Cities) {

        writer.write_csv_line(std::vector<double>({static_cast<double>(idx), city.x, city.y}));
        ++idx;

    }

    writer.change_file(csv_path + "Ex9_1_Solution.csv");
    writer.write_csv_line(Best_Fit.ReadSeq());

    Cities.clear();

    std::cout << "***********************************************************************************" << std::endl;

    // Randomly generate cities inside a square with unitary edge

    for (size_t i = 0; i < N_Cities; ++i) {

        Cities.push_back(City());

        Cities.back().x = Rnd.Rannyu();
        Cities.back().y = Rnd.Rannyu();

    }

    // Evolve population
    std::cout << "Starting TSP with cities randomly placed inside a square: " << std::endl;
    Pop.ChangeCities(Cities);
    Pop.Evolve();

    // Write out results

    writer.change_file(csv_path + "Ex9_2_Hist.csv");
    Hist = Pop.GetHist();

    for (size_t i = 0; i < Hist[0].size(); ++i) {

        writer.write_csv_line(std::vector<double>({Hist[0][i], Hist[1][i]}));

    }

    Best_Fit = Pop.GetBest();
    std::cout << "Evolution completed, best sequence and fit obtained: " << std::endl;
    Best_Fit.PrintData();

    writer.change_file(csv_path + "Ex9_2_Cities.csv");
    idx = 1;

    for (auto city : Cities) {

        writer.write_csv_line(std::vector<double>({static_cast<double>(idx), city.x, city.y}));
        ++idx;

    }

    writer.change_file(csv_path + "Ex9_2_Solution.csv");
    writer.write_csv_line(Best_Fit.ReadSeq());

    return 0;

}

void Input(uint32_t &N_Gens, uint32_t &N_Chroms, double &p, double &Mutate_Prob, double &Crossover_Prob) {

    std::ifstream ReadInput;
    ReadInput.open("params.dat");

    if (!ReadInput.is_open()) {

        std::cout << "No external configuration file found, using default parameters." << std::endl;
        std::cout << "***********************************************************************************" << std::endl;
        return;

    }

    std::cout << "Reading external parameters..." << std::endl;

    ReadInput >> N_Gens;
    ReadInput >> N_Chroms;
    ReadInput >> p;
    ReadInput >> Mutate_Prob;
    ReadInput >> Crossover_Prob;

    ReadInput.close();

    std::cout << "***********************************************************************************" << std::endl;
    

}