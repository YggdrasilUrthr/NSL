#include <iostream>

#include "Population.h"
#include "../Libs/CsvWriter/csvwriter.h"

const static std::string csv_path = "./CSV/";

int main() {

    const static uint32_t N_Cities = 34;
    const static uint32_t N_Gens = 100;
    const static uint32_t N_Chroms = 400;

    std::vector<City> Cities;

    for (size_t i = 0; i < N_Cities; ++i) {

        double theta = 2 * M_PI / N_Cities * i;
        Cities.push_back(City());

        Cities.back().x = cos(theta);
        Cities.back().y = sin(theta);

    }
    

    Population Pop(N_Chroms, N_Gens, Cities);
    std::vector<double> Mutate_Prob(MutationsNumber, 0.07);
    Pop.SetMutateProb(Mutate_Prob);
    Pop.SetCrossoverProb(0.5);

    Pop.Evolve();

    Chromosome Best_Fit = Pop.GetBest();
    Best_Fit.PrintData();

    csvwriter writer(csv_path + "Ex9_1_Cities.csv");
    uint32_t idx = 1;

    for (auto city : Cities) {

        writer.write_csv_line(std::vector<double>({static_cast<double>(idx), city.x, city.y}));
        ++idx;

    }

    writer.change_file(csv_path + "Ex9_1_Solution.csv");
    writer.write_csv_line(Best_Fit.ReadSeq());

    Random Rnd;
    Rnd.Init_Random_Gen();

    Cities.clear();

    for (size_t i = 0; i < N_Cities; ++i) {

        Cities.push_back(City());

        Cities.back().x = Rnd.Rannyu();
        Cities.back().y = Rnd.Rannyu();

    }

    Pop.ChangeCities(Cities);
    Pop.Evolve();

    Best_Fit = Pop.GetBest();
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