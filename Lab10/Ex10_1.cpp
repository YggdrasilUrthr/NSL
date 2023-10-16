#include <mpi.h>
#include <iostream>

#include "../Lab9/Population.h"
#include "../Libs/CsvWriter/csvwriter.h"

const static std::string csv_path = "./CSV/";

int main(int argc, char **argv) {
    
    int size, rank;
    
    const static uint32_t N_Cities = 34;
    const static uint32_t N_Gens = 100;
    const static uint32_t N_Chroms = 400;

    std::vector<City> Cities;
    Random Rnd;
    Rnd.Init_Random_Gen();

    for (size_t i = 0; i < N_Cities; ++i) {

        Cities.push_back(City());

        Cities.back().x = Rnd.Rannyu();
        Cities.back().y = Rnd.Rannyu();

    }

    csvwriter city_writer(csv_path + "Ex10_1_Cities.csv");
    uint32_t idx = 1;

    for (auto city : Cities) {

        city_writer.write_csv_line(std::vector<double>({static_cast<double>(idx), city.x, city.y}));
        ++idx;

    }

    std::vector<double> BestFits;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Population Pop(N_Chroms, N_Gens, Cities, rank + 1);
    std::vector<double> Mutate_Prob(MutationsNumber, 0.07);
    Pop.SetMutateProb(Mutate_Prob);
    Pop.SetCrossoverProb(0.5);

    Pop.Evolve();

    Chromosome Best_Fit = Pop.GetBest();
    //Best_Fit.PrintData();

    csvwriter writer(csv_path + "Ex10_1_Solution_" + std::to_string(rank) + ".csv");
    writer.write_csv_line(Best_Fit.ReadSeq());
    
    BestFits.push_back(Best_Fit.GetFit());

    MPI_Finalize();
    
    for (auto item : BestFits) {

        std::cout << "Fit: " << item << std::endl;

    }

    return 0;

}