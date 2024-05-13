#include <mpi.h>
#include <iostream>
#include <sstream>

#include "../Lab9/Population.h"
#include "../Libs/CsvWriter/csvwriter.h"

const static std::string csv_path = "./CSV/";

void Input(uint32_t &, uint32_t &, double &, double &, double &, uint32_t &, uint32_t &);
std::vector<City> LoadCities(std::string);

int main(int argc, char **argv) {
    
    // Initialize MPI and core variables    
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    const static uint32_t N_Cities = 50;                                    // Number of cities (chromosome sequence length)
    uint32_t N_Gens = 300;                                                  // Number of generations
    uint32_t N_Chroms = 1e3;                                                // Number of chromosomes per generation
    double p = 2.0;                                                         // Selection op p parameter, default 2.0
    double Mutate_Prob = 0.1;                                               // Mutation probability, default 10%
    double Crossover_Prob = 0.5;                                            // Crossover probability, default 50%
    uint32_t Migrate = 0;                                                   // Allow migration between continents?
    uint32_t N_Migr = 50;                                                   // Number of generations between migrations

    // Handle input parameters and load cities from external file

    Input(N_Gens, N_Chroms, p, Mutate_Prob, Crossover_Prob, Migrate, N_Migr);
    std::vector<City> Cities = LoadCities("American_capitals.dat");

    Random Rnd;
    Rnd.Init_Random_Gen();

    // Uncomment this if cities has to be randomly distributed inside a square
    /*std::vector<City> Cities;

    for (size_t i = 0; i < N_Cities; ++i) {

        Cities.push_back(City());

        Cities.back().x = Rnd.Rannyu();
        Cities.back().y = Rnd.Rannyu();

    }*/

    // All cores work with the same input => only one core needs to dump it
    if(rank == 0) {

        csvwriter city_writer(csv_path + "Ex10_1_Cities.csv");
        uint32_t idx = 1;

        for (auto city : Cities) {

            city_writer.write_csv_line(std::vector<double>({static_cast<double>(idx), city.x, city.y}));
            ++idx;

        }

    } 

    std::vector<double> BestFits;

    Population Pop(N_Chroms, N_Gens, Cities, rank + 1, false);

    // Let one core be verbose, so as to roughly monitor the execution

    if(rank == 0) {

        Pop.m_verbose = true;

    }

    std::vector<double> Mutate_Probs(MutationsNumber, Mutate_Prob);
    Pop.SetMutateProb(Mutate_Probs);
    Pop.SetCrossoverProb(Crossover_Prob);
    Pop.ChangeP(p);

    if(Migrate != 0) {

        Pop.ChangeGenNumber(N_Migr);                            // Evolve for the number of steps between migrations

        for (size_t i = 0; i < N_Gens / N_Migr; ++i) {

            Pop.Evolve(false);

            uint32_t NoMigrationRank = size + 1;                                            // Core that won't migrate (default n+1 => all cores migrate)      
            uint32_t FirstContinent[size / 2];                                              
            uint32_t SecondContinent[size / 2];

            if(rank == 0) {                                                                 // Only one core (rank = 0) governs the migrations 

                std::cout << "Migration " + std::to_string(i + 1) + "/" + std::to_string(N_Gens / N_Migr) + ": " << std::endl;

                std::vector<uint32_t> NotMigrated(size);
                std::iota(std::begin(NotMigrated), std::end(NotMigrated), 0);               // Cores that have not migrated yet

                // If the number of cores is odd => one randomly selected core does not migrate
                if(size % 2 != 0) {

                    NoMigrationRank = Rnd.GenInt(0, NotMigrated.size() - 1);
                    NotMigrated.erase(NotMigrated.begin() + NoMigrationRank);                                                  

                }

                // Randomly select two cores and exchange their best individuals => 
                // => repeat until one individual from each continent has migrated

                size_t idx = 0;

                do {
 
                    uint32_t m = Rnd.GenInt(0, NotMigrated.size() - 1);
                    FirstContinent[idx] = NotMigrated[m];
                    NotMigrated.erase(NotMigrated.begin() + m);
                    m = Rnd.GenInt(0, NotMigrated.size() - 1);
                    SecondContinent[idx] = NotMigrated[m];
                    NotMigrated.erase(NotMigrated.begin() + m);

                    std::cout << std::to_string(FirstContinent[idx]) + "->" + std::to_string(SecondContinent[idx]) + "\t";
                    ++idx;

                } while(!NotMigrated.empty());

                std::cout << std::endl;

            } 

            // Broadcast to all cores the prescribed migrations
            MPI_Bcast(FirstContinent, size / 2, MPI_INTEGER, 0, MPI_COMM_WORLD);     
            MPI_Bcast(SecondContinent, size / 2, MPI_INTEGER, 0, MPI_COMM_WORLD);           

            // Execute migrations 
            for (size_t i = 0; i < size / 2; ++i) {

                // First continents send then receive
                if(rank == FirstContinent[i]) {

                    uint32_t RecvRnk = SecondContinent[i];
                    MPI_Status Stat;

                    Chromosome * CurrentBest = Pop.GetModBest();
                    MPI_Send(CurrentBest->ReadSeq().data(), CurrentBest->ReadSeq().size(), MPI_INT, RecvRnk, 1, MPI_COMM_WORLD);
                    MPI_Recv(CurrentBest->ReadModSeq()->data(), CurrentBest->ReadModSeq()->size(), MPI_INT, RecvRnk, 1, MPI_COMM_WORLD, &Stat);

                } 

                // Second coninents receive then send
                if (rank == SecondContinent[i]) {

                    uint32_t SendRnk = FirstContinent[i];
                    MPI_Status Stat;

                    Chromosome Tmp = Pop.GetBest();
                    Chromosome * CurrentBest = Pop.GetModBest();
                    MPI_Recv(CurrentBest->ReadModSeq()->data(), CurrentBest->ReadModSeq()->size(), MPI_INT, SendRnk, 1, MPI_COMM_WORLD, &Stat);
                    MPI_Send(Tmp.ReadSeq().data(), Tmp.ReadSeq().size(), MPI_INT, SendRnk, 1, MPI_COMM_WORLD);

                }
            
            }    

        }

    } else {

        // If migration is not needed simply evolve
        Pop.Evolve();

    }

    // Print and save results
    Chromosome Best_Fit = Pop.GetBest();

    csvwriter writer(csv_path + "Ex10_1_Solution_" + std::to_string(rank) + ".csv");
    writer.write_csv_line(Best_Fit.ReadSeq());
    
    writer.change_file(csv_path + "Ex10_1_Hist_" + std::to_string(rank) + ".csv");
    
    std::vector<std::vector<double>> Hist = Pop.GetHist();

    for (size_t i = 0; i < Hist[0].size(); ++i) {

        writer.write_csv_line(std::vector<double>({Hist[0][i], Hist[1][i]}));

    }

    BestFits.push_back(Best_Fit.GetFit());
    
    for (auto item : BestFits) {

        std::cout << "Fit: " << item << std::endl;

    }

    MPI_Finalize();

    return 0;

}

void Input(uint32_t &N_Gens, uint32_t &N_Chroms, double &p, double &Mutate_Prob, double &Crossover_Prob, uint32_t &Migrate, uint32_t &N_Migr) {

    std::ifstream ReadInput;
    ReadInput.open("params.dat");

    if (!ReadInput.is_open()) {

        return;

    }

    ReadInput >> N_Gens;
    ReadInput >> N_Chroms;
    ReadInput >> p;
    ReadInput >> Mutate_Prob;
    ReadInput >> Crossover_Prob;

    ReadInput >> Migrate;
    ReadInput >> N_Migr;

    ReadInput.close();

}

std::vector<City> LoadCities(std::string filename) {

    std::ifstream CitiesInput;
    CitiesInput.open(filename);

    std::vector<City> Cities;

    if(!CitiesInput.is_open()) {

        return Cities;

    }

    std::string ReadLine;
    std::getline(CitiesInput, ReadLine);                // Discard first line (assumed to be a header like in American_capitals.dat)

    while(std::getline(CitiesInput, ReadLine)) {

        std::istringstream LineStream(ReadLine);
        std::string FieldValue;
        std::vector<std::string> ParsedLine;

        while(std::getline(LineStream, FieldValue, ',')) {

            ParsedLine.push_back(FieldValue);

        }
        
        
        Cities.push_back(City());

        Cities.back().x = std::stod(ParsedLine[2]);
        Cities.back().y = std::stod(ParsedLine[3]);

    }

    CitiesInput.close();

    return Cities;

}