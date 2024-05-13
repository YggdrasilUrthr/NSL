#include <list>
#include <numeric>
#include <utility>
#include <functional>
#include <iostream>

#include "Functions.h"
#include "Chromosome.h"
#include "../Libs/RanGen/random.h"

// Header for the population class

class Population {
    
    private:
    
        const uint32_t m_M;                                     // Number of chromosomes in population
        uint32_t m_N_Gen;                                       // Number of generations
        Random m_Rnd;                                           // Random generator

        std::vector<Chromosome> m_Chromosomes;                  // Population chromosomes
        std::vector<City> m_Cities;                             // Cities position
        std::vector<double> m_Mutate_Prob;                      // Probability of each mutation
        double m_Crossover_Prob;                                // Crossover probability
        uint32_t m_N_Cities;                                    // Chromosome sequence length

        std::vector<double> m_Best_Fit;                         // Best chromosome fit for each generation
        std::vector<double> m_Best_Half_Avg;                    // Average fit of the best half of population for each generation

        bool m_is_Ordered = false;                              // Are the population chromosomes ordered by fit?
        double m_p = 2;                                         // p parameter in selection operator

        void ClearHistory();                                    // Clears fitness vectors
        void InitChromosomes();                                 // Randomly fill the popultion
        void OrderPop();                                        // Order the population by fitness (lowest fitness first)
        void EvaluateFitness();                                 // Evaluate fitness for all chromosomes
        double EvalBestHalfAvg();                               // Evaluate population best half average
        Chromosome Select();                                    // Selection operator
        void Mutate(Chromosome &X);                             // Mutation operator
        void Crossover(Chromosome &X, Chromosome&Y);            // Crossover operator

    public:
    
        bool m_verbose = true;                                    // Print information during evolution

        // Constructors

        Population(uint32_t M, uint32_t N_Gen, std::vector<City> Cities, uint32_t PrimesLine = 1, bool verbose = true);
        Population(uint32_t M, uint32_t N_Gen, std::vector<City> Cities, 
                   std::vector<double> Mutate_Prob, double Crossover_Prob, uint32_t PrimesLine = 1, bool verbose = true);
        ~Population();

        // Change population properties

        void ChangeP(double p);
        void ChangeGenNumber(uint32_t N_Gen);
        void ChangeCities(std::vector<City> NewCities);
        void SetMutateProb(std::vector<double> Mutate_Prob);
        void SetCrossoverProb(double Crossover_Prob);

        // Evolve and get data

        void Evolve(bool ClearHist = true, bool WriteHist = true);      // Evolution operator
        Chromosome GetBest();                                           // Get the best-fitted chromosome
        Chromosome * GetModBest();                                      // Get a modifiable refernce to bestt-fitted chromosome
        std::vector<std::vector<double>> GetHist();                     // Get popultation history
        void PrintAll();                                                // Print all chromosome data

};