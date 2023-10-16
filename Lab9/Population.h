#include <list>
#include <numeric>
#include <utility>
#include <functional>
#include <iostream>

#include "Functions.h"
#include "Chromosome.h"
#include "../Libs/RanGen/random.h"

class Population {
    
    private:
    
        const uint32_t m_M;
        const uint32_t m_N_Gen;
        Random m_Rnd;

        std::vector<Chromosome> m_Chromosomes;
        std::vector<City> m_Cities;
        std::vector<double> m_Mutate_Prob;
        double m_Crossover_Prob;
        uint32_t m_N_Cities;

        bool m_is_Ordered = false;
        double m_p = 2;

        void InitChromosomes();
        void OrderPop();
        void EvaluateFitness();
        Chromosome Select();
        void Mutate(Chromosome &X);
        void Crossover(Chromosome &X, Chromosome&Y);

    public:
    
        Population(uint32_t M, uint32_t N_Gen, std::vector<City> Cities, uint32_t PrimesLine = 1);
        Population(uint32_t M, uint32_t N_Gen, std::vector<City> Cities, 
                   std::vector<double> Mutate_Prob, double Crossover_Prob, uint32_t PrimesLine = 1);
        ~Population();

        void ChangeP(double p);
        void ChangeCities(std::vector<City> NewCities);
        void SetMutateProb(std::vector<double> Mutate_Prob);
        void SetCrossoverProb(double Crossover_Prob);

        void Evolve();
        Chromosome GetBest();
        void PrintAll();

};