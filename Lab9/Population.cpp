#include "Population.h"

void Population::InitChromosomes() {

    std::vector<uint32_t> Seed(m_N_Cities, 0);
    std::iota(std::begin(Seed), std::end(Seed), 1); 

    m_Chromosomes.push_back(Chromosome(Seed));

    for (size_t i = 0; i < (m_M - 1); i++) {

        uint32_t a, b;

        a = m_Rnd.GenInt(1, m_N_Cities - 1);

        do {

            b = m_Rnd.GenInt(1, m_N_Cities - 1);
        
        } while(a == b);
        
        Permute(Seed, a, b);
        m_Chromosomes.push_back(Chromosome(Seed));

    }
    

}

Population::Population(uint32_t M, uint32_t N_Gen, std::vector<City> Cities, uint32_t PrimeLines) : 
    m_M(M), m_N_Gen(N_Gen), m_Cities(Cities) {

    m_Rnd.Init_Random_Gen(PrimeLines);
    m_N_Cities = m_Cities.size();
    InitChromosomes();

}

Population::Population(uint32_t M, uint32_t N_Gen, std::vector<City> Cities, 
                       std::vector<double> Mutate_Prob, double Crossover_Prob, uint32_t PrimeLines) : 
    m_M(M), m_N_Gen(N_Gen), m_Cities(Cities), m_Mutate_Prob(Mutate_Prob), m_Crossover_Prob(Crossover_Prob) {

    m_Rnd.Init_Random_Gen(PrimeLines);
    m_N_Cities = m_Cities.size();
    InitChromosomes();

}

Population::~Population() {}

void Population::ChangeCities(std::vector<City> NewCities) {

    m_Cities = NewCities;
    m_N_Cities = m_Cities.size();
    InitChromosomes();

}

void Population::ChangeP(double p) {

    m_p = p;

}

void Population::SetMutateProb(std::vector<double> Mutate_Prob) {

    if(Mutate_Prob.size() != MutationsNumber){

        return;

    }

    m_Mutate_Prob = Mutate_Prob;

}

void Population::SetCrossoverProb(double Crossover_Prob) {

    m_Crossover_Prob = Crossover_Prob;

}

void Population::EvaluateFitness() {

    for (auto &Chromosome : m_Chromosomes) {

        double fitness = 0;
        std::vector<uint32_t> seq = Chromosome.ReadSeq();

        for (size_t i = 0; i < m_N_Cities - 1; ++i) {

            fitness += L2_Distance(m_Cities[seq[i] - 1], m_Cities[seq[i + 1] - 1]);

        }

        fitness += L2_Distance(m_Cities[seq.front() - 1], m_Cities[seq.back() - 1]);

        Chromosome.SetFit(fitness);

    }

}

void Population::OrderPop() {

    std::sort(m_Chromosomes.begin(), m_Chromosomes.end(), [](const Chromosome &a, const Chromosome &b){

        return a.GetFit() < b.GetFit();

    });

    m_is_Ordered = true;

}

Chromosome Population::Select() {

    if(!m_is_Ordered) {

        OrderPop();

    }

    double r = m_Rnd.Rannyu();
    uint32_t j = floor(m_M * pow(r, m_p));

    return m_Chromosomes[j];

}

void Population::Mutate(Chromosome &X) {

    uint32_t n, m;

    std::vector<std::function<void(std::vector<uint32_t> &)>> MutateOps;

    do {

        n = m_Rnd.GenInt(1, m_N_Cities - 1);
        m = m_Rnd.GenInt(1, m_N_Cities - 1);    

    } while (n == m);

    MutateOps.push_back(std::bind(Permute, std::placeholders::_1, n, m));
    MutateOps.push_back(std::bind(M_Invert, std::placeholders::_1, m));

    do {

        n = m_Rnd.GenInt(1, m_N_Cities - 2);
        m = m_Rnd.GenInt(1, m_N_Cities - 2);

    } while (n == m);
    
    MutateOps.push_back(std::bind(N_Shift, std::placeholders::_1, n, m));

    do {

        n = m_Rnd.GenInt(2, floor((m_N_Cities - 1.0) / 2.0));

        if(m_N_Cities - n == 0) {

            m = 1;

        } else {

            m = m_Rnd.GenInt(1, m_N_Cities - 2 * n - 1); 

        }   

    } while (n == 0 && m == 0);

    MutateOps.push_back(std::bind(M_Permute, std::placeholders::_1, n, m));

    for(size_t i = 0; i < MutateOps.size(); ++i) {

        double r = m_Rnd.Rannyu();

        if(r < m_Mutate_Prob[i]) {

            MutateOps[i](*X.ReadModSeq());
            m_is_Ordered = false;

        }

    }
    

}

void Population::Crossover(Chromosome &X, Chromosome &Y) {

    double r = m_Rnd.Rannyu();

    if(r >= m_Crossover_Prob) {

        return;

    }

    uint32_t CutPnt = m_Rnd.GenInt(2, m_N_Cities - 1);

    std::vector<uint32_t> * X_Seq = X.ReadModSeq();
    std::vector<uint32_t> * Y_Seq = Y.ReadModSeq();

    std::vector<uint32_t> X_Idxs;
    std::vector<uint32_t> Y_Idxs;

    for(size_t i = CutPnt; i < m_N_Cities; ++i) {

        uint32_t CurrentXIdx = std::find(X_Seq->begin(), X_Seq->end(), Y_Seq->at(i)) - X_Seq->begin();
        uint32_t CurrentYIdx = std::find(Y_Seq->begin(), Y_Seq->end(), X_Seq->at(i)) - Y_Seq->begin();
        
        X_Idxs.push_back(CurrentXIdx);
        Y_Idxs.push_back(CurrentYIdx);

    }

    std::sort(X_Idxs.begin(), X_Idxs.end());
    std::sort(Y_Idxs.begin(), Y_Idxs.end());

    std::vector<uint32_t> X_Temp(m_N_Cities - CutPnt, 0);
    std::vector<uint32_t> Y_Temp(m_N_Cities - CutPnt, 0);

    for(size_t i = 0; i < X_Idxs.size(); ++i) {

        X_Temp[i] = Y_Seq->at(Y_Idxs[i]);
        Y_Temp[i] = X_Seq->at(X_Idxs[i]);

    }
    
    std::swap_ranges(X_Seq->begin() + CutPnt, X_Seq->end(), X_Temp.begin());
    std::swap_ranges(Y_Seq->begin() + CutPnt, Y_Seq->end(), Y_Temp.begin());

}

void Population::Evolve() {

    for (size_t i = 0; i < m_N_Gen; ++i) {

        std::vector<Chromosome> NewGen;
        EvaluateFitness();

        for (size_t j = 0; j < m_Chromosomes.size() / 2; ++j) {

            std::vector<Chromosome> Parents;

            for (size_t k = 0; k < 2; ++k) {

                Parents.push_back(Select());    

            }
            
            Crossover(Parents[0], Parents[1]);

            for (auto item : Parents) {

                Mutate(item);
                NewGen.push_back(item);

            }

        }

        std::swap(m_Chromosomes, NewGen);

    }    

}

Chromosome Population::GetBest() {

    if(!m_is_Ordered) {

        OrderPop();

    }

    return m_Chromosomes[0];

}

void Population::PrintAll() {

    EvaluateFitness();
    OrderPop();

    for (auto item : m_Chromosomes) {
    
        item.PrintData();

    }

    std::cout << "******************************" << std::endl;

}