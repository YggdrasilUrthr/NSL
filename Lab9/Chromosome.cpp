#include "Chromosome.h"

Chromosome::Chromosome(uint32_t N_Cities) : m_N_Cities(N_Cities) {

    m_Seq = std::vector<uint32_t>(m_N_Cities, 0);

}

Chromosome::Chromosome(std::vector<uint32_t> In_Seq) {

    m_N_Cities = m_Seq.size();

    if(CheckSeq(In_Seq)) {
        
        m_Seq = In_Seq;

    }

}

std::vector<uint32_t> Chromosome::ReadSeq() const {

    return m_Seq;

}

std::vector<uint32_t> * Chromosome::ReadModSeq() {

    return &m_Seq;

}

bool Chromosome::CheckSeq(std::vector<uint32_t> In_Seq) {

    for (uint32_t i = 1; i <= m_N_Cities; ++i) {

        if (std::count(In_Seq.begin(), In_Seq.end(), i) != 1) {

            return false;

        }

    }
    
    return true;

}

bool Chromosome::WriteSeq(std::vector<uint32_t> In_Seq) {

    bool Is_Okay = CheckSeq(In_Seq);

    if (Is_Okay) {

        m_Seq = In_Seq;

    }

    return Is_Okay;

}

void Chromosome::SetFit(double Fit) {

    m_Fit = Fit;

}

double Chromosome::GetFit() const {

    return m_Fit;

}

void Chromosome::PrintData() {

    for(auto id : m_Seq) {

        std::cout << id << " ";

    }
    
    std::cout << "\t" << m_Fit << std::endl;

}