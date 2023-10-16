#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <iostream>

class Chromosome {

    private:

        uint32_t m_N_Cities;
        std::vector<uint32_t> m_Seq;
        double m_Fit;
    
        bool CheckSeq(std::vector<uint32_t> In_Seq);

    public:

        std::vector<uint32_t> ReadSeq() const;
        std::vector<uint32_t> * ReadModSeq();
        bool WriteSeq(std::vector<uint32_t> In_Seq);
        void SetFit(double Fit);
        double GetFit() const;
        void PrintData();

        Chromosome(uint32_t N_Cities);
        Chromosome(std::vector<uint32_t> In_Seq);

};