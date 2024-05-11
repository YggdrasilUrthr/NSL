#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <iostream>

// Header for the chromosome class

class Chromosome {

    private:

        uint32_t m_N_Cities;                                    // Number of cities (alleles) in a chromosome
        std::vector<uint32_t> m_Seq;                            // Chromosome sequence (sequence of cities)
        double m_Fit;                                           // Chromosome fitness value
    
        bool CheckSeq(std::vector<uint32_t> In_Seq);            // Check if the argument sequence is valid

    public:

        std::vector<uint32_t> ReadSeq() const;                  // Return non-modifiable chromosome sequence
        std::vector<uint32_t> * ReadModSeq();                   // Return modifiable chromosome sequence
        bool WriteSeq(std::vector<uint32_t> In_Seq);            // Change chromosome sequence
        void SetFit(double Fit);                                // Assign fitness
        double GetFit() const;                                  // Read fitness
        void PrintData();                                       // Print the chromosome sequence

        Chromosome(uint32_t N_Cities);                          // Create a chromosome with N_Cities with an all-zero sequence
        Chromosome(std::vector<uint32_t> In_Seq);               // Create a chromosome from an existing cities sequence

};