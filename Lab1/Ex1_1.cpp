#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <array>
#include <math.h>

const static std::string rangen_lib_path = "../Libs/RanGen/";
const static std::string csv_path = "./CSV/";

double error(double avg, double avg2, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt((avg2 - avg * avg) / n);

    }

}

void init_random_gen(Random &rnd) {

    int seed[4];
    int p1, p2;
    std::ifstream Primes(rangen_lib_path + "Primes");
    
    if (Primes.is_open()){
        
        Primes >> p1 >> p2 ;
   
    } else std::cerr << "PROBLEM: Unable to open Primes" << std::endl;
   
    Primes.close();

    std::ifstream input(rangen_lib_path + "seed.in");
    std::string property;
    
    if(input.is_open()){
        
        while (!input.eof() ){
            
            input >> property;
            
            if( property == "RANDOMSEED" ){
            
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            
            }

        }
        
        input.close();
   
    } else std::cerr << "PROBLEM: Unable to open seed.in" << std::endl;

    rnd.SaveSeed();

}

bool in_interval(double lower, double upper, double num) {

    return (num < upper && num >= lower);

}

int main(){

    csvwriter writer(csv_path + "Ex1_1_1.csv");

    Random rnd;
    init_random_gen(rnd);

    const uint32_t M = 10000;
    const uint32_t N = 100;
    uint32_t L = M / N;

    std::vector<std::vector<double>> avg_1_data_vect;
    std::vector<std::vector<double>> avg_2_data_vect;

    for (size_t i = 0; i < N; ++i) {

        double sum_1 = 0;
        double sum_2 = 0;

        for (size_t j = 0; j < L; ++j) {

            double r = rnd.Rannyu();
            sum_1 += r; 
            sum_2 += pow((r - 0.5), 2);

        }
        
        double avg_1 = sum_1 / L;
        double avg_1_sqrd = pow(avg_1, 2);

        double avg_2 = sum_2 / L;
        double avg_2_sqrd = pow(avg_2, 2);

        avg_1_data_vect.push_back({avg_1, avg_1_sqrd});
        avg_2_data_vect.push_back({avg_2, avg_2_sqrd});

    }

    std::vector<std::array<double, 2>> sum_prog_vect_1;
    std::vector<std::array<double, 2>> data_out_1;
    std::vector<std::array<double, 2>> sum_prog_vect_2;
    std::vector<std::array<double, 2>> data_out_2;


    for (size_t i = 0; i < N; ++i) {

        sum_prog_vect_1.push_back({0.0, 0.0});
        data_out_1.push_back({0.0, 0.0});

        sum_prog_vect_2.push_back({0.0, 0.0});
        data_out_2.push_back({0.0, 0.0});

        for (size_t j = 0; j < i + 1; ++j) {

            sum_prog_vect_1[i][0] += avg_1_data_vect[j][0];
            sum_prog_vect_1[i][1] += avg_1_data_vect[j][1];
            sum_prog_vect_2[i][0] += avg_2_data_vect[j][0];
            sum_prog_vect_2[i][1] += avg_2_data_vect[j][1];

        }

        sum_prog_vect_1[i][0] /= (i + 1);
        sum_prog_vect_1[i][1] /= (i + 1);
        sum_prog_vect_2[i][0] /= (i + 1);
        sum_prog_vect_2[i][1] /= (i + 1);
        
        data_out_1[i][0] = sum_prog_vect_1[i][0] - 0.5;
        data_out_1[i][1] = error(sum_prog_vect_1[i][0], sum_prog_vect_1[i][1], i);
        data_out_2[i][0] = sum_prog_vect_2[i][0] - 1.0 / 12.0;
        data_out_2[i][1] = error(sum_prog_vect_2[i][0], sum_prog_vect_2[i][1], i);

    }

    for(auto item : data_out_1) {

        writer.write_csv_line<double>(std::vector<double>(item.begin(), item.end()));

    }    

    writer.change_file(csv_path + "Ex1_1_2.csv");

    
    for(auto item : data_out_2) {

        writer.write_csv_line<double>(std::vector<double>(item.begin(), item.end()));

    }    

    writer.change_file(csv_path + "Ex1_1_3.csv");

    const uint32_t Chi_M = 100;
    const uint32_t N_Throws = 10000;
    
    std::vector<double> chi_vect(100, 0);

    for (size_t i = 0; i < 100; ++i) {
        
        std::vector<double> count_vect(Chi_M, 0);

        for (size_t j = 0; j < N_Throws; ++j) {

            double r = rnd.Rannyu();

            for (size_t k = 0; k < Chi_M; ++k) {

                //std::cout << r << std::endl;
                count_vect[k] += in_interval(static_cast<double>(k) / Chi_M, static_cast<double>(k + 1) / Chi_M, r);   

            }        

        }

        for(auto item : count_vect) {

            chi_vect[i] += pow((item - static_cast<double>(N_Throws) / 100.0), 2) / (static_cast<double>(N_Throws) / 100.0);

        }

    }

    size_t idx = 0;

    for (auto item : chi_vect) {

        writer.write_csv_line(std::vector<double>({static_cast<double>(idx), item}));
        ++idx;

    }

    return 0;

}