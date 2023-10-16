#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <list>
#include <functional>
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

double ran_exp(Random &rnd, double lambda) {

    double y = rnd.Rannyu();
    double x = -1 / lambda * log(1 - y);

    return x;

}

double ran_cauchy_lorentz(Random &rnd, double Gamma, double mu) {

    double y = rnd.Rannyu();
    double x = Gamma * tan(M_PI * (y - mu - 0.5)); //!!!! VERIFY !!!!!

    return x;

}

std::vector<double> gen_sn(const uint32_t M, const uint32_t N, std::function<double()> ran_func) {

    std::vector<double> out_vect;

    for (size_t i = 0; i < M; ++i) {

        double s_n = 0;

        for (size_t j = 0; j < N; ++j) {

            s_n += ran_func();

        }

        out_vect.push_back(s_n / N);

    }
    
    return out_vect;

}

int main(){

    Random rnd;
    init_random_gen(rnd);

    const uint32_t M = 10000;
    const std::vector<uint32_t> Ns = {1, 2, 10, 100};
    
    std::list<std::vector<double>> s_exp;
    std::list<std::vector<double>> s_cauchy_lorentz;
    
    for (auto N : Ns) {

        s_exp.push_back(gen_sn(M, N, std::bind(ran_exp, rnd, 1)));
        s_cauchy_lorentz.push_back(gen_sn(M, N, std::bind(ran_cauchy_lorentz, rnd, 1, 0)));

    }

    csvwriter writer(csv_path + "Ex1_2_exp.csv");

    for (size_t idx; idx < M; ++idx) {

        std::vector<double> csv_line;

        for(auto item : s_exp) {
            
            csv_line.push_back(item[idx]);

        }
        
        writer.write_csv_line(csv_line);

    }

    writer.change_file(csv_path + "Ex1_2_cl.csv");

    for (size_t idx = 0; idx < M; ++idx) {
        
        std::vector<double> csv_line;

        for(auto item : s_cauchy_lorentz) {
            
            csv_line.push_back(item[idx]);

        }
        
        writer.write_csv_line(csv_line);

    }
    

    return 0;

}