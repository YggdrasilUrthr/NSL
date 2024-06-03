#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <list>
#include <functional>
#include <math.h>

const static std::string csv_path = "./CSV/";

double error(double avg, double avg2, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt((avg2 - avg * avg) / n);

    }

}

// Uniform random generator

double ran_lin(Random &rnd) {

    return rnd.Rannyu();

}

// Exponential random generator

double ran_exp(Random &rnd, double lambda) {

    double y = rnd.Rannyu();
    double x = -1 / lambda * log(1 - y);

    return x;

}

// Cauchy-Lorentx random generator

double ran_cauchy_lorentz(Random &rnd, double Gamma, double mu) {

    double y = rnd.Rannyu();
    double x = Gamma * tan(M_PI * (y - mu - 0.5));

    return x;

}

// This function computes the observable S_n as requested in the exercise for a random generator ran_func

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
    rnd.Init_Random_Gen();

    const uint32_t M = 10000;                                   // Number of throws
    const std::vector<uint32_t> Ns = {1, 2, 10, 100};           // N parameter
    
    std::list<std::vector<double>> s_lin;
    std::list<std::vector<double>> s_exp;
    std::list<std::vector<double>> s_cauchy_lorentz;
    
    for (auto N : Ns) {

        s_lin.push_back(gen_sn(M, N, std::bind(ran_lin, rnd)));
        s_exp.push_back(gen_sn(M, N, std::bind(ran_exp, rnd, 1)));
        s_cauchy_lorentz.push_back(gen_sn(M, N, std::bind(ran_cauchy_lorentz, rnd, 1, 0)));

    }

    // Save data

    csvwriter writer(csv_path + "Ex1_2_lin.csv");

    for (size_t idx; idx < M; ++idx) {

        std::vector<double> csv_line;

        for(auto item : s_lin) {
            
            csv_line.push_back(item[idx]);

        }
        
        writer.write_csv_line(csv_line);

    }

    writer.change_file(csv_path + "Ex1_2_exp.csv");

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