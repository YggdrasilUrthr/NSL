#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <array>
#include <math.h>

const static std::string csv_path = "./CSV/";

// Std. dev. (statistical error)

double error(double avg, double avg2, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt((avg2 - avg * avg) / n);

    }

}

// Is lower <= num < upper?

bool in_interval(double lower, double upper, double num) {

    return (num < upper && num >= lower);

}

int main(){

    csvwriter writer(csv_path + "Ex1_1_1.csv");

    Random rnd;
    rnd.Init_Random_Gen();

    const uint32_t M = 100000;          // Number of throws
    const uint32_t N = 100;             // Number of blocks
    uint32_t L = M / N;                 // Throws per block

    std::vector<std::vector<double>> avg_1_data_vect;
    std::vector<std::vector<double>> avg_2_data_vect;

    for (size_t i = 0; i < N; ++i) {

        double sum_1 = 0;
        double sum_2 = 0;

        // Generate L random numbers and add them to the two block accumulators

        for (size_t j = 0; j < L; ++j) {

            double r = rnd.Rannyu();
            sum_1 += r; 
            sum_2 += pow((r - 0.5), 2);

        }
        
        // Compute block average and averages squared

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

    // Compute progressive averages and errors

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

    // Write out results

    for(auto item : data_out_1) {

        writer.write_csv_line<double>(std::vector<double>(item.begin(), item.end()));

    }    

    writer.change_file(csv_path + "Ex1_1_2.csv");

    
    for(auto item : data_out_2) {

        writer.write_csv_line<double>(std::vector<double>(item.begin(), item.end()));

    }    

    writer.change_file(csv_path + "Ex1_1_3.csv");

    const uint32_t Chi_M = 100;                 // Number of repetitions of chi^2 test
    const uint32_t Chi_Intervals = 100;         // 100 bins for chi^2 histogram
    const uint32_t N_Throws = 10000;
    
    // Compute chi^2 and fill its histogram

    std::vector<double> chi_vect(Chi_M, 0);

    for (size_t i = 0; i < Chi_M; ++i) {
        
        std::vector<double> count_vect(Chi_Intervals, 0);

        // Fill the count vector

        for (size_t j = 0; j < N_Throws; ++j) {

            double r = rnd.Rannyu();

            for (size_t k = 0; k < Chi_Intervals; ++k) {

                count_vect[k] += in_interval(static_cast<double>(k) / Chi_M, static_cast<double>(k + 1) / Chi_M, r);   

            }        

        }

        // Fill chi^2 histogram based on count histogram

        for(auto item : count_vect) {

            chi_vect[i] += pow((item - static_cast<double>(N_Throws) / 100.0), 2) / (static_cast<double>(N_Throws) / 100.0);

        }

    }

    // Write out results

    size_t idx = 0;

    for (auto item : chi_vect) {

        writer.write_csv_line(std::vector<double>({static_cast<double>(idx), item}));
        ++idx;

    }

    return 0;

}