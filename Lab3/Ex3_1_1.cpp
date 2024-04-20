#include "../Libs/RanGen/random.h"
#include "../Libs/CsvWriter/csvwriter.h"

#include <cmath>

const static std::string csv_path = "./CSV/";

double error(double avg_unif, double avg2_unif, size_t n) {

    if(n == 0) {

        return 0;

    } else {

        return sqrt((avg2_unif - avg_unif * avg_unif) / n);

    }

}

double ran_gauss(Random &rnd, double mu, double sigma) {

    double rho = rnd.Rannyu();
    double theta = rnd.Rannyu() * 2 * M_PI;

    double r = sqrt(-2.0 * pow(sigma, 2) * log(1 - rho));
    double x = r * cos(theta) + mu;

    return x;

}

double max(double a, double b) {

    return (a <= b) ? b : a;

}

int main() {

    Random rnd;
    rnd.Init_Random_Gen();

    const uint32_t M = 10000;       // Number of repetitions
    const uint32_t N = 100;         // Number of blocks
    uint32_t L = M / N;

    const double S_0 = 100;
    const double T = 1;
    const double K = 100;
    const double r = 0.1;
    const double sigma = 0.25;

    std::vector<double> C_acc; 
    std::vector<double> C_acc2;

    std::vector<double> P_acc; 
    std::vector<double> P_acc2;

    for(size_t i = 0; i < N; ++i) {

        double C = 0;
        double P = 0;

        for(size_t j = 0; j < L; ++j) {

            double Z = ran_gauss(rnd, 0, 1);
            double S = S_0 * exp((r - pow(sigma, 2) / 2.0) * T + sigma * Z * sqrt(T));
            
            C += exp(-r * T) * max(0, (S - K));
            P += exp(-r * T) * max(0, (K - S));

        }
        
        C_acc.push_back(C / L);
        C_acc2.push_back(pow(C / L, 2));

        P_acc.push_back(P / L);
        P_acc2.push_back(pow(P / L, 2));

    }

    std::vector<double> C_s(N, 0);
    std::vector<double> C_s_2(N, 0);
    std::vector<double> P_s(N, 0);
    std::vector<double> P_s_2(N, 0);

    std::vector<std::vector<double>> data_out(N, std::vector<double>(4, 0));

    for (size_t i = 0; i < N; ++i) {

        for (size_t j = 0; j < i + 1; ++j) {

            C_s[i] += C_acc[j];
            C_s_2[i] += C_acc2[j];

            P_s[i] += P_acc[j];
            P_s_2[i] += P_acc2[j];

        }

        C_s[i] /= (i + 1);
        C_s_2[i] /= (i + 1);

        P_s[i] /= (i + 1);
        P_s_2[i] /= (i + 1);

        data_out[i][0] = C_s[i];
        data_out[i][1] = error(C_s[i], C_s_2[i], i);
        data_out[i][2] = P_s[i];
        data_out[i][3] = error(P_s[i], P_s_2[i], i);

    }

    csvwriter writer(csv_path + "Ex3_1_dir.csv");

    for(auto item : data_out) {

        writer.write_csv_line(item);

    }

    return 0;

}