#include <cmath>
#include <vector>
#include <signals.hpp>
#include "approximators.hpp"
#include <cmath>
#include "icecream.hpp"

int main(){
    int N = 100;
    NP_DSP::ONE_D::APPROX::RBFBasedWithNoTrain approximator;
    std::vector<double> data;

    for (int i = 0; i < N; i++){
        data.push_back(std::sin(((double)i) / 3.0));
    }

    approximator.loadData(data);

    std::vector<double> approximated_data;
    for (int i = 0; i < N * N; i++){
        approximated_data.push_back(approximator.compute<double>(static_cast<double>(i) / N));
    }

    matplot::plot(approximated_data);
    matplot::show();

    N = 10;

    std::vector<std::vector<double>> data2d;
    for (int i = 0; i < N; i++){
        data2d.push_back({});
        for (int j = 0; j < N; j++){
            double r = sqrt((i - N/2.0)*(i - N/2.0) + (j - N/2.0)*(j - N/2.0));
            
            data2d[i].push_back(std::cos(r) / sqrt(r + 1));
            //IC(i, j, r, data2d[i][j]);
        }
    }

    std::vector<std::vector<double>> x_for_load;
    std::vector<std::vector<double>> y_for_load;
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            x_for_load.push_back({});
            x_for_load[i * N + j].push_back(i);
            x_for_load[i * N + j].push_back(j);
            y_for_load.push_back({});
            y_for_load[i * N + j].push_back(data2d[i][j]);
        }
    }

    NP_DSP::ONE_D::APPROX::RBFBasedWithNoTrain approximator2d;
    approximator2d.kind = NP_DSP::ONE_D::APPROX::RBFKind::TPS;
    approximator2d.lambda_v = 1.0;
    approximator2d.loadNDData<decltype(x_for_load)>(x_for_load, y_for_load, 2, 1, N*N);

    using namespace matplot;
    auto [X, Y] = meshgrid(iota(0, 0.1, N));
    auto Z = transform(X, Y, [&](double x, double y) {
        std::vector<double> idx = {x, y};
        std::vector<double> result = {0.0};
        approximator2d.compute
            <decltype(idx), 
            decltype(result)
            >(idx, result);
        return result[0];
    });
    mesh(X, Y, Z);

    show();

    return 0;
}