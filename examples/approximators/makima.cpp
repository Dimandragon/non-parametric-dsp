#include <icecream.hpp>

#include <signals.hpp>
#include <approximators.hpp>
#include <cmath>
#include <vector>

int main(){
    NP_DSP::ONE_D::APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> approximator;
    auto size = 100;

    std::vector<double> y_data(size);
    std::vector<double> x_data(size);
    for (int i = 0; i < size; i++){
        y_data[i] = std::rand();
        x_data[i] = i;
    }

    approximator.loadData(x_data, y_data);

    std::vector<double> plotting_data;

    for (auto i = 0; i < size; i++){
        double add_ = 0.0;
        while (add_ < 1.0){
            auto sample = approximator.compute(i + add_);
            plotting_data.push_back(sample);
            add_ = add_ + 0.1;
        }
    }

    matplot::plot(plotting_data);
    matplot::show();

    return 0;
}