#include <icecream.hpp>

#include <signals.hpp>
#include <approximators.hpp>
#include <cmath>
#include <vector>

int main(){
    NP_DSP::ONE_D::APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> approximator;
    auto size = 100;
    double d_add = 0.2;

    std::vector<double> y_data(size);
    std::vector<double> x_data(size);
    for (int i = 0; i < size; i++){
        y_data[i] = std::rand();
        x_data[i] = i;
    }
    auto avg = 0.0;
    for (int i = 0; i < y_data.size(); i++){
        avg+= y_data[i];
    }
    for (int i = 0; i < y_data.size(); i++){
        y_data[i] /= avg;
    }

    approximator.loadData(x_data, y_data);

    std::vector<double> plotting_data;

    for (auto i = 0; i < size; i++){
        double add_ = 0.0;
        while (add_ < 1.0){
            auto sample = approximator.compute(i + add_);
            plotting_data.push_back(sample);
            add_ = add_ + d_add;
        }
    }
    matplot::hold(false);
    matplot::plot(plotting_data);
    matplot::hold(true);

    std::vector<double> plotting_data2, plotting_data3;
    for (auto i = 0; i < size; i++){
        double add_ = 0.0;
        while (add_ < 1.0){
            //auto sample = approximator.computeDerive(i + add_, idx / 10.0 + 0.1);
            auto sample = approximator.computeDerive(i + add_);
            plotting_data2.push_back(sample);
            add_ = add_ + d_add;
        }
    }
    matplot::plot(plotting_data2);

    for (auto i = 0; i < size; i++){
        double add_ = 0.0;
        while (add_ < 1.0){
            //auto sample = approximator.computeDerive(i + add_, idx / 10.0 + 0.1);
            auto sample = approximator.computeDerive(i + add_, 1.0);
            plotting_data3.push_back(sample);
            add_ = add_ + d_add;
        }
    }
    matplot::plot(plotting_data3);

    matplot::hold(false);
    matplot::show();

    std::vector<std::vector<double>> plotting_vecs(15);
    for (int idx = 0; idx < 15; idx++){
        for (auto i = 0; i < size; i++){
            double add_ = 0.0;
            while (add_ < 1.0){
                //auto sample = approximator.computeDerive(i + add_, idx / 10.0 + 0.1);
                auto sample = approximator.computeDerive(i + add_, idx / 14.0);
                plotting_vecs[idx].push_back(sample);
                add_ = add_ + d_add;
            }
        }
        matplot::plot(plotting_vecs[idx]);
        matplot::hold(true);
    }
    matplot::hold(false);
    matplot::show();

    return 0;
}