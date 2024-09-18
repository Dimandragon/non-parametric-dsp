#include <icecream.hpp>

#include <signals.hpp>
#include <approximators.hpp>
#include <cmath>
#include <vector>

int main(){
    auto size = 2000;

    std::vector<double> y_data(size);
    std::vector<double> x_data(size);
    for (int idx = 0; idx < 1000; idx++){
        NP_DSP::ONE_D::APPROX::PiecewiseCubicHermitePolynomialBasedWithNoTrain<std::vector<double>> approximator;
        
        y_data[0] = 0.0;
        x_data[0] = 0.0;
        for (int i = 1; i < size; i++){
            y_data[i] = x_data[i-1] + 3.14;
            x_data[i] = x_data[i-1] + std::rand() % 100 + 1;
        }

        approximator.loadData(x_data, y_data);

        //std::vector<double> plotting_data;
        bool flag = true;
        for (auto i = 0; i < x_data[size - 1]; i++){
            if (approximator.computeDerive(i) < 0.0){
                std::cout << "incorrect derive " << idx << " " << approximator.computeDerive(i) << std::endl;
                flag = false;
            }
            //auto sample = approximator.compute(i);
            //plotting_data.push_back(sample);
        }

        if (flag){
            std::cout << "complete " << idx << std::endl;
        }

        //matplot::plot(plotting_data);
        //matplot::show();
    }
}