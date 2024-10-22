#include <iostream>
#include <vector>

#include <matplot/matplot.h>
#include <utility_math.hpp>


int main(){
    auto x1 = 1500.;
    auto x2 = 2000.;
    auto x3 = 2500.;
    auto y1 = 200;
    auto y2 = -400;
    auto y3 = 500;
    
    NP_DSP::ONE_D::UTILITY_MATH::SquarePolynome polynome;
    polynome.solve(x1, x2, x3, y1, y2, y3);

    std::vector<double> data;
    for (int i = 0; i < 5000; i++){
        data.push_back(polynome.compute(i));
    }

    matplot::plot(data);
    matplot::show();
}