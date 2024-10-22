#include "icecream.hpp"

#include <cmath>
#include <vector>
#include <matplot/matplot.h>
#include <iostream>

#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <npdsp_config.hpp>

int main(){
    std::vector<double> mydata = {};
    NP_DSP::ONE_D::SimpleVecWrapper wrapper(mydata);
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal(wrapper);
    for (auto i = 0; i < 500; i++) {
        mydata.push_back(static_cast<double>(i));
        //mydata.push_back(std::rand());
        //mydata.push_back(std::sin(static_cast<double>(i) / 50) * 10 +
        //  std::sin(static_cast<double>(i) / 200) +
        //  std::sin(static_cast<double>(i) / 400) * 5) ;

    }
    signal.show(NP_DSP::ONE_D::PlottingKind::Interpolate);
    std::vector<double> mydata2 = {};
    
    int j = -(signal.size()*5);
    for (int i = -(signal.size()*5); i < static_cast<int>(signal.size()*5); i++){
        mydata2.push_back(signal.interpolate(static_cast<double>(i)/2.5, NP_DSP::ONE_D::SignalKind::Universal));
    }
    NP_DSP::ONE_D::SimpleVecWrapper wrapper2(mydata2);
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal2(wrapper2);
    
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, 
        true> signal1;
    //NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, 
    //    true> signal3;
    
    signal1.base->vec->push_back(0.0);
    for (int i = 1; i < 500; i++){
        signal1.base->vec->push_back(static_cast<double>(signal1[i-1] + std::rand() % 4));
    }


    std::vector<double> plotting_vector;
    for (int i = -1000; i < 1500; i++){
        plotting_vector.push_back(signal1.interpolate(i, NP_DSP::ONE_D::SignalKind::Monotone));
    }
    matplot::plot(plotting_vector);
    matplot::show();

    std::cout << "AAAAAAAAA " << signal1.size() << " " << signal1[signal1.size()-1] << std::endl;
    for (int i = -10; i <= signal1[signal1.size()-1] + 10; i++){
        auto temp = signal1.findMonotone(i, {}, {}, {}, 0.01);
        IC(i, temp);
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    for (int i = -10; i <= signal1[signal1.size()-1] + 10; i++){
        auto temp = signal1.findMonotone(i, {}, {}, {signal1.findMonotone(i-1, {}, {}, {}, 0.01)}, 0.01);
        IC(i, temp);
    }

    //signal2.show(NP_DSP::ONE_D::PlottingKind::Interpolate);
    return 0;
}