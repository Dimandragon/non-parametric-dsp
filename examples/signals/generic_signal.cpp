#include "icecream.hpp"

import <cmath>;
import <vector>;
import <matplot/matplot.h>;
import <iostream>;

import signals;
import npdsp_concepts;
import npdsp_config;

int main(){
    std::vector<double> mydata = {};
    NP_DSP::ONE_D::SimpleVecWrapper wrapper(mydata);
    NP_DSP::ONE_D::GenericSignal<decltype(wrapper), true> signal(wrapper);
    for (auto i = 0; i < 100; i++) {
        mydata.push_back(static_cast<double>(i));
    }

    std::vector<double> mydata2 = {};
    
    int j = -(signal.size()*5);
    for (int i = -(signal.size()*5); i < static_cast<int>(signal.size()*5); i++){
        mydata2.push_back(signal.interpolate(static_cast<double>(i)/2.5));
    }
    NP_DSP::ONE_D::SimpleVecWrapper wrapper2(mydata2);
    NP_DSP::ONE_D::GenericSignal<decltype(wrapper2), true> signal2(wrapper2);

    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/signals/images/signal2.svg");
    return 0;
}