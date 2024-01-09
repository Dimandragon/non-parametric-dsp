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
    NP_DSP::ONE_D::GenericSignal signal(wrapper);
    for (auto i = 0; i < 100; i++) {
        mydata.push_back(static_cast<double>(i));
    }
    for (auto i = 0; i < signal.size(); i++) {
        //std::cout << static_cast<double>(i) << std::endl;
    }
    
    //signal.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/signals/images/signal1.svg");
    //todo

    std::vector<double> mydata2 = {};
    
    int j = -(signal.size()*5);
    //int j = 0;
    //std::cout << j << std::endl;
    //std::cout << signal.size() << std::endl;
    for (int i = -(signal.size()*5); i < static_cast<int>(signal.size()*5); i++){
        //std::cout << "a" << ' ';
        mydata2.push_back(signal.interpolate(static_cast<double>(i)/2.5));
        
        //std::cout << mydata2[mydata2.size() - 1] << ' ';
    }
    NP_DSP::ONE_D::SimpleVecWrapper wrapper2(mydata2);
    NP_DSP::ONE_D::GenericSignal signal2(wrapper2);

    //std::cout << std::endl << signal2.findMonotone<int>(50, {0}, {1000}) << std::endl;
    
    IC(signal2.findMonotone<int>(50, {0}, {1500}));

    //std::string random_str = "AFdserwfaeugnrey4uanguf";
    //IC(random_str);

    //auto answer = signal2.findMonotone<int>(50, {0}, {999});

    //IC(random_str);

    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/signals/images/signal2.svg");
    return 0;
}