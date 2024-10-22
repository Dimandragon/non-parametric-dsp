#include "npdsp_concepts.hpp"
#include <filters.hpp>
#include <signals.hpp>
#include <cmath>

int main(){
    int size = 500;

    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;

    for (int i = 0; i < size; i++){
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0.0);
    }

    NP_DSP::ONE_D::FILTERS::LocalFilter<double, 
        NP_DSP::ONE_D::FILTERS::LocalFilteringType::MakimaInterpolationExtremums> 
            filter;
    filter.compute(signal1, signal2, nullptr);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < size; i++){
        signal1[i] -= signal2[i];
    }
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
}