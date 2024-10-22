#include <utility_math.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <cstdlib>
#include <filters.hpp>

int main() {
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;

    int size = 500;
    for(int i = 0; i < size; i++){
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0.0);
    }

    NP_DSP::ONE_D::FILTERS::MonoFreqFilters<double, 
        NP_DSP::ONE_D::FILTERS::MonoInstFreqFilteringType::SincPaddedFIR>filter;
    filter.is_low_pass = true;
    filter.freq = 0.2;
    filter.compute(signal1, signal2, signal3);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Spectre);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Spectre);
    return 0;
}