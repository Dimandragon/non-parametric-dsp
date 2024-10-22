#include <utility_math.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <cstdlib>
#include <filters.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT signal4;
    SignalT compute_buffer;
    //std::vector<std::complex<double>> resp;
    //std::vector<std::complex<double>> filter;
    std::vector<double> plotting_vector;

    //auto freq_slise = 50;

    NP_DSP::ONE_D::FILTERS::MonoFreqFilters<double,
            NP_DSP::ONE_D::FILTERS::MonoInstFreqFilteringType::Sinc> filter;

    for (int i = 0; i < 500; i++) {
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0.0);
        signal3.base->vec->push_back(0.0);
        compute_buffer.base->vec->push_back(0.0);
        signal4.base->vec->push_back(0.0);
        plotting_vector.push_back(0.0);
        //filter.push_back({0.0, 0.0});
    }

    //IC(signal1.size(), signal2.size(), signal3.size());
    auto ext_info = NP_DSP::ONE_D::UTILITY_MATH::circleExtend(signal1, *(signal2.base->vec), 0.009);
    NP_DSP::ONE_D::UTILITY_MATH::circleExtend(signal1, *(signal3.base->vec), 0.009);
    IC(signal1.size(), signal2.size(), signal3.size());

    filter.computeSincWithMonoFreq(signal2, signal3, ext_info.freq_idx, true);
    IC(signal1.size(), signal2.size(), signal3.size(), ext_info.pad);

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Spectre);

    for (int i = 0; i < signal1.size(); i++){
        compute_buffer[i] = signal3[i + ext_info.pad];
        signal1[i] = signal3[i];
    }
    compute_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);
    compute_buffer.show(NP_DSP::ONE_D::PlottingKind::Spectre);

    NP_DSP::ONE_D::FILTERS::MonoFreqFilters<double, 
        NP_DSP::ONE_D::FILTERS::MonoInstFreqFilteringType::SincPadded> filter2;

    filter2.freq = 0.2;
    filter2.is_low_pass = true;

    filter2.compute(signal1, signal4, signal3);

    signal4.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal4.show(NP_DSP::ONE_D::PlottingKind::Spectre);

    filter2.freq = 0.2;
    filter2.is_low_pass = false;

    filter2.compute(signal1, signal4, signal3);

    signal4.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal4.show(NP_DSP::ONE_D::PlottingKind::Spectre);

    return 0;
}