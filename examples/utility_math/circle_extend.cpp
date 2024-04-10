#include <utility_math.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <complex>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    std::vector<double> signal_extend;

    for (int i = 0; i < 500; i++){
        signal1.base->vec->push_back(std::rand());
        //signal2.base->vec->push_back(0.0);
    }

    NP_DSP::ONE_D::UTILITY_MATH::ftResamplingData res_data =
            NP_DSP::ONE_D::UTILITY_MATH::getResamplingSize(500, 0.009);

    auto pad = (res_data.new_size - 500) / 2;

    for (int i = 0; i < res_data.new_size; i++){
        signal2.base->vec->push_back(signal1.interpolate
            (static_cast<double>(i - pad), NP_DSP::ONE_D::SignalKind::Harmonic));
    }

    NP_DSP::ONE_D::UTILITY_MATH::circleExtend(signal1, signal_extend, 0.009);

    matplot::plot(signal_extend);
    matplot::hold(true);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    matplot::hold(false);

    return 0;
}