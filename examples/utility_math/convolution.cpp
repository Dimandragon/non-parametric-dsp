//
// Created by dmitry on 18.02.24.
//
#include <utility_math.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <complex>
#include <cstdlib>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;

    for (int i = 0; i < 50; i++) {
        signal1.base->vec->push_back(std::rand());
        signal3.base->vec->push_back(0.0);
    }
    for (int i = 0; i < 10; i++) {
        signal2.base->vec->push_back(0.1);
    }

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::fastConvolution(signal1, signal2, signal1);

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}