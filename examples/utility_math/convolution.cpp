//
// Created by dmitry on 18.02.24.
//
import utility_math;
import signals;
import npdsp_concepts;
import <vector>;
import <complex>;
import <cstdlib>;

int main(){
    auto signal1 = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});
    using SignalT = decltype(signal1);
    auto signal2 = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});
    auto signal3 = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});

    for (int i = 0; i < 50; i++) {
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal1.base)->vec->push_back(std::rand());
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal3.base)->vec->push_back(0.0);
    }
    for (int i = 0; i < 10; i++) {
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal2.base)->vec->push_back(0.1);
    }

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::fastConvolution(signal1, signal2, signal1);

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}