#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <inst_ampl_computers.hpp>
#include <vector>
#include <complex>
#include <cstdlib>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal2;
    //std::vector<std::complex<double>> specter;
    //std::vector<std::complex<double>> buffer;
    //std::vector<double> amplitude;
    //std::vector<double> freqency;

    size_t size = 200;

    for (int i = 0; i < size; i++) {
        //signal1.base->vec->push_back(std::rand());
        signal1.base->vec->push_back(std::sin(i / 2.0) * ((i / 4) % 10));
        signal2.base->vec->push_back(0.0);
        //specter.push_back({0.0, 0.0});
        //buffer.push_back({0.0, 0.0});
        //amplitude.push_back(0.0);
        //freqency.push_back(0.0);
    }

    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    double avg = 0.0;
    for (int i = 0; i < size; i++){
        avg += signal1[i];
    }
    avg = avg / size;
    for (int i = 0; i < size; i++){
        signal1[i] -= avg;
    }
    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::HilbertTransformBased
        <NP_DSP::ONE_D::UTILITY_MATH::HTKind::Mull> inst_ampl_computer;

    inst_ampl_computer.compute(signal1, signal2, nullptr);

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    return 0;
}