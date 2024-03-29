#include <signals.hpp>
#include <integrators.hpp>
#include <derivators.hpp>
#include <inst_freq_computers.hpp>
#include <cmath>
#include <vector>
#include <npdsp_concepts.hpp>
#include <string>
#include <inst_ampl_computers.hpp>
#include <phase_computers.hpp>
#include <filters.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Backward> derivator;

    auto size = 200;

    for (auto i = 0; i < size; i++){
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 50.) * i + i + std::sin(static_cast<double>(i) / 700.) * i);
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
            phase_computer;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)> 
                    inst_freq_computer
                        (integrator, derivator, phase_computer);

    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), decltype(derivator), decltype(inst_freq_computer),
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
            inst_ampl_computer(integrator, derivator, inst_freq_computer);

    static_assert(NP_DSP::ONE_D::is_inst_ampl_computer<decltype(inst_ampl_computer), double>);
    
    NP_DSP::ONE_D::FILTERS::InstAmplNormalizator<double, decltype(integrator), decltype(derivator), 
        decltype(inst_ampl_computer)> inst_ampl_normalizer(integrator, derivator, inst_ampl_computer);

    auto inst_ampl_normalizer2 = NP_DSP::ONE_D::FILTERS::InstAmplNormalizatorUsingInstFreq
        <double, decltype(integrator), decltype(derivator), decltype(inst_ampl_computer),
            false, false, decltype(signal1)>(integrator, derivator, inst_ampl_computer);


    inst_ampl_computer.compute(signal1, signal2, &signal3);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    derivator.compute(signal1, compute_buffer, nullptr);

    auto avg = 0.0;
    auto b = signal1[0];

    for (int i = 0; i < size; i++){
        avg += signal2[i] / size;
    }
    //compute_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for (int i = 0; i < size; i++){
        compute_buffer[i] /= (signal2[i] / avg);
    }
    //compute_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);

    integrator.compute(compute_buffer, signal3, nullptr);

    //signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < size; i++){
        signal3[i] += b;
    }

    //signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);

    inst_freq_computer.compute(signal1, compute_buffer, &signal2);
    inst_ampl_normalizer2.inst_freq = &compute_buffer;
    inst_ampl_normalizer2.compute(signal1, signal2, signal3);

    IC(*(signal2.base->vec));
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    inst_ampl_computer.compute(signal2, signal1, &signal3);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}