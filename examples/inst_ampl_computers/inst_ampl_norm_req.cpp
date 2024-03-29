#include <signals.hpp>
#include <integrators.hpp>
#include <derivators.hpp>
#include <inst_freq_computers.hpp>
#include <phase_computers.hpp>
#include <inst_ampl_computers.hpp>
#include <cmath>
#include <vector>
#include <npdsp_concepts.hpp>
#include <string>
#include <filters.hpp>

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <ctime>

#include<icecream.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT inst_freq_buffer;
    SignalT inst_freq_mode_buffer;
    SignalT mode;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
            phase_computer;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)> 
                    inst_freq_computer
                        (integrator, derivator, phase_computer);

    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq<double, decltype(integrator), decltype(derivator),
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
            inst_ampl_computer(integrator, derivator, inst_freq_buffer);


    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                non_opt_filter(integrator);

    NP_DSP::ONE_D::FILTERS::InstAmplNormalizatorNaive<double, decltype(inst_ampl_computer)
        , decltype(inst_freq_buffer), decltype(non_opt_filter)> 
        naive_normalizer(inst_ampl_computer, non_opt_filter);

    

    auto size = 200;
    for (auto i = 0; i < size; i++){
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
        inst_freq_buffer.base->vec->push_back(0);
        inst_freq_mode_buffer.base->vec->push_back(0);
        mode.base->vec->push_back(0);
    }

    inst_freq_computer.compute(signal1, inst_freq_buffer, &signal2);
    non_opt_filter.compute(signal1, signal2, &inst_freq_buffer);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for (int i = 0; i < size; i++){
        signal2[i] = signal1[i] - signal2[i];
    }
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    naive_normalizer.inst_freq = &inst_freq_buffer;
    double avg = naive_normalizer.computeGetAvg(signal2, signal3, compute_buffer);
    IC(avg);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for(;;){
        naive_normalizer.computeWithExternalAvg(signal3, signal2, compute_buffer, avg);
        signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
        naive_normalizer.computeWithExternalAvg(signal2, signal3, compute_buffer, avg);
        signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
    }   

    return 0;
}