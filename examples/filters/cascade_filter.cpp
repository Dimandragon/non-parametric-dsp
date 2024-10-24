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

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT inst_freq_buffer;
    SignalT mode;
    SignalT compute_buffer;

    int size = 50;

    for (int i = 0; i < size; i++){
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0.);
        compute_buffer.base->vec->push_back(0.);
        inst_freq_buffer.base->vec->push_back(0.);
        mode.base->vec->push_back(0.);
    }

    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
            phase_computer;
    
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)> 
                    inst_freq_computer
                        (integrator, derivator, phase_computer);

    NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                non_opt_filter(integrator);

    NP_DSP::ONE_D::FILTERS::InstFreqBased<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                non_opt_filter2(integrator);

    NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
        decltype(non_opt_filter2)> cascade_filter(non_opt_filter, non_opt_filter2);

    
    

    for (int i = 0; i < 1000; i++){
        inst_freq_computer.compute(signal1, inst_freq_buffer, &compute_buffer);
        cascade_filter.compute(signal1, signal2, &inst_freq_buffer);

        for (int j = 0; j < size; j++){
            mode[j] = signal1[j] - signal2[j];
        }

        signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
        signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
        mode.show(NP_DSP::ONE_D::PlottingKind::Simple);

        for (int j = 0; j < size; j++){
            signal1[j] = signal2[j];
        }
    }
    cascade_filter.compute(signal1, signal2, &inst_freq_buffer);

    for (int i = 0; i < size; i++){
        mode[i] = signal1[i] - signal2[i];
    }

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    mode.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}