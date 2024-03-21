#include <signals.hpp>
#include <integrators.hpp>
#include <derivators.hpp>
#include <inst_freq_computers.hpp>
#include <cmath>
#include <vector>
#include <npdsp_concepts.hpp>
#include <string>
#include <filters.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    for (auto i = 0; i < 500; i++){
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 50) * 10 +
        //    std::sin(static_cast<double>(i) / 200) +
        //    std::sin(static_cast<double>(i) / 400) * 5) ;
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
            phase_computer;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
        NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage, decltype(phase_computer)>
            inst_freq_computer (integrator, derivator, phase_computer);
    
    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
        decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average> non_opt_filter;

    NP_DSP::ONE_D::FILTERS::CascadeFilter<double, decltype(non_opt_filter), 
        decltype(non_opt_filter)> cascade_filter(non_opt_filter, non_opt_filter);

    non_opt_filter.period_muller = 1.6;

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    
    
    for(;;){
        inst_freq_computer.compute(signal1, signal3, &signal2);
        signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
        non_opt_filter.compute(signal1, signal2, &signal3);
        //signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
        for(int i = 0; i < signal1.size(); i++){
            double bf = signal1[i];
            signal1[i] = signal2[i]; 
            signal2[i] = bf - signal1[i];  
        }
        signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
        signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    }
}