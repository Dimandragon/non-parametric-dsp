#include <utility_math.hpp>
#include <signals.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <cstdlib>
#include <filters.hpp>
#include <derivators.hpp>
#include <integrators.hpp>
#include <inst_freq_computers.hpp>
#include <phase_computers.hpp>

int main() {
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT signal4;

    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)> 
        phase_computer;
    
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ComputedOnPhase
        <double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage>
            inst_freq_computer;


    int size = 500;
    for(int i = 0; i < size; i++){
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0.0);
        signal3.base->vec->push_back(0.0);
        signal4.base->vec->push_back(0.0);
    }

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::FILTERS::SincResLocalFilterWithResReq<double, decltype(phase_computer), decltype(inst_freq_computer)> filter;
    filter.phase_computer = &phase_computer;
    filter.inst_freq_computer = &inst_freq_computer;
    //filter.is_low_pass = true;
    filter.locality_coeff = 5.;
    filter.period_muller = 1.05;
    filter.debug = false;
    

    filter.compute(signal1, signal2, &signal3);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < signal1.size(); i++){
        signal2[i] = signal1[i] - signal2[i];
    }
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    return 0;
}