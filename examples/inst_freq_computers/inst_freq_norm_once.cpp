#include "derivators.hpp"
#include "inst_freq_computers.hpp"
#include "npdsp_concepts.hpp"
#include "signals.hpp"
#include "cmath"
#include "phase_computers.hpp"
#include "derivators.hpp"
#include "integrators.hpp"
#include <vector>
#include "approximators.hpp"

int main(){
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Backward> derivator;

    using namespace NP_DSP;
    using namespace NP_DSP::ONE_D;

    int size = 5000;
    int res_coeff = 50;

    GenericSignal<SimpleVecWrapper<double>, true> signal;
    GenericSignal<SimpleVecWrapper<double>, true> resampled_signal;
    GenericSignal<SimpleVecWrapper<double>, true> phase;
    GenericSignal<SimpleVecWrapper<double>, true> inst_freq;
    GenericSignal<SimpleVecWrapper<double>, true> compute_buffer;
    std::vector<double> freq_conv;

    APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> approximator;

    std::vector<double> y_approx_picks;
    std::vector<double> x_approx_picks;

    for (int i = 0; i < size / res_coeff + 2; i++){
        y_approx_picks.push_back(std::rand());
        x_approx_picks.push_back(i * res_coeff);
    }
    approximator.loadData(x_approx_picks, y_approx_picks);

    for (int i = 0; i < size; i++){
        signal.base->vec->push_back(approximator.compute(i));
        phase.base->vec->push_back(0.0);
        inst_freq.base->vec->push_back(0.0);
        compute_buffer.base->vec->push_back(0.0);
        resampled_signal.base->vec->push_back(0.0);
    }

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
            phase_computer;
    
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ComputedOnPhase
        <double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_freq_computer(integrator, derivator);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
        <NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>
                inst_freq_computer2;

    while(true){
        phase_computer.compute(signal, phase, &compute_buffer);
        inst_freq_computer.compute(phase, inst_freq, &compute_buffer);
        //inst_freq_computer2.compute(signal, inst_freq, &compute_buffer);

        signal.show(PlottingKind::Simple);
        phase.show(PlottingKind::Simple);
        inst_freq.show(PlottingKind::Simple);

        double base_inst_freq = 1.0/
            INST_FREQ_COMPUTERS::instFreqNormOnce
            (signal, resampled_signal, inst_freq, freq_conv);

        phase_computer.compute(resampled_signal, phase, &compute_buffer);
        inst_freq_computer.compute(phase, inst_freq, &compute_buffer);
        //inst_freq_computer2.compute(resampled_signal, inst_freq, &compute_buffer);

        resampled_signal.show(PlottingKind::Simple);
        phase.show(PlottingKind::Simple);
        inst_freq.show(PlottingKind::Simple);

        base_inst_freq = 1.0/
            INST_FREQ_COMPUTERS::instFreqNormOnce
            (resampled_signal, signal, inst_freq, freq_conv);
    }
    

    return 0;
}