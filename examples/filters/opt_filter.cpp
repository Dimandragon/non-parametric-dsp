import signals;
import integrators;
import derivators;
import inst_freq_computers;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;
import filters;
import phase_computers;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    for (auto i = 0; i < 5000; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 50) * 100000 +
            std::sin(static_cast<double>(i) / 200) +
            std::sin(static_cast<double>(i) / 400) * 1000) ;
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
            phase_computer;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
        NP_DSP::ONE_D::PHASE_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage, decltype(phase_computer)>
            inst_freq_computer (integrator, derivator, phase_computer);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased<NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear> inst_freq_computer2;
    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
        decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average> non_opt_filter;


    NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilter<double, decltype(non_opt_filter), decltype(inst_freq_computer2),
        decltype(phase_computer)> filter(non_opt_filter, inst_freq_computer2, phase_computer);

    filter.compute(signal1, signal2, &signal3);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
}