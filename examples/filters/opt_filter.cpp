import signals;
import integrators;
import derivators;
import inst_freq_computers;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;
import filters;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<SignalT, SignalT, NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<SignalT, SignalT, NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    for (auto i = 0; i < 5000; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 50) * 100000 +
            std::sin(static_cast<double>(i) / 200) +
            std::sin(static_cast<double>(i) / 100) * 1000) ;
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) * 2.) * 1000);
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    NP_DSP::GENERAL::Nil nil;
    
    NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilter<decltype(signal1), decltype(signal2), 
        decltype(signal3), NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased, 
            decltype(integrator), decltype(derivator), NP_DSP::ONE_D::FILTERS::InstFreqComputerKind::extremums_based,
                NP_DSP::ONE_D::FILTERS::PhaseComputingKind::extremums_based_non_opt>
                    filter(integrator, derivator);

    filter.compute(signal1, signal2, signal3);

    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::FILTERS::OptPeriodBasedFilter<decltype(signal1), decltype(signal2),
        decltype(signal3), NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
            decltype(integrator), decltype(derivator), NP_DSP::ONE_D::FILTERS::InstFreqComputerKind::phase_based_time_average,
                NP_DSP::ONE_D::FILTERS::PhaseComputingKind::extremums_based_non_opt>
                    filter2(integrator, derivator);

    filter2.compute(signal1, signal2, signal3);

    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
}