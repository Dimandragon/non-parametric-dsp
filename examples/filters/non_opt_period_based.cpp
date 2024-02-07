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
    
    for (auto i = 0; i < 500; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 5000. * static_cast<double>(i)) * 100000 + std::sin(static_cast<double>(i) / 20000. * static_cast<double>(i)) * 100000);
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) * 2.) * 1000);
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    NP_DSP::GENERAL::Nil nil;

    NP_DSP::ONE_D::InstFreqComputers::ExtremumsBased
        <SignalT, SignalT, 
            NP_DSP::ONE_D::InstFreqComputers::ExtremumsBasedComputeInstFreqKind::Linear>
                inst_freq_computer1;
    
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    inst_freq_computer1.compute(signal1, signal2, {});
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < signal2.size(); i++){
        signal2[i] = signal2[i];
    }
    
    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<decltype(signal1), decltype(signal3), 
        decltype(signal2), NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased, 
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average> 
                filter (signal2, integrator);
    
    filter.compute(signal1, signal3, nil);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<decltype(signal1), decltype(signal3), 
        decltype(signal2), NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased, 
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average> 
                filter2 (signal2, integrator);

    filter2.compute(signal1, signal3, nil);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
}