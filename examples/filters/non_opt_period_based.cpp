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
    auto signal1 = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});
    using SignalT = decltype(signal1);
    auto signal2 = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});
    auto signal3 = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});
    auto compute_buffer = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});
    NP_DSP::ONE_D::INTEGRATORS::Riman<double, NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<double, NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    for (auto i = 0; i < 500; i++){
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal1.base)->vec->push_back(std::sin(static_cast<double>(i) / 5000. * static_cast<double>(i)) * 100000 + std::sin(static_cast<double>(i) / 20000. * static_cast<double>(i)) * 100000);
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) * 2.) * 1000);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal2.base)->vec->push_back(0);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal3.base)->vec->push_back(0);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(compute_buffer.base)->vec->push_back(0);
    }

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
        <double,
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>
                inst_freq_computer;
    
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    inst_freq_computer.compute(signal1, signal2, nullptr);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    
    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, NP_DSP::ONE_D::FILTERS::FilteringType::Median,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                filter (integrator);
    
    filter.compute(signal1, signal3, &signal2);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for (int i = 0; i < signal3.size(); i++) {
        signal3[i] = 0.0;
    }

    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, NP_DSP::ONE_D::FILTERS::FilteringType::Median,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average> 
                filter1 (integrator);

    filter1.compute(signal1, signal3, &signal2);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for (int i = 0; i < signal3.size(); i++) {
        signal3[i] = 0.0;
    }

/*
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::DerivativeBased<SignalT, SignalT, decltype(integrator), decltype(derivator), NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental> 
        inst_freq_computer1(integrator, derivator);
    inst_freq_computer1.compute(signal1, signal2, compute_buffer);
    
    filter.compute(signal1, signal3, nil);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);

    filter1.compute(signal1, signal3, nil);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
*/
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::DerivativeBased<double,
        decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_freq_computer2(integrator, derivator);
    inst_freq_computer2.compute(signal1, signal2, &compute_buffer);
    
    filter.compute(signal1, signal3, &signal2);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for (int i = 0; i < signal3.size(); i++) {
        signal3[i] = 0.0;
    }

    filter1.compute(signal1, signal3, &signal2);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for (int i = 0; i < signal3.size(); i++) {
        signal3[i] = 0.0;
    }

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::DerivativeBased<double,
        decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage> 
                inst_freq_computer3(integrator, derivator);
    inst_freq_computer3.compute(signal1, signal2, &compute_buffer);
    
    filter.compute(signal1, signal3, &signal2);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for (int i = 0; i < signal3.size(); i++) {
        signal3[i] = 0.0;
    }

    filter1.compute(signal1, signal3, &signal2);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
}