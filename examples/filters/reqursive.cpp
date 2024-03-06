import signals;
import integrators;
import derivators;
import inst_freq_computers;
import phase_computers;
import inst_ampl_computers;
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
    SignalT inst_freq_buffer;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    //NP_DSP::ONE_D::
    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
            phase_computer;

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
            phase_computer_for_mode;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)> inst_freq_computer
                    (integrator, derivator, phase_computer);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer_for_mode)> inst_freq_computer_for_mode
                    (integrator, derivator, phase_computer_for_mode);
    
    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeBasedUsingExternalInstFreq
        <double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_ampl_computer;

    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::ValueBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                non_opt_filter(integrator);

    NP_DSP::ONE_D::FILTERS::RecursiveFilterInstAmplChanges<double, decltype(non_opt_filter),
        decltype(inst_freq_computer), decltype(phase_computer), decltype(inst_ampl_computer),
            decltype(inst_freq_computer_for_mode), decltype(phase_computer_for_mode)> filter
                (non_opt_filter, inst_freq_computer, phase_computer, inst_ampl_computer,
                    inst_freq_computer_for_mode, phase_computer_for_mode);


    for (auto i = 0; i < 50; i++){
        signal1.base->vec->push_back(std::rand());
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }
    inst_freq_computer.variability = 0.5;
    inst_freq_computer.compute(signal1, signal3, &compute_buffer);
    //non_opt_filter.compute(signal1, signal2, &signal3);

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    while (true){
        for(int i = 0; i < 25000; i++){
            non_opt_filter.compute(signal1, signal2, &signal3);
            signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
            //signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
        
            for (auto j = 0; j < 50; j++){
                signal1[j] = signal2[j];
            }
        }
        inst_freq_computer.compute(signal1, signal3, &compute_buffer);
        
        //signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
    }
}