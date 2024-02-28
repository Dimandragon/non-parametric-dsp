import signals;
import integrators;
import derivators;
import inst_freq_computers;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;
import phase_computers;

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
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal1.base)->vec->push_back(std::sin(static_cast<double>(i) / 5000. * static_cast<double>(i)));
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) * 2.) * 1000);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal2.base)->vec->push_back(0);
        //signal3.base->vec->push_back(static_cast<double>(i) / 2400. / std::numbers::pi );
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal3.base)->vec->push_back(0);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(compute_buffer.base)->vec->push_back(0);
    }
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    NP_DSP::ONE_D::PHASE_COMPUTERS::ArctgScaledToExtremums<double,
        NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(integrator), decltype(derivator)> phase_computer(integrator, derivator);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)>
                inst_freq_computer(integrator, derivator, phase_computer);
   
    inst_freq_computer.compute(signal1, signal2, &compute_buffer);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);


    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::DeriveAverage, decltype(phase_computer)>
            inst_freq_computer2(integrator, derivator, phase_computer);
    inst_freq_computer2.compute(signal1, signal2, &compute_buffer);
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal1_1.svg");
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased<double, decltype(integrator), decltype(derivator),
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::Momental, decltype(phase_computer)>
            inst_freq_computer3(integrator, derivator, phase_computer);
    inst_freq_computer3.compute(signal1, signal2, &compute_buffer);
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal1_1.svg");
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
        <double,
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>
                inst_freq_computer1;
    

    inst_freq_computer1.compute(signal1, signal2, nullptr);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
            phase_computer2;

    NP_DSP::GENERAL::Nil nil;
    
    phase_computer2.compute(signal1, signal3, nullptr);
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ComputedOnPhase
        <double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_freq_computer4(integrator, derivator);

    inst_freq_computer4.compute(signal3, signal2, nullptr);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
}