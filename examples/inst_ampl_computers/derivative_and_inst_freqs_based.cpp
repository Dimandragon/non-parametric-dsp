import signals;
import integrators;
import derivators;
import inst_freq_computers;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;
import inst_ampl_computers;
import phase_computers;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<SignalT, SignalT, NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<SignalT, SignalT, NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    for (auto i = 0; i < 5000; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 500.) * i);
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) * 2.) * 1000);
        signal2.base->vec->push_back(0);
        //signal3.base->vec->push_back(static_cast<double>(i) / 2400. / std::numbers::pi );
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
        <SignalT, SignalT, 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>
                inst_freq_computer;

    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<SignalT, SignalT, 
        SignalT, decltype(integrator), decltype(derivator), decltype(inst_freq_computer), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_ampl_computer(integrator, derivator, inst_freq_computer);

    inst_ampl_computer.compute(signal1, signal2, signal3);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::PHASE_COMPUTERS::ArctgScaledToExtremums<SignalT, SignalT, SignalT,
        decltype(integrator), decltype(derivator)> phase_computer(integrator, derivator);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <SignalT, SignalT, SignalT, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)>
                inst_freq_computer2(integrator, derivator, phase_computer);
    
    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<SignalT, SignalT, 
        SignalT, decltype(integrator), decltype(derivator), decltype(inst_freq_computer2), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_ampl_computer2(integrator, derivator, inst_freq_computer2);

    inst_ampl_computer2.compute(signal1, signal2, signal3);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt<SignalT, SignalT>
            phase_computer3;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <SignalT, SignalT, SignalT, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer3)>
                inst_freq_computer3(integrator, derivator, phase_computer3);

    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<SignalT, SignalT,
        SignalT, decltype(integrator), decltype(derivator), decltype(inst_freq_computer3),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_ampl_computer3(integrator, derivator, inst_freq_computer3);

    inst_ampl_computer3.compute(signal1, signal2, signal3);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::PHASE_COMPUTERS::ArctgScaledToExtremumsSquare<SignalT, SignalT, SignalT,
    decltype(integrator), decltype(derivator),
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
            phase_computer4(integrator, derivator);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <SignalT, SignalT, SignalT, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer4)>
                inst_freq_computer4(integrator, derivator, phase_computer4);

    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<SignalT, SignalT,
        SignalT, decltype(integrator), decltype(derivator), decltype(inst_freq_computer4),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_ampl_computer4(integrator, derivator, inst_freq_computer4);

    inst_ampl_computer2.compute(signal1, signal2, signal3);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
}