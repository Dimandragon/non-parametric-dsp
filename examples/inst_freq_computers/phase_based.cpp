import signals;
import integrators;
import derivators;
import inst_freq_computers;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<SignalT, SignalT, NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<SignalT, SignalT, NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    for (auto i = 0; i < 500; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 5000. * static_cast<double>(i)));
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) * 2.) * 1000);
        signal2.base->vec->push_back(0);
        //signal3.base->vec->push_back(static_cast<double>(i) / 2400. / std::numbers::pi );
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    
    NP_DSP::ONE_D::InstFreqComputers::PhaseBased
        <SignalT, SignalT, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::InstFreqComputers::InstFreqDerivativeBasedKind::TimeAverage>
                inst_freq_computer(integrator, derivator);
   
    inst_freq_computer.compute(signal1, signal2, compute_buffer);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    

    NP_DSP::ONE_D::InstFreqComputers::PhaseBased<SignalT, SignalT, decltype(integrator), decltype(derivator), 
        NP_DSP::ONE_D::InstFreqComputers::InstFreqDerivativeBasedKind::DeriveAverage> 
        inst_freq_computer2(integrator, derivator);
    inst_freq_computer2.compute(signal1, signal2, compute_buffer);
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal1_1.svg");
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::InstFreqComputers::PhaseBased<SignalT, SignalT, decltype(integrator), decltype(derivator), 
        NP_DSP::ONE_D::InstFreqComputers::InstFreqDerivativeBasedKind::Momental> 
        inst_freq_computer3(integrator, derivator);
    inst_freq_computer3.compute(signal1, signal2, compute_buffer);
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal1_1.svg");
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    
    NP_DSP::ONE_D::InstFreqComputers::ExtremumsBased
        <SignalT, SignalT, 
            NP_DSP::ONE_D::InstFreqComputers::ExtremumsBasedComputeInstFreqKind::Linear>
                inst_freq_computer1;
    

    inst_freq_computer1.compute(signal1, signal2, {});
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
}