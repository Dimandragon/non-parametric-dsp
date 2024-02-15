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
    NP_DSP::ONE_D::INTEGRATORS::Riman<SignalT, SignalT, NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<SignalT, SignalT, NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    for (auto i = 0; i < 500; i++){
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 5000. * static_cast<double>(i)));
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 10));
        signal2.base->vec->push_back(0);
        //signal3.base->vec->push_back(0.5 / std::numbers::pi / 10.);
        //signal3.base->vec->push_back(0.5 / std::numbers::pi);
        //compute_buffer.base->vec->push_back(0);
    }
    
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
        <SignalT, SignalT, 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>
                inst_freq_computer1;
    

    inst_freq_computer1.compute(signal1, signal2, {});
    
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal1.svg");
    //signal3.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal3.svg");
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal2_1.svg");

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
        <SignalT, SignalT, 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Simple>
                inst_freq_computer2;
    inst_freq_computer2.compute(signal1, signal2, {});
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal2_2.svg");
}