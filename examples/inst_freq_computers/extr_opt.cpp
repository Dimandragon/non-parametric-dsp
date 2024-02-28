#include <icecream.hpp>

import signals;
import integrators;
import derivators;
import inst_freq_computers;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;

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
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal1.base)->vec->push_back(std::sin(static_cast<double>(i) / 5000. * static_cast<double>(i)) * 100000);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal2.base)->vec->push_back(0.0);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal3.base)->vec->push_back(i / 10000. / std::numbers::pi );
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(compute_buffer.base)->vec->push_back(0);
    }
    
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PeriodAndExtremumsBasedExternal
        <double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_freq_computer(integrator, derivator);
    inst_freq_computer.approx_order_coeff = 1.0;
   
    inst_freq_computer.compute(signal1, signal2, &compute_buffer);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    IC(* static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal2.base)->vec);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    inst_freq_computer.approx_order_coeff = 0.5;
    inst_freq_computer.compute(signal1, signal2, &compute_buffer);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    inst_freq_computer.approx_order_coeff = 0.25;
    inst_freq_computer.compute(signal1, signal2, &compute_buffer);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    inst_freq_computer.approx_order_coeff = 0.1;
    inst_freq_computer.compute(signal1, signal2, &compute_buffer);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    inst_freq_computer.approx_order_coeff = 0.05;
    inst_freq_computer.compute(signal1, signal2, &compute_buffer);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    inst_freq_computer.approx_order_coeff = 0.025;
    inst_freq_computer.compute(signal1, signal2, &compute_buffer);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    inst_freq_computer.approx_order_coeff = 0.01;
    inst_freq_computer.compute(signal1, signal2, &compute_buffer);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::DerivativeBased<double, decltype(integrator), decltype(derivator), NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
        inst_freq_computer2(integrator, derivator);
    inst_freq_computer2.compute(signal1, signal2, &compute_buffer);
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/inst_freq_computers/images/signal1_1.svg");
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBased
        <double, NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ExtremumsBasedComputeInstFreqKind::Linear>
            inst_freq_computer1;
    

    inst_freq_computer1.compute(signal1, signal2, nullptr);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
}

//todo