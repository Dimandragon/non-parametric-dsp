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
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;

    auto size = 5000;

    for (auto i = 0; i < size; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 50.) * i + i + std::sin(static_cast<double>(i) / 700.) * i);
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
            phase_computer;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)> 
                    inst_freq_computer
                        (integrator, derivator, phase_computer);

    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::DerivativeAndInstFreqBased<double, decltype(integrator), decltype(derivator), decltype(inst_freq_computer),
        NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
            inst_ampl_computer(integrator, derivator, inst_freq_computer);

    static_assert(NP_DSP::ONE_D::is_inst_ampl_computer<decltype(inst_ampl_computer), double>);
    
    NP_DSP::ONE_D::INST_AMPL_COMPUTERS::InstAmplNormalizator<double, decltype(integrator), decltype(derivator), 
        decltype(inst_ampl_computer)> inst_ampl_normalizer(integrator, derivator, inst_ampl_computer);

    inst_ampl_computer.compute(signal1, signal2, &signal3);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    //signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    derivator.compute(signal1, compute_buffer, nullptr);

    auto avg = 0.0;
    auto b = signal1[0];

    for (int i = 0; i < size; i++){
        avg += signal2[i] / size;
    }
    //compute_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);
    for (int i = 0; i < size; i++){
        compute_buffer[i] /= (signal2[i] / avg);
    }
    //compute_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);

    integrator.compute(compute_buffer, signal3, nullptr);

    //signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < size; i++){
        signal3[i] += b;
    }

    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);

    inst_ampl_normalizer.compute(signal1, signal2, signal3);

    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}