import phase_computers;
import signals;
import integrators;
import derivators;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;
import <complex>;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<SignalT, SignalT, NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<SignalT, SignalT, NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;
    
    for (auto i = 0; i < 50; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 3. /* *static_cast<double>(i)*/) + std::sin(static_cast<double>(i) / 6. /* *static_cast<double>(i)*/));
        //signal1.base->vec->push_back(std::sin(static_cast<double>(i) * 2.) * 1000);
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    auto phase_computer = 
        NP_DSP::ONE_D::PHASE_COMPUTERS::ArctgScaledToExtremums<SignalT, SignalT, SignalT,
            NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(integrator),
                decltype(derivator)>(integrator, derivator);
    phase_computer.tile_size = 25;

    phase_computer.compute(signal1, signal2, signal3);

    //for (int i = 494; i < 500; i++){
        //signal2[i] = 0.0;
    //}

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    //last tile
    
    //for (auto i = 0; i < 500; i++){
        //(*signal1.base->vec)[i] = std::sin(static_cast<double>(i) / 5000. * static_cast<double>(i)) * 100000;
    //}
}