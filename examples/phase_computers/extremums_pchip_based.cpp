#include <phase_computers.hpp>
#include <signals.hpp>
#include <integrators.hpp>
#include <derivators.hpp>
#include <cmath>
#include <vector>
#include <npdsp_concepts.hpp>
#include <string>
#include <complex>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    SignalT signal3;
    SignalT compute_buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;

    for (auto i = 0; i < 500; i++){
        signal1.base->vec->push_back(std::cos(i * 3.14 / (i / 100 + 1)));
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedUsingPCHIP<double,
            NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple,
                decltype(derivator)> phase_computer;
    phase_computer.derivator = derivator;

    phase_computer.compute(signal1, signal2, &signal3);

    //for (int i = 0; i < signal1.size(); i++){
    //    IC(signal1[i], signal2[i], i);
    //}
    for(int i = 0; i < signal1.size(); i++){
        signal3[i] = phase_computer.derive(i) / (std::numbers::pi * 2.0);
    }

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal3.show(NP_DSP::ONE_D::PlottingKind::Simple);
}