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

    for (auto i = 0; i < 50; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 3.) + std::sin(static_cast<double>(i) / 6.));
        signal2.base->vec->push_back(0);
        signal3.base->vec->push_back(0);
        compute_buffer.base->vec->push_back(0);
    }

    auto phase_computer = 
        NP_DSP::ONE_D::PHASE_COMPUTERS::ArctgScaledToExtremums<double,
            NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(integrator),
                decltype(derivator)>(integrator, derivator);
    phase_computer.tile_size = 25;

    phase_computer.compute(signal1, signal2, &signal3);

    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);
}