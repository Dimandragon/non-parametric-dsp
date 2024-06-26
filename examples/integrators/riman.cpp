#include <signals.hpp>
#include <integrators.hpp>
#include <cmath>
#include <vector>
#include <npdsp_concepts.hpp>
#include <string>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    for (auto i = 0; i < 5000; i++){
        signal1.base->vec->push_back(std::sin(static_cast<double>(i) / 100.0));
        signal2.base->vec->push_back(0);
    }
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/integrators/images/signal1.svg");
    integrator.compute(signal1, signal2, nullptr);
    //signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/integrators/images/signal1_2.svg");
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/integrators/images/signal2.svg");
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByAverage> integrator2;
    integrator2.compute(signal1, signal2, nullptr);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/integrators/images/signal2_2.svg");
    return 0;
}