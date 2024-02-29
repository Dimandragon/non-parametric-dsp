import signals;
import derivators;
import <cmath>;
import <vector>;
import npdsp_concepts;
import <string>;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2;

    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Forward> derivator;
    for (auto i = 0; i < 5000; i++){
        signal1.base->vec->push_back(-std::cos(static_cast<double>(i) / 100.0) * 100 + 100);
        signal2.base->vec->push_back(0);
    }
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/derivators/images/signal1.svg");
    NP_DSP::GENERAL::Nil nil;
    derivator.compute(signal1, signal2, &nil);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/derivators/images/signal1_2.svg");
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/derivators/images/signal2.svg");
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator2;
    derivator2.compute(signal1, signal2, &nil);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/derivators/images/signal2_2.svg");
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Backward> derivator3;
    derivator3.compute(signal1, signal2, &nil);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple, "/home/dmitry/projects/non-parametric-dsp/examples/derivators/images/signal2_3.svg");
    return 0;
}