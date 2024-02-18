import utility_math;
import signals;
import npdsp_concepts;
import <vector>;
import <complex>;
import <cstdlib>;

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<std::complex<double>>, true> signal1;
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<std::complex<double>>, true> signal2;
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal_plotting;

    for (int i = 0; i < 64; i++) {
        signal1.base->vec->push_back({static_cast<double>(std::rand()), 0.0});
        signal_plotting.base->vec->push_back(signal1[i].real());
        signal2.base->vec->push_back({0.0, 0.0});
    }

    NP_DSP::ONE_D::UTILITY_MATH::fftc2c(*signal1.base->vec, *signal2.base->vec);
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < 64; i++) {
        signal_plotting[i] = signal2[i].real();
    }
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::ifftc2c(*signal2.base->vec, *signal1.base->vec);

    for (int i = 0; i < 64; i++) {
        signal_plotting[i] = signal1[i].real();
    }
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);



    NP_DSP::ONE_D::UTILITY_MATH::fftc2c<decltype(signal1), decltype(signal2), double>(signal1, signal2);
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < 64; i++) {
        signal_plotting[i] = signal2[i].real();
    }
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::ifftc2c<decltype(signal2), decltype(signal1), double>(signal2, signal1);

    for (int i = 0; i < 64; i++) {
        signal_plotting[i] = signal1[i].real();
    }
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);
}