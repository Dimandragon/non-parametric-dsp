import utility_math;
import signals;
import npdsp_concepts;
import <vector>;
import <complex>;
import <cstdlib>;

int main(){
    auto signal1 = NP_DSP::ONE_D::GenericSignal<std::complex<double>, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<std::complex<double>>>{});
    using SignalT = decltype(signal1);
    auto signal2 = NP_DSP::ONE_D::GenericSignal<std::complex<double>, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<std::complex<double>>>{});
    auto signal_plotting = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});

    for (int i = 0; i < 64; i++) {
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<std::complex<double>> *>(signal1.base)->vec->push_back({static_cast<double>(std::rand()), 0.0});
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(signal_plotting.base)->vec->push_back(signal1[i].real());
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<std::complex<double>> *>(signal2.base)->vec->push_back({0.0, 0.0});
    }

    using WT = NP_DSP::ONE_D::SimpleVecWrapper<std::complex<double>> *;

    NP_DSP::ONE_D::UTILITY_MATH::fftc2c(*static_cast<WT>(signal1.base)->vec, *static_cast<WT>(signal2.base)->vec);
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < 64; i++) {
        signal_plotting[i] = signal2[i].real();
    }
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::ifftc2c(*static_cast<WT>(signal2.base)->vec, *static_cast<WT>(signal1.base)->vec);

    for (int i = 0; i < 64; i++) {
        signal_plotting[i] = signal1[i].real();
    }
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);



    NP_DSP::ONE_D::UTILITY_MATH::fftc2c<double>(signal1, signal2);
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < 64; i++) {
        signal_plotting[i] = signal2[i].real();
    }
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::ifftc2c<double>(signal2, signal1);
    for (int i = 0; i < 64; i++) {
        signal_plotting[i] = signal1[i].real();
    }
    signal_plotting.show(NP_DSP::ONE_D::PlottingKind::Simple);
}