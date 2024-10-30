#include "npdsp_concepts.hpp"
#include <filters.hpp>
#include <signals.hpp>
#include <cmath>
#include <approximators.hpp>
#include <cmath>

int main(){
    int size = 100;
    double res_ratio = 1;

    int phase_shift_count = 15;

    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> signal1;
    using SignalT = decltype(signal1);
    SignalT signal2, signal_ref;

    for (int i = 0; i < size; i++){
        signal_ref.base->vec->push_back(std::rand());
    }

    size *= res_ratio;

    for (int i = 0; i < size; i++){
        signal1.base->vec->push_back(signal_ref.interpolate(double(i) /res_ratio, NP_DSP::ONE_D::SignalKind::Universal));
        signal2.base->vec->push_back(0.0);
    }

    NP_DSP::ONE_D::FILTERS::LocalFilter<double, 
        NP_DSP::ONE_D::FILTERS::LocalFilteringType::MakimaInterpolationExtremums> 
            filter;

    filter.extremums_rotation_kind_e = NP_DSP::ONE_D::PHASE_SHIFTERS::RotateKind::NaiveFTFracDir;
    filter.debug = true;
    filter.phase_shifts = {};
    for (int i = 0; i < phase_shift_count; i++){
        filter.phase_shifts.push_back(1.0 / phase_shift_count * i * std::numbers::pi);
    }
    filter.compute(signal1, signal2, nullptr);
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
    signal2.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < size; i++){
        signal1[i] -= signal2[i];
    }
    signal1.show(NP_DSP::ONE_D::PlottingKind::Simple);
}