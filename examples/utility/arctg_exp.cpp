#include <icecream.hpp>

import <vector>;
import signals;
import derivators;
import <cmath>;
import npdsp_concepts;
import filters;
import integrators;

int main(){
    NP_DSP::GENERAL::Nil nil;
    auto data = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});
    using DataT = decltype(data);
    auto out = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});
    auto buffer = NP_DSP::ONE_D::GenericSignal<double, true>
        (NP_DSP::GENERAL::Tag<NP_DSP::ONE_D::SimpleVecWrapper<double>>{});

    for (int i = 0; i < 100; i++) {
        //IC(data.base);
        //IC(static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(data.base)->vec);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(data.base)->vec->push_back(static_cast<double>(i) * 2.0 + std::cos(i));
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(out.base)->vec->push_back(0.0);
        static_cast<NP_DSP::ONE_D::SimpleVecWrapper<double> *>(buffer.base)->vec->push_back(0.159 / 2.);
    }

    NP_DSP::ONE_D::DERIVATORS::FinniteDifference
        <double, NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central>
            derivator;

    //IC(data.base, out.base);
    derivator.compute(data, out, nullptr);

    for (int i = 0; i < data.size(); i++) {
        out[i] = std::atan(out[i]);
    }

    data.show(NP_DSP::ONE_D::PlottingKind::Simple);
    out.show(NP_DSP::ONE_D::PlottingKind::Simple);

    derivator.compute(out, data, nullptr);
    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < data.size(); i++) {
        data[i] = std::abs(data[i]);
    }
    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter
        <double, NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
            NP_DSP::ONE_D::INTEGRATORS::Riman<double, NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint>,
                NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                filter;

    auto * buffer_ptr = &buffer;
    IC(buffer_ptr);
    filter.compute(data, out, buffer_ptr);
    out.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}