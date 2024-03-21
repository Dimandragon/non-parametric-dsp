#include <icecream.hpp>

#include <vector>
#include <signals.hpp>
#include <derivators.hpp>
#include <cmath>
#include <npdsp_concepts.hpp>
#include <filters.hpp>
#include <integrators.hpp>

int main(){
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;
    using SignalT = decltype(data);
    SignalT out;
    SignalT buffer;
    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central> derivator;

    for (int i = 0; i < 100; i++) {
        data.base->vec->push_back(static_cast<double>(i) * 2.0 + std::cos(i));
        out.base->vec->push_back(0.0);
        buffer.base->vec->push_back(0.159 / 2.);
    }

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
            NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint>,
                NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                filter;

    auto * buffer_ptr = &buffer;
    IC(buffer_ptr);
    filter.compute(data, out, buffer_ptr);
    out.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}