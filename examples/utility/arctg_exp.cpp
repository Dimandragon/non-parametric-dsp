import <vector>;
import signals;
import derivators;
import <cmath>;
import npdsp_concepts;
import filters;
import integrators;

int main(){
    NP_DSP::GENERAL::Nil nil;
    NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> data;
    using DataT = decltype(data);
    DataT out;
    DataT buffer;

    for (int i = 0; i < 100; i++) {
        data.base->vec->push_back(static_cast<double>(i) * 2.0 + std::cos(i));
        out.base->vec->push_back(0.0);
        buffer.base->vec->push_back(0.159 / 2.);
    }

    NP_DSP::ONE_D::DERIVATORS::FinniteDifference
        <DataT, DataT, NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Central>
            derivator;

    derivator.compute(data, out, nil);

    for (int i = 0; i < data.size(); i++) {
        out[i] = std::atan(out[i]);
    }

    data.show(NP_DSP::ONE_D::PlottingKind::Simple);
    out.show(NP_DSP::ONE_D::PlottingKind::Simple);

    derivator.compute(out, data, nil);
    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for (int i = 0; i < data.size(); i++) {
        data[i] = std::abs(data[i]);
    }
    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter
        <DataT, DataT, DataT, NP_DSP::ONE_D::FILTERS::FilteringType::DerivativeBased,
            NP_DSP::ONE_D::INTEGRATORS::Riman<DataT, DataT, NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint>,
                NP_DSP::ONE_D::FILTERS::InstFreqKind::Average>
                filter;

    filter.compute(data, out, buffer);
    out.show(NP_DSP::ONE_D::PlottingKind::Simple);

    return 0;
}