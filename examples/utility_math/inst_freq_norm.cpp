#include <cmath>
#include <icecream.hpp>

#include <signals.hpp> 
#include <utility_math.hpp>
#include <inst_freq_computers.hpp>
#include <utility_math.hpp>
#include <derivators.hpp>
#include <integrators.hpp>
#include <npdsp_concepts.hpp>
#include <vector>
#include <complex>
#include <filters.hpp>

int main(){
    using DataT = NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true>; 
    DataT data;
    DataT out; 
    DataT inst_freq;
    DataT inst_freq2;
    DataT inst_freq3;
    DataT buffer;
    std::vector<double> freq_conv;
    std::vector<double> freq_conv_image;

    for(int i = 0; i < 500; i++) {
        //data.base->vec->push_back(std::sin(i / 8.0) * 5 + std::sin(i * i));
        data.base->vec->push_back(std::sin(static_cast<double>(i * i) / 5000.));
        out.base->vec->push_back(0.);
        inst_freq.base->vec->push_back(0.);
        buffer.base->vec->push_back(0.);
        freq_conv.push_back(1.);
        inst_freq2.base->vec->push_back(0.0);
        inst_freq3.base->vec->push_back(0.0);
    }
    
    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Backward> derivator;
    
    
    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
            phase_computer;

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::PhaseBased
        <double, decltype(integrator), decltype(derivator), 
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage,
                decltype(phase_computer)> 
                    inst_freq_computer
                        (integrator, derivator, phase_computer);

    inst_freq_computer.compute(data, inst_freq, &buffer);

    inst_freq.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::FILTERS::NonOptPeriodBasedFilter<double, 
        NP_DSP::ONE_D::FILTERS::FilteringType::AverageBased,
            decltype(integrator), NP_DSP::ONE_D::FILTERS::InstFreqKind::Average> non_opt_filter;

    non_opt_filter.period_muller = 1.0;

    for (int i = 0; i < data.size(); i++){
        inst_freq2[i] =  inst_freq[i];
        inst_freq3[i] =  inst_freq[i];
    }
    non_opt_filter.compute(inst_freq3, inst_freq, &inst_freq2);

    inst_freq.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqNorm(data, out, inst_freq, freq_conv, freq_conv_image);

    out.show(NP_DSP::ONE_D::PlottingKind::Simple);

    //inst_freq_computer.compute(out, inst_freq, &buffer);

    //inst_freq.show(NP_DSP::ONE_D::PlottingKind::Simple);

    for(int i = 0; i < data.size(); i++){
        data[i] = out[i];
    }

    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqNorm(data, out, inst_freq, freq_conv, freq_conv_image);

    out.show(NP_DSP::ONE_D::PlottingKind::Simple);

    inst_freq_computer.compute(out, inst_freq, &buffer);

    inst_freq.show(NP_DSP::ONE_D::PlottingKind::Simple);


    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::backInstFreqNorm(out, data, freq_conv);

    IC(data.size());
    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::resampling(data, *(out.base->vec), 1000);

    out.show(NP_DSP::ONE_D::PlottingKind::Simple);

    NP_DSP::ONE_D::UTILITY_MATH::resampling(out, *(data.base->vec), 500);

    data.show(NP_DSP::ONE_D::PlottingKind::Simple);

    //NP_DSP::ONE_D::INST_FREQ_COMPUTERS::backInstFreqNorm(out, data, freq_conv);

    //IC(data.size());
    //data.show(NP_DSP::ONE_D::PlottingKind::Simple);
}