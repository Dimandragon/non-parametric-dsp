#include "approximators.hpp"
#include <cstdint>
#include <vector>
#include <cmath>
#include "npdsp_concepts.hpp"
#include "signals.hpp"
#include "phase_computers.hpp"
#include "derivators.hpp"
#include "integrators.hpp"
#include "inst_freq_computers.hpp"

template<typename T>
void computeExtremums(T & signal, std::vector<int64_t> & extremums){
    extremums.clear();
    extremums.push_back(0);
    for (int i = 1; i < signal.size() - 1; i++) {
        if ((signal[i] >= signal[i - 1] &&
             signal[i] > signal[i + 1]) ||
            (signal[i] > signal[i - 1] &&
             signal[i] >= signal[i + 1]) ||
            (signal[i] <= signal[i - 1] &&
             signal[i] < signal[i + 1]) ||
            (signal[i] < signal[i - 1] &&
             signal[i] <= signal[i + 1])) {
            extremums.push_back(i);
        }
    }
    extremums.push_back(static_cast<int>(signal.size() - 1));
}

int main(){
    using namespace NP_DSP;
    using namespace NP_DSP::ONE_D;

    int size = 5000;
    int res_coeff = 5;

    NP_DSP::ONE_D::INTEGRATORS::Riman<NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator;
    NP_DSP::ONE_D::DERIVATORS::FinniteDifference<NP_DSP::ONE_D::DERIVATORS::FinniteDifferenceType::Backward> derivator;

    NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsBasedNonOpt
        <double, NP_DSP::ONE_D::PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
            phase_computer;
    
    NP_DSP::ONE_D::INST_FREQ_COMPUTERS::ComputedOnPhase
        <double, decltype(integrator), decltype(derivator),
            NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind::TimeAverage>
                inst_freq_computer(integrator, derivator);

    GenericSignal<SimpleVecWrapper<double>, true> signal;
    GenericSignal<SimpleVecWrapper<double>, true> resampled_signal;
    GenericSignal<SimpleVecWrapper<double>, true> phase;
    GenericSignal<SimpleVecWrapper<double>, true> inst_freq;
    GenericSignal<SimpleVecWrapper<double>, true> compute_buffer;

    auto prepare_data = [&](){
        APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> approximator;
        std::vector<double> data_y_approx_picks;
        std::vector<double> data_x_approx_picks;

        for (int i = 0; i < size / res_coeff + 2; i++){
            data_y_approx_picks.push_back(std::rand());
            data_x_approx_picks.push_back(i * res_coeff);
        }
        approximator.loadData(data_x_approx_picks, data_y_approx_picks);

        for (int i = 0; i < size; i++){
            signal.base->vec->push_back(approximator.compute(i));
            resampled_signal.base->vec->push_back(0.0);
            phase.base->vec->push_back(0.0);
            inst_freq.base->vec->push_back(0.0);
            compute_buffer.base->vec->push_back(0.0);
        }
    };
    prepare_data();

    auto showSignal = [&](){
        phase_computer.compute(signal, phase, &compute_buffer);
        inst_freq_computer.compute(phase, inst_freq, &compute_buffer);
        signal.show(PlottingKind::Simple);
        phase.show(PlottingKind::Simple);
        inst_freq.show(PlottingKind::Simple);
    };

    auto showResSignal = [&](){
        phase_computer.compute(resampled_signal, phase, &compute_buffer);
        inst_freq_computer.compute(phase, inst_freq, &compute_buffer);
        resampled_signal.show(PlottingKind::Simple);
        phase.show(PlottingKind::Simple);
        inst_freq.show(PlottingKind::Simple);
    };

/*---------------------------------------------------------------*/

    /*std::vector<int64_t> extremums;
    computeExtremums(signal, extremums);
    double period = static_cast<double>(size - 1) / static_cast<double>(extremums.size() - 1);

    IC(0.5 / period);

    std::vector<double> extremums_old_idx;
    std::vector<double> extremums_new_idx;
    for (int i = 0; i < extremums.size(); i++){
        extremums_old_idx.push_back(extremums[i]);
        extremums_new_idx.push_back(i * period);
    }

    APPROX::PiecewiseCubicHermitePolynomialBasedWithNoTrain<std::vector<double>> idx_approx;
    idx_approx.loadData(extremums_new_idx, extremums_old_idx);

    APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> signal_approx;
    signal_approx.loadData(*(signal.base->vec));

    for (int i = 0; i < size; i++){
        resampled_signal[i] = signal_approx.compute(idx_approx.compute(i));
        //IC(i, idx_approx.compute(i));
    }
    */
    APPROX::PiecewiseCubicHermitePolynomialBasedWithNoTrain<std::vector<double>> idx_approx;
    std::vector<double> extremums;
    INST_FREQ_COMPUTERS::instFreqNormExtrBased(signal, resampled_signal, idx_approx, extremums);
    showSignal();
    showResSignal();

    INST_FREQ_COMPUTERS::backInstFreqNormExtrBased(resampled_signal, signal, idx_approx);
    showSignal();
    /*auto val_expr = [&](double idx){
        return idx_approx.compute(idx);
    };
    auto size_expr = [&](){
        return size;
    };
    ExpressionWrapper<double, double, decltype(val_expr), GENERAL::Nil, decltype(size_expr), false> expr_wrapper(val_expr, size_expr);
    GenericSignalRExpr<decltype(expr_wrapper)> resampled_idx_signal(expr_wrapper);

    APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> signal_res_approx;

    signal_res_approx.loadData(*(resampled_signal.base->vec));
    
    showSignal();
    for(int i = 0; i < size; i++){
        //IC(i, resampled_idx_signal.findMonotone(i, {}, {}, {i}, 0.01));
        signal[i] = signal_res_approx.compute(resampled_idx_signal.findMonotone(i, {}, {}, {i}, 0.01));
    }

    showResSignal();

    showSignal();
    */


    return 0;
}