#pragma once

#include "approximators.hpp"
#include "inst_freq_computers.hpp"
#include "matplot/freestanding/plot.h"
#include <cstddef>
#include <icecream.hpp>

#include <algorithm>
#include <complex>
#include <concepts>
#include <inst_ampl_computers.hpp>
#include <integrators.hpp>
#include <npdsp_concepts.hpp>
#include <signals.hpp>
#include <string>
#include <utility>
#include <utility_math.hpp>
#include <vector>

namespace NP_DSP::ONE_D::FILTERS {
enum class InstFreqKind { Average, Double };

void rotateExtremums(const std::vector<double> & extremums, 
  std::vector<double> & rotated_extremums, double phase_shift)
{
  APPROX::PiecewiseCubicHermitePolynomialBasedWithNoTrain<std::vector<double>> approx;
  rotated_extremums.clear();
  std::vector<double> phase_x;
  std::vector<double> phase_y;
  for (int i = 0; i < extremums.size(); i++){
    phase_x.push_back(extremums[i]);
    phase_y.push_back(i * std::numbers::pi);
  }
  approx.loadData(phase_y, phase_x);
  
  rotated_extremums.push_back(extremums[0]);

  for (int i = 0; i < phase_y.size() - 1; i++){
    double temp = approx.compute(phase_y[i] + phase_shift);
    //rotated_extremums.push_back(approx.compute(phase_y[i] + phase_shift));
    if (rotated_extremums[rotated_extremums.size() - 1] == temp){
      continue;
    }
    rotated_extremums.push_back(temp);
  }
  if (rotated_extremums[rotated_extremums.size() - 1] != 
      extremums[extremums.size() - 1]){
    rotated_extremums.push_back(extremums[extremums.size() - 1]);
  }
}

    enum class FilteringType { DerivativeBased, ValueBased, AverageBased, Median, ValueBasedSmart, DerivativeBasedSmart };

    enum class MaskKind {Gaussian, Triangle, Flat};

    

    template<Signal MaskT>
    int generateConvMask(double inst_freq, double period_muller,
                         MaskT & mask, MaskKind kind, double pow){
        //todo
        double period = 1.0 / inst_freq * period_muller;
        size_t mask_size = period;
        mask.base->vec->clear();
        for (int i = 0; i < mask_size; i++){
            mask.base->vec->push_back(0.0);
        }
        if (kind == MaskKind::Gaussian){
            double width = period * period_muller;
            double c = width / 2.35482;
            double a = 1.0 / (2.506628275 * c);
            //double b = width / 2.0;

            for (int i = 0; i < 2; i++){
                mask.base->vec->push_back(0.0);
            }

            int width_i = static_cast<int>(width) + 2;
            double b = static_cast<double>(width_i) / 2.0 - 0.5;

            for (int i = 0; i < width_i; i++){
                mask[i] = a * std::pow(std::numbers::e, -((i - b) * (i - b)  / (2 * c * c)));
            }
            return width_i;
        }
        else if (kind == MaskKind::Triangle){
            double temp = mask_size / 2.0;
            double mn = UTILITY_MATH::powFact(temp, pow);
            mn = mn * 2.0;
            double mn_temp = mn;
            mn = 1.0 / (mn_temp + (mask_size % 2) * powf(mask_size / 2 + 1, pow));
            if (mask_size % 2 == 1) {
                temp = temp + 1;
            }
            for (int i = 0; i < mask.size(); i++){
                //mask[i].imag() = 0.0;
                mask[i] = 0.0;
            }
            for (int i = 0; i < temp; i++) {
                mask[i] = mn * std::pow((i + 1.0), pow);
                //mask[mask_size - i - 1].real() = mn * std::pow((i + 1.0), pow);
            }
            double sum = 0.0;
            for (int i = 0; i < mask.size(); i++){
                sum = sum + mask[i];
            }
            for (int i = 0; i < mask.size(); i++){
                mask[i] = mask[i]/sum;
            }
        }
        else if (kind == MaskKind::Flat){
            double val = 1.0 / mask.size();
            for (int i = 0; i < mask_size; i++){
                mask[i] = val;
            }
        }
        double avg = 0.0;

        for(int i = 0; i < mask.size(); i++){
            avg += mask[i];
        }

        for (int i = 0; i < mask.size(); i++){
            mask[i] /= avg;
        }

        return mask.size();
    }
    
    template<typename U, FilteringType filtering_type_k,
        Integrator<U> IntegratorT, InstFreqKind inst_freq_k>
    struct NonOptPeriodBasedFilter {
        using MaskT = GenericSignal<SimpleVecWrapper<std::complex<double>>, true>;\
        using AdditionalDataType = SignalPrototype<double>;
        using IdxType = size_t;

        IntegratorT integrator;
        constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
        constexpr static bool is_filter = true;

        double step = 0.05;
        double period_muller = 1.0;

        NonOptPeriodBasedFilter(IntegratorT integrator) {
            this->integrator = integrator;
        }

        NonOptPeriodBasedFilter() {
        }

        //data and inst freq must be not monotone
        template<Signal DataType, Signal OutType, typename InstFreqType>
        void compute(const DataType& data, OutType& out, const InstFreqType * inst_freq) {
            using T = typename OutType::SampleType;
            if constexpr (filtering_type_k == FilteringType::DerivativeBased) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    auto size_expr = [&]() {
                        return data.size();
                    };

                    auto val_expression = [&](size_t idx) {
                        //IC(inst_freq);
                        return (data.interpolate(idx + 0.5 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic) - data.interpolate(
                                    idx - 0.5 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) *
                               inst_freq->interpolate(idx, SignalKind::Stohastic);
                    };

                    ExpressionWrapper<T, size_t,
                        decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                            expr_wrapper(val_expression, size_expr);

                    const GenericSignal<decltype(expr_wrapper), false> expr_signal(expr_wrapper);

                    INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator_new;
                    //todo fix its big crucher
                    integrator_new.compute(expr_signal, out, nullptr);

                    double avg = 0.0;
                    for (int i = 0; i < data.size(); i++){
                        avg += out[i];
                    }
                    avg = avg / data.size();

                    for (int i = 0; i < data.size(); i++){
                        out[i] -= avg;
                    }
                } else if constexpr (inst_freq_kind == InstFreqKind::Double) {
                    auto inst_freq_first_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].first;
                    };
                    auto inst_freq_second_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].second;
                    };
                    auto size_expr = [&]() {
                        return data.size();
                    };
                    ExpressionWrapper<T, size_t,
                        decltype(inst_freq_first_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_first_expr_wrapper(inst_freq_first_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_first_expr_wrapper), false> inst_freq_first(inst_freq_first_expr_wrapper);

                    ExpressionWrapper<T, size_t,
                        decltype(inst_freq_second_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_second_expr_wrapper(inst_freq_second_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_second_expr_wrapper), false> inst_freq_second (inst_freq_second_expr_wrapper);

                    auto val_expression = [&](typename DataType::IdxType idx) {
                        return (data.interpolate(idx + 0.5 / inst_freq_second.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic) -
                                data.interpolate(idx - 0.5 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) /
                               (0.5 / inst_freq_second.interpolate(idx, SignalKind::Stohastic) + 0.5 / inst_freq_first.interpolate(idx, SignalKind::Stohastic));
                    };

                    ExpressionWrapper<T, size_t,
                                decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                            expr_wrapper(val_expression, size_expr);

                    GenericSignal<decltype(expr_wrapper), false> expr_signal(expr_wrapper);

                    INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator_new;

                    //todo fix its big crucher
                    integrator_new.compute(expr_signal, out, nullptr);

                    double avg = 0.0;
                    for (int i = 0; i < data.size(); i++){
                        avg += out[i];
                    }
                    avg = avg / data.size();

                    for (int i = 0; i < data.size(); i++){
                        out[i] -= avg;
                    }
                } else {
                    /*std::unreachable();*/
                }
            } else if constexpr (filtering_type_k == FilteringType::ValueBased) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    auto val_expression = [&](typename DataType::IdxType idx) {
                        return (data.interpolate(idx + 0.25 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)
                                + data.interpolate(idx - 0.25 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) / 2;
                    };
                    for (auto i = 0; i < data.size(); i++) {
                        out[i] = val_expression(i);
                    }
                } else if constexpr (inst_freq_kind == InstFreqKind::Double) {
                    auto inst_freq_first_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].first;
                    };
                    auto inst_freq_second_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].second;
                    };
                    auto size_expr = [&]() {
                        return data.size();
                    };
                    ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_first_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_first_expr_wrapper(inst_freq_first_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_first_expr_wrapper), false> inst_freq_first(
                        inst_freq_first_expr_wrapper);

                    ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_second_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_second_expr_wrapper(inst_freq_second_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_second_expr_wrapper), false> inst_freq_second(
                        inst_freq_second_expr_wrapper);

                    auto val_expression = [&](typename DataType::IdxType idx) {
                        return (data.interpolate(idx + 0.25 / inst_freq_second.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic) +
                                data.interpolate(idx - 0.25 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) / 2.0;
                        /*return UTILITY_MATH::linearInterpolate<double, double>
                        ({idx - 0.25 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), 
                            data.interpolate(idx - 0.25 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)},
                            {idx + 0.25 / inst_freq_second.interpolate(idx, SignalKind::Stohastic),
                            data.interpolate(idx + 0.25 / inst_freq_second.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)},
                            idx);*/
                    };
                    for (auto i = 0; i < data.size(); i++) {
                        out[i] = val_expression(i);
                    }
                    
                } else {
                    /*std::unreachable();*/
                }
            } else if constexpr (filtering_type_k == FilteringType::ValueBasedSmart) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    auto val_expression = [&](typename DataType::IdxType idx) {
                        for (int i = idx; i >= 0; i--){
                            //find 
                            if (i + 0.25 / inst_freq->interpolate(i, SignalKind::Stohastic) <= idx){
                                std::pair<double, double> point1 = {static_cast<double>(i), i + 0.25 / inst_freq->interpolate(i, SignalKind::Stohastic)};
                                std::pair<double, double> point2 = {static_cast<double>(i+1), i + 1. + 0.25 / inst_freq->interpolate(i+1, SignalKind::Stohastic)};
                                auto dy = (point2.second - point1.second);
                                auto dx = (point2.first - point1.first);
                                double idx_new;

                                if (dy == 0) {
                                    idx_new = (point1.first + point2.first) / 2;
                                }
                                idx_new = point1.first + dx * (idx - point1.second) / dy;
                                
                                auto dy1_new = data.interpolate(idx_new, SignalKind::Stohastic) - 
                                    data.interpolate(idx_new - 0.25 / inst_freq->interpolate(idx_new, SignalKind::Stohastic), SignalKind::Stohastic);
                                auto dx1_new = 0.25 / inst_freq->interpolate(idx_new, SignalKind::Stohastic);

                                auto dx2_new = 0.25 / inst_freq->interpolate(idx_new, SignalKind::Stohastic);
                                auto dy2_new = dy1_new;
                                return data.interpolate(idx_new, SignalKind::Stohastic) + dy2_new / dx1_new * dx2_new ;
                            }
                        }
                        
                        
                        return (data.interpolate(idx + 0.25 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)
                                + data.interpolate(idx - 0.25 / inst_freq->interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) / 2;
                    };
                    for (auto i = 0; i < data.size(); i++) {
                        out[i] = val_expression(i);
                    }
                } else if constexpr (inst_freq_kind == InstFreqKind::Double) {
                    auto inst_freq_first_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].first;
                    };
                    auto inst_freq_second_val_expression = [&](typename DataType::IdxType idx) {
                        return (*inst_freq)[idx].second;
                    };
                    auto size_expr = [&]() {
                        return data.size();
                    };
                    ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_first_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_first_expr_wrapper(inst_freq_first_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_first_expr_wrapper), false> inst_freq_first(
                        inst_freq_first_expr_wrapper);

                    ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType,
                                decltype(inst_freq_second_val_expression), GENERAL::Nil, decltype(size_expr), false>
                            inst_freq_second_expr_wrapper(inst_freq_second_val_expression, size_expr);
                    GenericSignal<decltype(inst_freq_second_expr_wrapper), false> inst_freq_second(
                        inst_freq_second_expr_wrapper);

                    auto val_expression = [&](typename DataType::IdxType idx) {
                        for (int i = idx; i >= 0; i--){
                            //find 
                            if (i + 0.25 / inst_freq_second.interpolate(i, SignalKind::Stohastic) <= idx){
                                std::pair<double, double> point1 = {static_cast<double>(i), i + 0.25 / inst_freq_second.interpolate(i, SignalKind::Stohastic)};
                                std::pair<double, double> point2 = {static_cast<double>(i+1), i + 1. + 0.25 / inst_freq_second.interpolate(i+1, SignalKind::Stohastic)};
                                auto dy = (point2.second - point1.second);
                                auto dx = (point2.first - point1.first);
                                double idx_new;

                                if (dy == 0) {
                                    idx_new = (point1.first + point2.first) / 2;
                                }
                                idx_new = point1.first + dx * (idx - point1.second) / dy;
                                
                                auto dy1_new = data.interpolate(idx_new, SignalKind::Stohastic) - 
                                    data.interpolate(idx_new - 0.25 / inst_freq_first.interpolate(idx_new, SignalKind::Stohastic), SignalKind::Stohastic);
                                auto dx1_new = 0.25 / inst_freq_first.interpolate(idx_new, SignalKind::Stohastic);

                                auto dx2_new = 0.25 / inst_freq_second.interpolate(idx_new, SignalKind::Stohastic);
                                auto dy2_new = dy1_new;
                                return data.interpolate(idx_new, SignalKind::Stohastic) + dy2_new / dx1_new * dx2_new ;
                            }
                        }

                        return (data.interpolate(idx + 0.25 / inst_freq_second.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic) +
                                data.interpolate(idx - 0.25 / inst_freq_first.interpolate(idx, SignalKind::Stohastic), SignalKind::Stohastic)) / 2.0;
                    };
                    for (auto i = 0; i < data.size(); i++) {
                        out[i] = val_expression(i);
                    }
                } else {
                    /*std::unreachable();*/
                }
            }
            else if constexpr (filtering_type_k == FilteringType::AverageBased) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = 0;
                        int iters_count = 0;
                        for (double j = -0.5 * period_muller / (*inst_freq)[i]; j <= 0.5 * period_muller / (*inst_freq)[i]; j = j + step * period_muller / (*inst_freq)[i]) {
                            iters_count++;
                            out[i] += data.interpolate(i + j, SignalKind::Stohastic);
                        }
                        out[i] = out[i] / iters_count;
                    }
                } else if (inst_freq_kind == InstFreqKind::Double) {
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = 0;
                        int iters_count = 0;
                        for (double j = -0.5 * period_muller / (*inst_freq)[i].first; j <= 0.5 * period_muller / (*inst_freq)[i].second; 
                            j = j + (step / (*inst_freq)[i].first + step / (*inst_freq)[i].second) * period_muller / 2.0) {
                            iters_count++;
                            out[i] += data.interpolate(i + j, SignalKind::Stohastic);// * step;
                        }
                        out[i] = out[i] / iters_count;
                    }
                }
            } else if constexpr (filtering_type_k == FilteringType::Median) {
                if constexpr (inst_freq_kind == InstFreqKind::Average) {
                    std::vector<T> buffer;
                    for (int i = 0; i < data.size(); i++) {
                        buffer.clear();
                        for (double j = -0.5 / (*inst_freq)[i]; j <= 0.5 / (*inst_freq)[i]; j = j + step / (*inst_freq)[i]) {
                            buffer.push_back(data.interpolate(i + j, SignalKind::Stohastic));
                        }
                        std::sort(buffer.begin(), buffer.end());
                        if (buffer.size() % 2 == 0) {
                            out[i] = buffer[buffer.size() / 2];
                        } else {
                            out[i] = (buffer[buffer.size() / 2] + buffer[buffer.size() / 2 + 1]) / 2.0;
                        }
                    }
                } else if (inst_freq_kind == InstFreqKind::Double) {
                    std::vector<T> buffer;
                    for (int i = 0; i < data.size(); i++) {
                        buffer.clear();
                        //auto step_left = step*(1./inst_freq[i].left + 1./inst_freq[i].right)/1./inst_freq[i].left;
                        for (double j = -0.5 / (*inst_freq)[i].first; j < 0; j = j + step / (*inst_freq)[i].first) {
                            buffer.push_back(data.interpolate(i + j, SignalKind::Stohastic));
                        }
                        for (double j = 0.0; j < 0.5 / (*inst_freq)[i].second; j = j + step / (*inst_freq)[i].second) {
                            buffer.push_back(data.interpolate(i + j, SignalKind::Stohastic));
                        }
                        std::sort(buffer.begin(), buffer.end());
                        if (buffer.size() % 2 == 0) {
                            out[i] = buffer[buffer.size() / 2];
                        } else {
                            out[i] = (buffer[buffer.size() / 2] + buffer[buffer.size() / 2 + 1]) / 2.0;
                        }
                        //todo left and right steps
                    }
                }
            } else {
                /*std::unreachable();*/
            }
        }
    };

    enum class NonLocalFilteringType { Conv, Sinc, SincLocal };
    template<typename U, NonLocalFilteringType kind_e>
    struct NonLocalNonOptFiltering {
        constexpr static bool is_filter = true;
        size_t conv_filter_len = 0;
        double width_muller = 5.0;

        template<Signal DataType, Signal OutType, Signal BufferType>
        void compute(const DataType& data, OutType& out, const BufferType & buffer) {
            if constexpr (kind_e == NonLocalFilteringType::Conv) {
                //buffer is conv mask
                UTILITY_MATH::fastConvolution(data, buffer, out);
                std::vector<typename OutType::SampleType> copy_buffer;

                for (int i = 0; i < conv_filter_len / 2; i++){
                    copy_buffer.push_back(out[i]);
                }
                for (int i = 0; i < data.size() - conv_filter_len / 2; i++){
                    out[i] = out[i + conv_filter_len / 2];
                }
                for (int i = data.size() - conv_filter_len / 2; i < data.size(); i++){
                    out[i] = out[i + conv_filter_len / 2 - data.size()];
                }
                //todo test
            }
            else if constexpr (kind_e == NonLocalFilteringType::Sinc) {
                std::vector<std::complex<double>> buffer1;
                std::vector<std::complex<double>> buffer2;

                for (int i = 0; i < data.size(); i++) {
                    buffer1.push_back({data[i], 0.0});
                    buffer2.push_back({0.0, 0.0});
                }

                UTILITY_MATH::fftc2c(buffer1, buffer2);

                for (int i = 0; i < UTILITY_MATH::getUniqueSpecterSamplesCount(data.size()); i++) {
                    buffer2[i] = buffer2[i] * buffer[i];
                    if (i != 0 && data.size() - i > i) {
                        buffer2[data.size() - i] = buffer2[data.size() - i] * buffer[i];
                    }
                }

                UTILITY_MATH::ifftc2c(buffer2, buffer1);

                for (int i = 0; i < data.size(); i++) {
                    out[i] = buffer1[i].real();
                }
            }
            else if constexpr (kind_e == NonLocalFilteringType::SincLocal){
                /*std::unreachable();*/
            }
        }

        template<Signal DataType, Signal OutType>
        void computeSincWithMonoFreq(const DataType& data, OutType& out, int freq_idx, bool low_pass) {
            if constexpr (kind_e == NonLocalFilteringType::Sinc){
                std::vector<std::complex<double>> buffer1;
                std::vector<std::complex<double>> buffer2;

                for (int i = 0; i < data.size(); i++) {
                    buffer1.push_back({data[i], 0.0});
                    buffer2.push_back({0.0, 0.0});
                }

                UTILITY_MATH::fftc2c(buffer1, buffer2);

                if (low_pass){
                    for(int i = freq_idx; i < UTILITY_MATH::getUniqueSpecterSamplesCount(data.size()); i++){
                        buffer2[i] = buffer2[i] * 0.0;
                        if (i != 0 && data.size() - i > i) {
                            buffer2[data.size() - i] = buffer2[data.size() - i] * 0.0;
                        }
                    }
                }
                else{
                    for(int i = 0; i < freq_idx && i < UTILITY_MATH::getUniqueSpecterSamplesCount(data.size()); i++){
                        buffer2[i] = buffer2[i] * 0.0;
                        if (i != 0 && data.size() - i > i) {
                            buffer2[data.size() - i] = buffer2[data.size() - i] * 0.0;
                        }
                    }
                }

                UTILITY_MATH::ifftc2c(buffer2, buffer1);

                for (int i = 0; i < data.size(); i++) {
                    out[i] = buffer1[i].real();
                }
            }
            else if constexpr (kind_e == NonLocalFilteringType::SincLocal){
                std::vector<std::complex<double>> buffer1;
                std::vector<std::complex<double>> buffer2;
                std::vector<std::complex<double>> filter_non_local_afr;
                std::vector<std::complex<double>> filter;

                for (int i = 0; i < data.size(); i++) {
                    buffer1.push_back({data[i], 0.0});
                    buffer2.push_back({0.0, 0.0});
                    filter_non_local_afr.push_back({1.0, 1.0});
                    filter.push_back({0.0, 0.0});
                }

                if (low_pass){
                    for(int i = freq_idx; i <= UTILITY_MATH::getUniqueSpecterSamplesCount(data.size()); i++){
                        filter_non_local_afr[i] = filter_non_local_afr[i] * 0.0;
                        if (i != 0 && data.size() - i > i) {
                            filter_non_local_afr[data.size() - i] = filter_non_local_afr[data.size() - i] * 0.0;
                        }
                    }
                }
                else{
                    for(int i = 0; i < freq_idx && i <= UTILITY_MATH::getUniqueSpecterSamplesCount(data.size()); i++){
                        filter_non_local_afr[i] = filter_non_local_afr[i] * 0.0;
                        if (i != 0 && data.size() - i > i) {
                            filter_non_local_afr[data.size() - i] = filter_non_local_afr[data.size() - i] * 0.0;
                        }
                    }
                }
                /*std::vector<double> plotting_vector;
                for(int i = 0; i < data.size(); i++){
                    plotting_vector.push_back(filter_non_local_afr[i].real());
                }
                matplot::plot(plotting_vector);
                matplot::hold(true);
                for(int i = 0; i < data.size(); i++){
                    plotting_vector[i] = (filter_non_local_afr[i].imag());
                }
                matplot::plot(plotting_vector);
                matplot::hold(false);
                matplot::show();*/

                UTILITY_MATH::ifftc2c(filter_non_local_afr, filter);

                /*for(int i = 0; i < data.size(); i++){
                    plotting_vector[i] = (filter[i].real());
                }
                matplot::plot(plotting_vector);
                matplot::hold(true);
                for(int i = 0; i < data.size(); i++){
                    plotting_vector[i] = (filter[i].imag());
                }
                matplot::plot(plotting_vector);
                matplot::hold(false);
                matplot::show();*/

                double freq = UTILITY_MATH::getFreqByIdx(data.size(), freq_idx);
                double width = 1.0 / freq * width_muller;
                double c = width / 2.35482;
                double a = 1.0 / (2.506628275 * c);
                double b = 0;

                for(int i = 0; i <= UTILITY_MATH::getUniqueSpecterSamplesCount(data.size()); i++){
                    filter[i] = filter[i] * a * std::pow(std::numbers::e, -((i - b) * (i - b)  / (2 * c * c)));
                    if (i != 0 && data.size() - i > i) {
                        filter[data.size() - i] = filter[data.size() - i] * a * std::pow(std::numbers::e, -((i - b) * (i - b)  / (2 * c * c)));
                    }
                }

                double sum_re = 0.0;
                double sum_im = 0.0;
                for (int i = 0; i < data.size(); i++){
                    sum_re += filter[i].real();
                    sum_im += filter[i].imag();
                }

                for (int i = 0; i < data.size(); i++){
                    auto re= filter[i].real() / sum_re * data.size();
                    auto im= filter[i].imag() / sum_im * data.size();
                    filter[i] = {re, im};
                }

                /*for(int i = 0; i < data.size(); i++){
                    plotting_vector[i] = (filter[i].real());
                }
                matplot::plot(plotting_vector);
                matplot::hold(true);
                for(int i = 0; i < data.size(); i++){
                    plotting_vector[i] = (filter[i].imag());
                }
                matplot::plot(plotting_vector);
                matplot::hold(false);
                matplot::show();*/

                UTILITY_MATH::fftc2c(filter, filter_non_local_afr);
                /*for(int i = 0; i < data.size(); i++){
                    plotting_vector[i] = (filter_non_local_afr[i].real());
                }
                matplot::plot(plotting_vector);
                matplot::hold(true);
                for(int i = 0; i < data.size(); i++){
                    plotting_vector[i] = (filter_non_local_afr[i].imag());
                }
                matplot::plot(plotting_vector);
                matplot::hold(false);
                matplot::show();*/

                UTILITY_MATH::fftc2c(buffer1, buffer2);

                for(int i = 0; i < data.size(); i++){
                    buffer2[i] *= filter_non_local_afr[i];
                }

                UTILITY_MATH::ifftc2c(buffer2, buffer1);

                for (int i = 0; i < data.size(); i++) {
                    out[i] = buffer1[i].real();
                }
            }
        }
    };

    template<typename U>
    struct SincResFilter{
        constexpr static bool is_filter = true;
        GenericSignal<SimpleVecWrapper<U>, true> buffer1;
        GenericSignal<SimpleVecWrapper<U>, true> buffer2;
        double freq;
        bool is_low_pass;
        NonLocalNonOptFiltering<double,
                NonLocalFilteringType::Sinc> filter;

        template<Signal DataT, Signal OutT>
        void compute(const DataT & data, OutT & out, std::nullptr_t null){
            auto ext_info = UTILITY_MATH::circleExtend(data, *(buffer1.base->vec), freq);
            UTILITY_MATH::circleExtend(data, *(buffer2.base->vec), freq);
            //IC(ext_info.freq_idx, freq, ext_info.pad);
            filter.computeSincWithMonoFreq(buffer1, buffer2, ext_info.freq_idx, is_low_pass);
            for (int i = 0; i < data.size(); i++){
                out[i] = buffer2[i + ext_info.pad];
            }
        }
    };

    template<typename U>
    struct SincResLocalFilter{
        constexpr static bool is_filter = true;
        GenericSignal<SimpleVecWrapper<U>, true> buffer1;
        GenericSignal<SimpleVecWrapper<U>, true> buffer2;
        double freq;
        bool is_low_pass;
        double locality_coeff = 5.0;
        NonLocalNonOptFiltering<double,
                NonLocalFilteringType::SincLocal> filter;

        template<Signal DataT, Signal OutT>
        void compute(const DataT & data, OutT & out, std::nullptr_t null){
            //std::vector<std::complex>
            auto ext_info = UTILITY_MATH::circleExtend(data, *(buffer1.base->vec), freq);
            UTILITY_MATH::circleExtend(data, *(buffer2.base->vec), freq);
            //IC(ext_info.freq_idx, freq, ext_info.pad);
            filter.width_muller = locality_coeff;
            filter.computeSincWithMonoFreq(buffer1, buffer2, ext_info.freq_idx, true);
            if (is_low_pass){
                for (int i = 0; i < data.size(); i++){
                    out[i] = buffer2[i + ext_info.pad];
                }
            }
            else{
                for (int i = 0; i < data.size(); i++){
                    out[i] = data[i] - buffer2[i + ext_info.pad];
                }
            }
        }
    };

    template<typename U, InstAmplComputer<U> InstAmplComputerT, Signal InstFreqT, 
        Filter<U> FilterT, Integrator<U> IntegratorT, Derivator<U> DerivatorT>
    struct InstAmplNormalizatorNaiveDer{
        InstAmplComputerT * inst_ampl_computer;

        double min_ampl = 10.0;

        DerivatorT derivator;
        IntegratorT integrator;
        InstFreqT * inst_freq;
        FilterT * filter;

        GenericSignal<SimpleVecWrapper<U>, true> buffer2;

        InstAmplNormalizatorNaiveDer(InstAmplComputerT & inst_ampl_computer, FilterT & filter){
            this->inst_ampl_computer = &inst_ampl_computer;
            this->filter = &filter;
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        void compute(const DataT & data, OutT & out, ComputerBufferT & compute_buffer){
            if constexpr (InstAmplComputerT::is_used_external_inst_freq){
                inst_ampl_computer->inst_freq = inst_freq;
            }
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }
            
            inst_ampl_computer->compute(data, buffer2, &compute_buffer);
            //auto old_muller = filter->period_muller;
            //filter->period_muller = 2.0;
            filter->compute(buffer2, out, inst_freq);
            derivator.compute(data, compute_buffer, nullptr);
            //filter->period_muller = old_muller;

            auto b = data[0];
            auto avg = 0.0; 

            for (int i = 0; i < size; i++){
                avg += out[i] / size;
            }
            for (int i = 0; i < size; i++){
                if (out[i] > min_ampl){
                    compute_buffer[i] = compute_buffer[i] / (out[i] / avg);
                }
                else{
                    compute_buffer[i] = compute_buffer[i]  / (min_ampl / avg);
                }
            }

            integrator.compute(compute_buffer, out, nullptr);
            for (int i = 0; i < size; i++){
                out[i] += b;
            }
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        double computeGetAvg(const DataT & data, OutT & out, ComputerBufferT & compute_buffer){
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }
            
            inst_ampl_computer->compute(data, buffer2, &compute_buffer);
            //auto old_muller = filter->period_muller;
            //filter->period_muller = 2.0;
            filter->compute(buffer2, out, inst_freq);
            derivator.compute(data, compute_buffer, nullptr);
            //filter->period_muller = old_muller;

            auto b = data[0];
            auto avg = 0.0;
            for (int i = 0; i < size; i++){
                avg += out[i] / size;
            }
            for (int i = 0; i < size; i++){
                if (out[i] > min_ampl){
                    compute_buffer[i] = compute_buffer[i] / (out[i] / avg);
                }
                else{
                    compute_buffer[i] = compute_buffer[i] / (min_ampl / avg);
                }
            }
            integrator.compute(compute_buffer, out, nullptr);
            for (int i = 0; i < size; i++){
                out[i] += b;
            }

            return avg;
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        void computeWithExternalAvg(const DataT & data, OutT & out, ComputerBufferT & compute_buffer, double avg){
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }
            auto b = data[0];
            
            inst_ampl_computer->compute(data, buffer2, &compute_buffer);
            //auto old_muller = filter->period_muller;
            //filter->period_muller = 2.0;
            filter->compute(buffer2, out, inst_freq);
            derivator.compute(data, compute_buffer, nullptr);
            //filter->period_muller = old_muller;
            //out.show(PlottingKind::Simple);
            
            for (int i = 0; i < size; i++){
                if (out[i] > min_ampl){
                    compute_buffer[i] = compute_buffer[i] / (out[i] / avg);
                }
                else{
                    compute_buffer[i] = compute_buffer[i]  / (min_ampl / avg);
                }
            }
            integrator.compute(compute_buffer, out, nullptr);
            for (int i = 0; i < size; i++){
                out[i] += b;
            }
        }
    };

    
    template<typename U, InstAmplComputer<U> InstAmplComputerT, typename InstFreqT, Filter<U> FilterT>
    struct InstAmplNormalizatorNaive{
        InstAmplComputerT * inst_ampl_computer;

        double min_ampl = 10.0;

        InstFreqT * inst_freq;
        FilterT * filter;

        GenericSignal<SimpleVecWrapper<U>, true> buffer2;

        InstAmplNormalizatorNaive(InstAmplComputerT & inst_ampl_computer, FilterT & filter){
            this->inst_ampl_computer = &inst_ampl_computer;
            this->filter = &filter;
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        void compute(const DataT & data, OutT & out, ComputerBufferT & compute_buffer){
            if constexpr (InstAmplComputerT::is_used_external_inst_freq){
                inst_ampl_computer->inst_freq = inst_freq;
            }
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }
            
            inst_ampl_computer->compute(data, buffer2, &compute_buffer);
            //auto old_muller = filter->period_muller;
            //filter->period_muller = 2.0;
            filter->compute(buffer2, out, inst_freq);
            //filter->period_muller = old_muller;

            auto avg = 0.0;
            for (int i = 0; i < size; i++){
                avg += out[i] / size;
            }
            for (int i = 0; i < size; i++){
                if (out[i] > min_ampl){
                    compute_buffer[i] = data[i] / (out[i] / avg);
                }
                else{
                    compute_buffer[i] = data[i]  / (min_ampl / avg);
                }
            }

            for (int i = 0; i < size; i++){
                out[i] = compute_buffer[i];
            }
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        double computeGetAvg(const DataT & data, OutT & out, ComputerBufferT & compute_buffer){
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }
            
            inst_ampl_computer->compute(data, buffer2, &compute_buffer);
            //auto old_muller = filter->period_muller;
            //filter->period_muller = 2.0;
            filter->compute(buffer2, out, inst_freq);
            //filter->period_muller = old_muller;

            auto avg = 0.0;
            for (int i = 0; i < size; i++){
                avg += out[i] / size;
            }
            for (int i = 0; i < size; i++){
                if (out[i] > min_ampl){
                    compute_buffer[i] = data[i] / (out[i] / avg);
                }
                else{
                    compute_buffer[i] = data[i]  / (min_ampl / avg);
                }
            }
            for (int i = 0; i < size; i++){
                out[i] = compute_buffer[i];
            }

            return avg;
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        void computeWithExternalAvg(const DataT & data, OutT & out, ComputerBufferT & compute_buffer, double avg){
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }
            
            inst_ampl_computer->compute(data, buffer2, &compute_buffer);
            //auto old_muller = filter->period_muller;
            //filter->period_muller = 2.0;
            filter->compute(buffer2, out, inst_freq);
            //filter->period_muller = old_muller;
            //out.show(PlottingKind::Simple);
            
            for (int i = 0; i < size; i++){
                if (out[i] > min_ampl){
                    compute_buffer[i] = data[i] / (out[i] / avg);
                }
                else{
                    compute_buffer[i] = data[i]  / (min_ampl / avg);
                }
            }
            for (int i = 0; i < size; i++){
                out[i] = compute_buffer[i];
            }
        }
    };

    
    template<typename U, typename InstAmplNormalizer, typename InstFreqT, Filter<U> FilterT>
    struct InstAmplNormalizatorNaiveReqursive{
        InstAmplNormalizer * single_normalizer;
        InstFreqT * inst_freq;
        FilterT * filter;

        NP_DSP::ONE_D::GenericSignal<NP_DSP::ONE_D::SimpleVecWrapper<double>, true> mode;
        int iters_count = 50;

        InstAmplNormalizatorNaiveReqursive(InstAmplNormalizer & single_normalizer, 
            InstFreqT & inst_freq, FilterT & filter){
            this->single_normalizer = &single_normalizer;
            this->inst_freq = &inst_freq;
            this->filter = &filter;
        }

        template <Signal DataT, Signal OutT, Signal ComputeBufferT>
        void compute(const DataT & data, OutT & out, ComputeBufferT * compute_buffer){
            if (mode.size() != data.size()){
                mode.base->vec->clear();
                for (int i = 0; i < data.size(); i++){
                    mode.base->vec->push_back(0.0);
                }
            }
            
            filter->compute(data, mode, inst_freq);
            //mode.show(PlottingKind::Simple);
            for (int i = 0; i < data.size(); i++){
                mode[i] = data[i] - mode[i];
            }
            //mode.show(PlottingKind::Simple);
            single_normalizer->inst_freq = inst_freq;
            double avg = single_normalizer->computeGetAvg(mode, out, *compute_buffer);
            //out.show(PlottingKind::Simple);
            for (int i = 0; i < iters_count / 2; i++){
                single_normalizer->computeWithExternalAvg(out, mode, *compute_buffer, avg);
                //mode.show(PlottingKind::Simple);
                single_normalizer->computeWithExternalAvg(mode, out, *compute_buffer, avg);
                //mode.show(PlottingKind::Simple);
            }
        }
    };

    template<typename U, Integrator<U> Integrator, Derivator<U> Derivator, InstAmplComputer<U> InstAmplComputerT>
    struct InstAmplNormalizator{
        InstAmplComputerT * inst_ampl_computer;
        Derivator * derivator;
        Integrator * integrator;

        double min_ampl = 1.0;

        GenericSignal<SimpleVecWrapper<U>, true> buffer2;

        InstAmplNormalizator(Integrator & integrator, Derivator & derivator, InstAmplComputerT & inst_ampl_computer){
            this->inst_ampl_computer = &inst_ampl_computer;
            this->derivator = &derivator;
            this->integrator = &integrator;
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        void compute(const DataT & data, OutT & out, ComputerBufferT & compute_buffer){
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }
            
            inst_ampl_computer->compute(data, out, &compute_buffer);

            derivator->compute(data, compute_buffer, nullptr);
            auto avg = 0.0;
            auto b = data[0];
            for (int i = 0; i < size; i++){
                avg += out[i] / size;
            }
            for (int i = 0; i < size; i++){
                if (out[i] > min_ampl){
                    compute_buffer[i] /= (out[i] / avg);
                }
                else{
                    compute_buffer[i] /= (min_ampl / avg);
                }
            }
            integrator->compute(compute_buffer, out, nullptr);
            for (int i = 0; i < size; i++){
                out[i] += b;
            }
        }
    };

    template<typename U, Integrator<U> Integrator, Derivator<U> Derivator, 
        InstAmplComputer<U> InstAmplComputerT, bool is_double, 
            bool remove_big_der_before_derivating, Signal InstFreqT>
    struct InstAmplNormalizatorUsingInstFreq{
        InstAmplComputerT * inst_ampl_computer;
        Derivator * derivator;
        Integrator * integrator;
        InstFreqT * inst_freq;

        double min_ampl = 1.0;

        GenericSignal<SimpleVecWrapper<U>, true> buffer2;

        InstAmplNormalizatorUsingInstFreq(Integrator & integrator, Derivator & derivator, InstAmplComputerT & inst_ampl_computer){
            this->inst_ampl_computer = &inst_ampl_computer;
            this->derivator = &derivator;
            this->integrator = &integrator;
        }

        template<Signal DataT, Signal OutT, Signal ComputerBufferT>
        void compute(const DataT & data, OutT & out, ComputerBufferT & compute_buffer){
            auto size = data.size();
            if (buffer2.size() != size){
                buffer2.base->vec->clear();
                for(int i = 0; i < size; i++){
                    buffer2.base->vec->push_back(0);
                }
            }

            if constexpr (remove_big_der_before_derivating){
                /*std::vector<U> buffer;
                
                for (int i = 0; i < data.size(); i++){
                    double big_der;
                    if constexpr (is_double){
                        big_der = (data.interpolate
                            (i + 0.5/inst_freq_computer[i].second, SignalKind::Universal) -
                            data.interpolate(i - 0.5/inst_freq_computer[i].first, SignalKind::Universal)) 
                            * inst_freq_computer[i];
                    }
                    else{
                        big_der = (data.interpolate
                            (i + 0.5/inst_freq_computer[i], SignalKind::Universal) -
                            data.interpolate(i - 0.5/inst_freq_computer[i], SignalKind::Universal)) 
                            * inst_freq_computer[i];
                    }
                    buffer.push_back(big_der);
                    out[i] = data[i] - big_der;
                }
                auto b = out[0];
                derivator->compute(out, compute_buffer, nullptr);
                //todo resolve this conflict
                GenericSignal<SimpleVecWrapper<double>, true> computer_buffer2;
                for (int i = 0; i < data.size(); i++){
                    computer_buffer2.base->vec->push_back(0.0);
                }

                inst_ampl_computer->compute(data, out, &computer_buffer2);
                auto avg = 0.0;
                for (int i = 0; i < size; i++){
                    avg += out[i] / size;
                }
                for (int i = 0; i < size; i++){
                    if (out[i] > min_ampl){
                        compute_buffer[i] /= (out[i] / avg);
                    }
                }
                
                integrator->compute(compute_buffer, out, nullptr);
                for (int i = 0; i < size; i++){
                    out[i] += buffer[i];
                }*/
            }
            else{
                std::vector<U> buffer;
                inst_ampl_computer->compute(data, out, &compute_buffer);
                derivator->compute(data, compute_buffer, nullptr);
                //IC(*(data.base->vec));
                //data.show(PlottingKind::Simple);
                //IC(*(out.base->vec));
                //out.show(PlottingKind::Simple);
                //IC(*(compute_buffer.base->vec));
                //compute_buffer.show(PlottingKind::Simple);

                for (int i = 0; i < data.size(); i++){
                    double big_der;
                    if constexpr (is_double){
                        big_der = (data.interpolate
                            (i + 0.5/(*inst_freq)[i].second, SignalKind::Universal) -
                            data.interpolate(i - 0.5/(*inst_freq)[i].first, SignalKind::Universal)) 
                            / (0.5 / (*inst_freq)[i].first + 0.5 / (*inst_freq)[i].second);
                    }
                    else{
                        big_der = (data.interpolate
                            (i + 0.5/(*inst_freq)[i], SignalKind::Universal) -
                            data.interpolate(i - 0.5/(*inst_freq)[i], SignalKind::Universal)) 
                            * (*inst_freq)[i];
                    }
                    buffer.push_back(big_der);
                    compute_buffer[i] = compute_buffer[i] - big_der;
                }
                //IC(*(inst_freq->base->vec));
                //IC(buffer);
                //matplot::plot(buffer);
                //matplot::show();

                //IC(*(compute_buffer.base->vec));
                //compute_buffer.show(PlottingKind::Simple);

                auto avg = 0.0;
                auto b = data[0];
                for (int i = 0; i < size; i++){
                    if (out[i] > min_ampl){
                        avg += out[i];
                    }
                    else{
                        avg += min_ampl;
                    }
                }
                /*auto size_norm = 0;
                for (int i = 0; i < size; i++){
                    if (out[i] > min_ampl){
                        size_norm++;
                    }
                }*/
                avg /= size;
                for (int i = 0; i < size; i++){
                    if (out[i] > min_ampl){
                        compute_buffer[i] /= (out[i] / avg);
                    }
                    else{
                        compute_buffer[i] /= (min_ampl / avg);
                    }
                }
                //IC(*(compute_buffer.base->vec));
                //compute_buffer.show(PlottingKind::Simple);
                for (int i = 0; i < data.size(); i++){
                    compute_buffer[i] += buffer[i]; 
                }
                //IC(*(compute_buffer.base->vec));
                //compute_buffer.show(PlottingKind::Simple); 
                integrator->compute(compute_buffer, out, nullptr);
                //IC(*(out.base->vec));
                //out.show(PlottingKind::Simple);
                for (int i = 0; i < size; i++){
                    out[i] += b;
                }
                //IC(*(out.base->vec));
                //out.show(PlottingKind::Simple);
            }
        }
    };

    enum class PhaseComputingKind { extremums_based_non_opt, arctg_scaled };

    
    template<typename U, Filter<U> FilterT, InstFreqComputer<U> InstFreqComputerT,
        PhaseComputer<U> PhaseComputerT, InstFreqComputer<U> InstFreqComputerForModeT, 
        PhaseComputer<U> PhaseComputerForModeT>
    //INST_FREQ_COMPUTERS::InstFreqDerivativeBasedKind inst_freq_k>
    struct OptPeriodBasedFilter {
        using AdditionalDataType = SignalPrototype<U>;

        //SimpleVecWrapper<T> mode_base;
        //SimpleVecWrapper<T> inst_freq_buffer2_base;
        size_t max_iter_number = 1000;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;
        BufferT mode;
        BufferT inst_freq_buffer2;

        double error_threshold = 0.001;
        size_t iter_number = 0;
        size_t good_iter_number = 0;
        size_t true_iter_number = 0;

        std::vector<double> inst_freq_cache;
        double error_old;

        //constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
        constexpr static bool is_filter = true;

        FilterT * filter;
        InstFreqComputerT * inst_freq_computer;
        PhaseComputerT * phase_computer;
        InstFreqComputerForModeT * inst_freq_computer_for_mode;
        PhaseComputerForModeT * phase_computer_for_mode;

        OptPeriodBasedFilter(FilterT & filter,
            InstFreqComputerT & inst_freq_computer,
            PhaseComputerT & phase_computer, 
            InstFreqComputerForModeT & inst_freq_computer_for_mode, 
            PhaseComputerForModeT & phase_computer_for_mode) {
            this->filter = &filter;
            this->inst_freq_computer = &inst_freq_computer;
            this->phase_computer = &phase_computer;
            this->inst_freq_computer_for_mode = &inst_freq_computer_for_mode;
            this->phase_computer_for_mode = &phase_computer_for_mode;
        }

        ~OptPeriodBasedFilter() {
            //delete mode;
            //delete inst_freq_buffer2;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        bool computeIter(const DataType& data, OutType& out, InstFreqType& inst_freq_buffer) {
            using T = typename OutType::SampleType;
            //IC(iter_number);
            iter_number++;
            true_iter_number++;
            if (true_iter_number > max_iter_number){
                for (int i = 0; i < data.size(); i++) {
                    inst_freq_buffer[i] = inst_freq_cache[i];
                }
                filter->compute(data, out, &inst_freq_buffer);
                return false;
            }
            filter->compute(data, out, &inst_freq_buffer);

            //after computing filtering
            if (mode.size() != data.size()) {
                static_cast<SimpleVecWrapper<T> *>(mode.base)->vec->clear();
                for (auto i = 0; i < data.size(); i++) {
                    static_cast<SimpleVecWrapper<T> *>(mode.base)->vec->push_back(data[i] - out[i]);
                }
            } else {
                for (auto i = 0; i < data.size(); i++) {
                    mode[i] = data[i] - out[i];
                }
            }
            if (inst_freq_buffer2.size() != inst_freq_buffer.size()) {
                static_cast<SimpleVecWrapper<T> *>(inst_freq_buffer2.base)->vec->clear();
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    static_cast<SimpleVecWrapper<T> *>(inst_freq_buffer2.base)->vec->push_back(0.0);
                }
            }

            if (inst_freq_computer_for_mode->is_phase_based()) {
                //todo
                /*std::unreachable();*/
            }
            else {
                if constexpr (std::is_same_v<typename InstFreqComputerForModeT::AdditionalDataType, GENERAL::Nil>) {
                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, nullptr);
                } else {
                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, &out);
                }
            }

            double error = UTILITY_MATH::signalsL2NormedDistance
                    <double, InstFreqType, decltype(inst_freq_buffer2)>
                    (inst_freq_buffer, inst_freq_buffer2);

            //double error_old_new = UTILITY_MATH::signalsL2Distance
            //<double, InstFreqType, decltype(inst_freq_cache)>
            //(inst_freq_buffer, inst_freq_cache);
            //IC(error, error_old);
            if (error_old < error) {
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 5 &&  iter_number > 50) {
                    for (int i = 0; i < data.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }
                    filter->compute(data, out, &inst_freq_buffer);
                    return false;
                }
                if (iter_number - good_iter_number > 4) {
                    for (int i = 0; i < inst_freq_buffer.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }

                    iter_number = good_iter_number;
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = data[i] - mode[i];
                    }
                    //IC(error);
                    auto res = computeIter(data, out, inst_freq_buffer);
                    for (int i = 0; i < data.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }
                    filter->compute(data, out, &inst_freq_buffer);
                    return res;
                    error = error_old;
                }
            
            } else {
                //IC(error, good_iter_number, iter_number, true_iter_number);
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 4) {
                    for (int i = 0; i < data.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }
                    filter->compute(data, out, &inst_freq_buffer);
                    return false;
                }
                
                good_iter_number = iter_number;
                error_old = error;
            }

            if (error < error_threshold) {
                for (int i = 0; i < data.size(); i++) {
                    inst_freq_buffer[i] = inst_freq_cache[i];
                }
                filter->compute(data, out, &inst_freq_buffer);
                return false;
            }
            if (iter_number == good_iter_number) {
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    inst_freq_cache[i] = inst_freq_buffer[i];
                }
            }

            for (int i = 0; i < inst_freq_buffer.size(); i += std::rand()%3) {
                //= i + (std::rand() + 1) % (inst_freq_buffer.size()/5)){
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                inst_freq_buffer[i] =
                        (inst_freq_buffer2[i] * coeff2 + inst_freq_buffer[i] * coeff1)
                        / (coeff1 + coeff2);
            }

            for (int i = 0; i < data.size(); i++) {
                out[i] = data[i] - mode[i];
            }

            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        void compute(const DataType& data, OutType& out, InstFreqType * inst_freq_buffer) {
            iter_number = 0;
            true_iter_number = 0;
            good_iter_number = 0;
            using T = typename OutType::SampleType;
            //inst_freq_computer->compute(data, inst_freq_buffer, out);
            if (inst_freq_computer->is_phase_based()) {
                //todo
                /*std::unreachable();*/
            }
            if constexpr (std::is_same_v<typename InstFreqComputerT::AdditionalDataType, GENERAL::Nil>) {
                inst_freq_computer->compute(data, *inst_freq_buffer, nullptr);
            } else {
                inst_freq_computer->compute(data, *inst_freq_buffer, &out);
            }

            filter->compute(data, out, inst_freq_buffer);

            //after computing filtering
            if (mode.size() != data.size()) {
                static_cast<SimpleVecWrapper<T> *>(mode.base)->vec->clear();
                for (auto i = 0; i < data.size(); i++) {
                    static_cast<SimpleVecWrapper<T> *>(mode.base)->vec->push_back(data[i] - out[i]);
                }
            } else {
                for (auto i = 0; i < data.size(); i++) {
                    mode[i] = data[i] - out[i];
                }
            }
            if (inst_freq_buffer2.size() != inst_freq_buffer->size()) {
                static_cast<SimpleVecWrapper<T> *>(inst_freq_buffer2.base)->vec->clear();
                for (int i = 0; i < inst_freq_buffer->size(); i++) {
                    static_cast<SimpleVecWrapper<T> *>(inst_freq_buffer2.base)->vec->push_back(0.0);
                }
            }


            //inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, out);
            if constexpr (std::is_same_v<typename InstFreqComputerForModeT::AdditionalDataType, GENERAL::Nil>) {
                inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, nullptr);
            } else {
                inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, &out);
            }

            double error = UTILITY_MATH::signalsL2NormedDistance
                    <double, InstFreqType, decltype(inst_freq_buffer2)>
                    (*inst_freq_buffer, inst_freq_buffer2);
            error_old = error;
            //IC(error);

            if (error < error_threshold) {
                //for (int i = 0; i < data.size(); i++) {
                    //(*inst_freq_buffer)[i] = inst_freq_cache[i];
                //}
                filter->compute(data, out, inst_freq_buffer);
                return;
            }
            inst_freq_cache.clear();
            for (int i = 0; i < inst_freq_buffer->size(); i++) {
                inst_freq_cache.push_back((*inst_freq_buffer)[i]);
            }
            for (int i = 0; i < inst_freq_buffer->size(); i += std::rand() % 3) {
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                (*inst_freq_buffer)[i] =
                        (inst_freq_buffer2[i] * coeff2 + (*inst_freq_buffer)[i] * coeff1)
                        / (coeff1 + coeff2);
            }
            while (computeIter(data, out, *inst_freq_buffer)) {
            }
        }
    };


    
    template<typename U, Filter<U> FilterT, InstFreqComputer<U> InstFreqComputerT,
        PhaseComputer<U> PhaseComputerT, InstFreqComputer<U> InstFreqComputerForModeT, 
        PhaseComputer<U> PhaseComputerForModeT>
    struct OptPeriodBasedFilterInstFreqDouble {
        using AdditionalDataType = SignalPrototype<U>;


        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;
        BufferT mode;
        GenericSignal<SimpleVecWrapper<std::pair<U, U>>, true> inst_freq_buffer2;

        size_t max_iter_number = 1000;
        double error_threshold = 0.1;
        double error_treshhold_muller = 0.01;
        size_t iter_number = 0;
        size_t good_iter_number = 0;
        size_t true_iter_number = 0;

        std::vector<std::pair<double, double>> inst_freq_cache;
        double error_old;

        constexpr static bool is_filter = true;

        FilterT * filter;
        InstFreqComputerT * inst_freq_computer;
        PhaseComputerT * phase_computer;
        InstFreqComputerForModeT * inst_freq_computer_for_mode;
        PhaseComputerForModeT * phase_computer_for_mode;

        OptPeriodBasedFilterInstFreqDouble(FilterT & filter,
            InstFreqComputerT & inst_freq_computer,
            PhaseComputerT & phase_computer, 
            InstFreqComputerForModeT & inst_freq_computer_for_mode, 
            PhaseComputerForModeT & phase_computer_for_mode) {
            this->filter = &filter;
            this->inst_freq_computer = &inst_freq_computer;
            this->phase_computer = &phase_computer;
            this->inst_freq_computer_for_mode = &inst_freq_computer_for_mode;
            this->phase_computer_for_mode = &phase_computer_for_mode;
        }

        ~OptPeriodBasedFilterInstFreqDouble() {
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        bool computeIter(const DataType& data, OutType& out, InstFreqType& inst_freq_buffer) {
            using T = typename OutType::SampleType;
            iter_number++;
            true_iter_number++;

            if (true_iter_number > max_iter_number){
                for (int i = 0; i < data.size(); i++) {
                    inst_freq_buffer[i] = inst_freq_cache[i];
                }
                filter->compute(data, out, &inst_freq_buffer);
                return false;
            }

            filter->compute(data, out, &inst_freq_buffer);

            if (mode.size() != data.size()) {
                (mode.base)->vec->clear();
                for (auto i = 0; i < data.size(); i++) {
                    (mode.base)->vec->push_back(data[i] - out[i]);
                }
            } else {
                for (auto i = 0; i < data.size(); i++) {
                    mode[i] = data[i] - out[i];
                }
            }
            if (inst_freq_buffer2.size() != inst_freq_buffer.size()) {
                (inst_freq_buffer2.base)->vec->clear();
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    (inst_freq_buffer2.base)->vec->push_back({0.0, 0.0});
                }
            }

            if (inst_freq_computer_for_mode->is_phase_based()) {
                /*std::unreachable();*/
            }
            else {
                if constexpr (std::is_same_v<typename InstFreqComputerForModeT::AdditionalDataType, GENERAL::Nil>) {
                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, nullptr);
                } else {
                    inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, &out);
                }
            }

            double error = UTILITY_MATH::signalsL2NormedDistanceDouble
                    <double, InstFreqType, decltype(inst_freq_buffer2)>
                    (inst_freq_buffer, inst_freq_buffer2);

            if (error_old < error) {
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 5 && iter_number > 50) {
                    for (int i = 0; i < data.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }
                    filter->compute(data, out, &inst_freq_buffer);

                    return false;
                }
                if (iter_number - good_iter_number > 4) {
                    for (int i = 0; i < inst_freq_buffer.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }

                    iter_number = good_iter_number;
                    for (int i = 0; i < data.size(); i++) {
                        out[i] = data[i] - mode[i];
                    }
                    auto res = computeIter(data, out, inst_freq_buffer);
                    //for (int i = 0; i < data.size(); i++) {
                        //inst_freq_buffer[i] = inst_freq_cache[i];
                    //}
                    //filter->compute(data, out, &inst_freq_buffer);
                    return res;
                    error = error_old;
                }
            
            } else {
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 4) {
                    for (int i = 0; i < data.size(); i++) {
                        inst_freq_buffer[i] = inst_freq_cache[i];
                    }
                    filter->compute(data, out, &inst_freq_buffer);
                    return false;
                }
                
                good_iter_number = iter_number;
                error_old = error;
            }

            if (error < error_threshold) {
                for (int i = 0; i < data.size(); i++) {
                    out[i] = data[i] - mode[i];
                }
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    inst_freq_cache[i] = inst_freq_buffer[i];
                }
                return false;
            }

            if (iter_number == good_iter_number) {
                for (int i = 0; i < inst_freq_buffer.size(); i++) {
                    inst_freq_cache[i] = inst_freq_buffer[i];
                }
            }

            for (int i = 0; i < inst_freq_buffer.size(); i += std::rand() % 3) {
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                inst_freq_buffer[i].first =
                        (inst_freq_buffer2[i].first * coeff1 + inst_freq_buffer[i].first * coeff2)
                        / (coeff1 + coeff2);
                inst_freq_buffer[i].second =
                        (inst_freq_buffer2[i].second * coeff1 + inst_freq_buffer[i].second * coeff2)
                        / (coeff1 + coeff2);
            }

            for (int i = 0; i < data.size(); i++) {
                out[i] = data[i] - mode[i];
            }

            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        void compute(const DataType& data, OutType& out, InstFreqType * inst_freq_buffer) {
            iter_number = 0;
            true_iter_number = 0;
            good_iter_number = 0;
            using T = typename OutType::SampleType;
            if (inst_freq_computer->is_phase_based()) {
                /*std::unreachable();*/
            }
            if constexpr (std::is_same_v<typename InstFreqComputerT::AdditionalDataType, GENERAL::Nil>) {
                inst_freq_computer->compute(data, *inst_freq_buffer, nullptr);
            } else {
                inst_freq_computer->compute(data, *inst_freq_buffer, &out);
            }

            filter->compute(data, out, inst_freq_buffer);

            if (mode.size() != data.size()) {
                (mode.base)->vec->clear();
                for (auto i = 0; i < data.size(); i++) {
                    (mode.base)->vec->push_back(data[i] - out[i]);
                }
            } else {
                for (auto i = 0; i < data.size(); i++) {
                    mode[i] = data[i] - out[i];
                }
            }
            if (inst_freq_buffer2.size() != inst_freq_buffer->size()) {
                (inst_freq_buffer2.base)->vec->clear();
                for (int i = 0; i < inst_freq_buffer->size(); i++) {
                    (inst_freq_buffer2.base)->vec->push_back({0.0, 0.0});
                }
            }

            if constexpr (std::is_same_v<typename InstFreqComputerForModeT::AdditionalDataType, GENERAL::Nil>) {
                inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, nullptr);
            } else {
                inst_freq_computer_for_mode->compute(mode, inst_freq_buffer2, &out);
            }

            double der_avg = 0.0;
            for (int i = 1; i < data.size(); i++){
                der_avg += std::abs((*inst_freq_buffer)[i].first) + std::abs((*inst_freq_buffer)[i].second);
            }
            der_avg = der_avg / data.size();

            error_threshold = der_avg * error_treshhold_muller;
            //todo


            double error = UTILITY_MATH::signalsL2NormedDistanceDouble
                    <double, InstFreqType, decltype(inst_freq_buffer2)>
                    (*inst_freq_buffer, inst_freq_buffer2);
            error_old = error;
            
            if (error < error_threshold) {
                for (int i = 0; i < data.size(); i++) {
                    out[i] = data[i] - mode[i];
                }
                return;
            }
            inst_freq_cache.clear();
            for (int i = 0; i < inst_freq_buffer->size(); i++) {
                inst_freq_cache.push_back((*inst_freq_buffer)[i]);
            }
            for (int i = 0; i < inst_freq_buffer->size(); i++) {
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                (*inst_freq_buffer)[i].first =
                        (inst_freq_buffer2[i].first * coeff1 + (*inst_freq_buffer)[i].first * coeff2)
                        / (coeff1 + coeff2);
                (*inst_freq_buffer)[i].second =
                        (inst_freq_buffer2[i].second * coeff1 + (*inst_freq_buffer)[i].second * coeff2)
                        / (coeff1 + coeff2);
            }
            while (computeIter(data, out, *inst_freq_buffer)) {
            }
        }
    };


    
    template<typename U, Integrator<U> IntegratorT, Derivator<U> DerivatorT, Filter<U> FilterT, 
        InstFreqComputer<U> InstFreqComputerT, 
            PhaseComputer<U> PhaseComputerT, InstAmplComputer<U> InstAmplComputerT,
                InstFreqComputer<U> InstFreqComputerForModeT, 
                    PhaseComputer<U> PhaseComputerForModeT>
    struct RecursiveFilterInstAmplChanges{
        using AdditionalDataType = SignalPrototype<U>;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;

        BufferT prediction_mode;
        BufferT prediction_mode_inst_freq;
        BufferT prediction_mode_ampl;
        BufferT new_result;

        static BufferT & prediction_phase(){
            static BufferT * prediction_mode_phase = new BufferT;
            return *prediction_mode_phase;
        }

        //BufferT computer_buffer;

        std::vector<double> result_cache;
        
        constexpr static bool is_filter = true;

        double error;
        double error_treshhold = 0.1;
        double error_treshhold_muller = 0.01;

        size_t good_iter_number = 0;
        size_t true_iter_number = 0;
        size_t iter_number = 0;
        size_t max_iter_number = 1000;

        double error_old = 0.0;


        FilterT * filter;
        InstFreqComputerT * inst_freq_computer;
        InstFreqComputerForModeT * inst_freq_computer_for_mode;
        PhaseComputerT * phase_computer;
        PhaseComputerForModeT * phase_computer_for_mode;
        InstAmplComputerT * inst_ampl_computer;
        InstAmplNormalizator<double, IntegratorT, DerivatorT, InstAmplComputerT> 
            * inst_ampl_normalizer;

        RecursiveFilterInstAmplChanges(IntegratorT & integrator, DerivatorT & derivator, 
            FilterT & filter, InstFreqComputerT & inst_freq_computer,
                PhaseComputerT & phase_computer, InstAmplComputerT & inst_ampl_computer, 
                    InstFreqComputerForModeT & inst_freq_computer_for_mode,
                        PhaseComputerForModeT & phase_computer_for_mode){
            this->filter = &filter;
            this->inst_freq_computer = &inst_freq_computer;
            this->inst_freq_computer_for_mode = &inst_freq_computer_for_mode;
            this->phase_computer = &phase_computer;
            this->phase_computer_for_mode = &phase_computer_for_mode;
            inst_ampl_normalizer = new InstAmplNormalizator
                <double, IntegratorT, DerivatorT, InstAmplComputerT> (integrator, derivator, inst_ampl_computer);
        }

        ~RecursiveFilterInstAmplChanges(){
            delete inst_ampl_normalizer;
        }

        template<Signal DataType, Signal OutType, Signal OldResultType>
        bool computeIter(const DataType& data, OutType& result_buffer, OldResultType* computer_buffer){
            iter_number++;
            true_iter_number++;

            if (true_iter_number > max_iter_number){
                for(int i = 0; i < data.size(); i++){
                    result_buffer[i] = result_cache[i];
                }
                return false;
            }

            using T = typename OutType::SampleType;

            if(result_cache.size()!=data.size()){
                result_cache.clear();
                for(auto i = 0; i < data.size(); i++){
                    result_cache.push_back(0.0);
                }
            }
            
            if(prediction_mode.size() != data.size()){
                prediction_mode.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_inst_freq.size() != data.size()){
                prediction_mode_inst_freq.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_inst_freq.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_ampl.size() != data.size()){
                prediction_mode_ampl.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_ampl.base->vec->push_back(0.0);
                }
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_inst_freq[i] = 0.0;
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_ampl[i] = 0.0;
            }
            for (auto i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = 0.0;
            }
            for(auto i = 0; i < data.size(); i++){
                prediction_mode[i] = data[i] - result_buffer[i];
            }

            if constexpr (InstFreqComputerForModeT::is_phase_based()){
                if (prediction_phase().size()!= data.size){
                    prediction_phase().base->vec->clear();
                    for (auto i = 0; i < data.size(); i++){
                        prediction_phase().base->vec->push_back(0.0);
                    }
                }
                for (auto i = 0; i < data.size(); i++){
                    prediction_phase()[i] = 0.0;
                }
                phase_computer->compute(prediction_mode, prediction_phase(), computer_buffer);
                inst_freq_computer_for_mode->compute(prediction_phase(), prediction_mode_inst_freq, computer_buffer);
            }
            else{
                inst_freq_computer_for_mode->compute(prediction_mode, prediction_mode_inst_freq, computer_buffer);
            }
            
            
            inst_ampl_normalizer->compute(prediction_mode, *computer_buffer, prediction_mode_ampl);

            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = result_buffer[i] + (*computer_buffer)[i];
            }

            filter->compute(*computer_buffer, prediction_mode, &prediction_mode_inst_freq);
            error = UTILITY_MATH::signalsL2NormedDistance
                    <double, decltype(result_buffer), decltype(prediction_mode)>
                    (result_buffer, prediction_mode);

            if (error_old <= error){
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    return false;
                }
                if(iter_number - good_iter_number > 9){
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    iter_number = good_iter_number;
                    return computeIter(data, result_buffer, computer_buffer);
                }
            }
            else if(error < error_old){
                //IC(error_old, error, error_treshhold);
                error_old = error;
                good_iter_number = iter_number;
                
                for (int i = 0; i < data.size(); i++){
                    result_cache[i] = result_buffer[i];
                }
                if(static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    return false;
                }

            }

            for(int i = 0; i < data.size(); i += std::rand() % 3){
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                result_buffer[i] = (result_buffer[i] * coeff1 + prediction_mode[i] * coeff2) / (coeff1 + coeff2);
            }
            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqT>
        void compute(const DataType& data, OutType& result_buffer, InstFreqT* inst_freq){
            iter_number = 0;
            good_iter_number = 0;
            true_iter_number = 0;
            
            //this->inst_freq_computer->compute(data, *inst_freq, &result_buffer);
            this->filter->compute(data, result_buffer, inst_freq);

            double der_avg = 0.0;
            auto sum = 0.0;
            for (int i = 1; i < data.size(); i++){
                der_avg += std::abs((data[i] - result_buffer[i])-(data[i-1] - result_buffer[i-1])) / data.size();
                sum += std::abs(data[i]);
            }
            error_old = sum * data.size();
            der_avg = der_avg / data.size();

            error_treshhold = der_avg * error_treshhold_muller;
            
            computeIter(data, result_buffer, inst_freq);

            while(error > error_treshhold){
                if(!computeIter(data, result_buffer, inst_freq)) break;
            }
        }   
    };



    
    template<typename U, Filter<U> FilterFirstT, 
        Filter<U> FilterSeconT>
    struct CascadeFilter{
        using AdditionalDataType = SignalPrototype<double>;
        using IdxType = size_t;
        constexpr static bool is_filter = true;

        FilterFirstT * filter_first;
        FilterSeconT * filter_second;

        double period_muller = 1.0; //todo remove cratch

        CascadeFilter(FilterFirstT & filter_first,
            FilterSeconT & filter_second){

            this->filter_first = &filter_first;
            this->filter_second = &filter_second;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqType>
        void compute(const DataType& data, OutType& out, InstFreqType * inst_freq){
            using T = typename OutType::SampleType;
            
            filter_first->compute(data, out, inst_freq);
            auto mode_val = [&](size_t idx){
                return data[idx] - out[idx];
            };
            auto mode_size = [&](){
                return data.size();
            };
            ExpressionWrapper<T, size_t, 
                decltype(mode_val), GENERAL::Nil, decltype(mode_size), false> 
                    mode_expr(mode_val, mode_size);
            GenericSignal<decltype(mode_expr), false> mode(mode_expr);

            GenericSignal<SimpleVecWrapper<T>, true> buffer;

            for(size_t i = 0 ; i < data.size(); i++){
                buffer.base->vec->push_back(0.);
            }

            filter_second->compute(mode, buffer, inst_freq);

            for (int i = 0; i < data.size(); i++){
                out[i] = data[i] - (mode[i] - buffer[i]);
            }
        }
    };

    
    template<typename U, Integrator<U> IntegratorT, Derivator<U> DerivatorT, Filter<U> FilterT, 
        InstAmplComputer<U> InstAmplComputerT>
    struct RecursiveFilterInstAmplChangesWithConstInstFreq{
        using AdditionalDataType = SignalPrototype<U>;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;

        BufferT prediction_mode;
        BufferT prediction_mode_inst_freq;
        BufferT prediction_mode_ampl;
        BufferT new_result;

        static BufferT & prediction_phase(){
            static BufferT * prediction_mode_phase = new BufferT;
            return *prediction_mode_phase;
        }

        //BufferT computer_buffer;

        std::vector<double> result_cache;
        
        constexpr static bool is_filter = true;

        double error;
        double error_treshhold = 0.1;
        double error_treshhold_muller = 0.01;

        size_t good_iter_number = 0;
        size_t true_iter_number = 0;
        size_t iter_number = 0;
        size_t max_iter_number = 10000;

        double error_old = 0.0;



        FilterT * filter;
        InstAmplComputerT * inst_ampl_computer;
        InstAmplNormalizatorNaive<double, InstAmplComputerT, BufferT, FilterT> 
            * inst_ampl_normalizer;

        RecursiveFilterInstAmplChangesWithConstInstFreq(IntegratorT & integrator, DerivatorT & derivator, 
            FilterT & filter, InstAmplComputerT & inst_ampl_computer){
            this->filter = &filter;
            inst_ampl_normalizer = new InstAmplNormalizatorNaive
                <double, InstAmplComputerT, BufferT, FilterT> (inst_ampl_computer, filter);
        }

        ~RecursiveFilterInstAmplChangesWithConstInstFreq(){
            delete inst_ampl_normalizer;
        }

        template<Signal DataType, Signal OutType, Signal OldResultType>
        bool computeIter(const DataType& data, OutType& result_buffer, OldResultType* computer_buffer){
            iter_number++;
            true_iter_number++;
            using T = typename OutType::SampleType;

            if (true_iter_number > max_iter_number){
                for(int i = 0; i < data.size(); i++){
                    result_buffer[i] = result_cache[i];
                }
                iter_number = good_iter_number;
                for(int i = 0; i < data.size(); i++){
                    (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                }
                return false;
            }
            //result_buffer.show(PlottingKind::Simple);

            if(result_cache.size()!=data.size()){
                result_cache.clear();
                for(auto i = 0; i < data.size(); i++){
                    result_cache.push_back(0.0);
                }
            }
            
            if(prediction_mode.size() != data.size()){
                prediction_mode.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_inst_freq.size() != data.size()){
                prediction_mode_inst_freq.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_inst_freq.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_ampl.size() != data.size()){
                prediction_mode_ampl.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_ampl.base->vec->push_back(0.0);
                }
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_inst_freq[i] = 0.0;
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_ampl[i] = 0.0;
            }
            for(auto i = 0; i < data.size(); i++){
                prediction_mode[i] = data[i] - result_buffer[i];
            }

            //prediction_mode.show(PlottingKind::Simple);

            for(int i = 0; i < data.size(); i++){
                prediction_mode_inst_freq[i] = (*computer_buffer)[i];
            }
            //auto label = 0;
            //IC(++label);
            //todo set inst_freq to inst_ampl_normalizer
            inst_ampl_normalizer->inst_freq = &prediction_mode_inst_freq;
            inst_ampl_normalizer->compute(prediction_mode, *computer_buffer, prediction_mode_ampl);

            //computer_buffer->show(PlottingKind::Simple);

            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = result_buffer[i] + (*computer_buffer)[i];
            }
            //IC(++label);
            //computer_buffer->show(PlottingKind::Simple);
            filter->compute(*computer_buffer, prediction_mode, &prediction_mode_inst_freq);
            //prediction_mode.show(PlottingKind::Simple);
            //IC(++label);
            error = UTILITY_MATH::signalsL2NormedDistance
                    <double, decltype(result_buffer), decltype(prediction_mode)>
                    (result_buffer, prediction_mode);
            //IC(error);

            if (error_old <= error){
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return false;
                }
                if(iter_number - good_iter_number > 9){
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    iter_number = good_iter_number;
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return computeIter(data, result_buffer, computer_buffer);
                }
            }
            else if(error < error_old){
                //IC(error_old, error, error_treshhold);
                error_old = error;
                good_iter_number = iter_number;
                
                for (int i = 0; i < data.size(); i++){
                    result_cache[i] = result_buffer[i];
                }
                if(static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        (*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return false;
                }
            }

            for(int i = 0; i < data.size(); i += std::rand() % 3){
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                result_buffer[i] = (result_buffer[i] * coeff1 + prediction_mode[i] * coeff2) / (coeff1 + coeff2);
            }
            //result_buffer.show(PlottingKind::Simple);
            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = prediction_mode_inst_freq[i];
            }
            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqT>
        void compute(const DataType& data, OutType& result_buffer, InstFreqT* inst_freq){
            iter_number = 0;
            good_iter_number = 0;
            true_iter_number = 0;
            
            //this->inst_freq_computer->compute(data, *inst_freq, &result_buffer);
            this->filter->compute(data, result_buffer, inst_freq);

            double der_avg = 0.0;
            auto sum = 0.0;
            for (int i = 1; i < data.size(); i++){
                der_avg += std::abs((data[i] - result_buffer[i])-(data[i-1] - result_buffer[i-1])) / data.size();
                sum += std::abs(data[i]);
            }
            error_old = sum * data.size();
            der_avg = der_avg / data.size();

            error_treshhold = der_avg * error_treshhold_muller;
            
            computeIter(data, result_buffer, inst_freq);

            while(error > error_treshhold){
                if(!computeIter(data, result_buffer, inst_freq)) break;
            }
        }   
    };

    
    template<typename U, Integrator<U> IntegratorT, Derivator<U> DerivatorT, Filter<U> FilterT, 
        InstAmplComputer<U> InstAmplComputerT>
    struct RecursiveFilterInstAmplChangesWithConstInstFreqDouble{
        //todo
        using AdditionalDataType = SignalPrototype<U>;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;

        BufferT prediction_mode;
        GenericSignal<SimpleVecWrapper<std::pair<double, double>>, true> prediction_mode_inst_freq;
        BufferT prediction_mode_ampl;
        BufferT new_result;

        static BufferT & prediction_phase(){
            static BufferT * prediction_mode_phase = new BufferT;
            return *prediction_mode_phase;
        }

        //BufferT computer_buffer;

        std::vector<double> result_cache;
        
        
        constexpr static bool is_filter = true;

        double error;
        double error_treshhold = 0.1;
        double error_treshhold_muller = 0.01;

        size_t good_iter_number = 0;
        size_t true_iter_number = 0;
        size_t iter_number = 0;
        size_t max_iter_number = 1000;

        double error_old = 0.0;


        FilterT * filter;
        InstAmplComputerT * inst_ampl_computer;
        InstAmplNormalizatorNaive
            <double, InstAmplComputerT, decltype(prediction_mode_inst_freq), FilterT> 
                * inst_ampl_normalizer;

        /*InstAmplNormalizatorNaiveReqursive
            <double, 
            InstAmplNormalizatorNaive
            <double, InstAmplComputerT, decltype(prediction_mode_inst_freq), 
                FilterT>, decltype(prediction_mode_inst_freq), 
                    FilterT> * inst_ampl_normalizer;*/
        //(*naive_normalizer, inst_freq_buffer, non_opt_filter);

        RecursiveFilterInstAmplChangesWithConstInstFreqDouble(IntegratorT & integrator, DerivatorT & derivator, 
            FilterT & filter, InstAmplComputerT & inst_ampl_computer){
            this->filter = &filter;
            inst_ampl_normalizer = new InstAmplNormalizatorNaive
                <double, InstAmplComputerT, decltype(prediction_mode_inst_freq), FilterT> 
                    (inst_ampl_computer, filter);
            /*inst_ampl_normalizer = new InstAmplNormalizatorNaiveReqursive
            <double, 
            InstAmplNormalizatorNaive
            <double, InstAmplComputerT, decltype(prediction_mode_inst_freq), 
                FilterT>, decltype(prediction_mode_inst_freq), 
                    FilterT>
                        (*naive_normalizer, prediction_mode_inst_freq, filter);*/
        }

        ~RecursiveFilterInstAmplChangesWithConstInstFreqDouble(){
            delete inst_ampl_normalizer;
        }

        template<Signal DataType, Signal OutType, Signal OldResultType>
        bool computeIter(const DataType& data, OutType& result_buffer, OldResultType* computer_buffer){
            iter_number++;
            true_iter_number++;
            using T = typename OutType::SampleType;

            if (true_iter_number > max_iter_number){
                for(int i = 0; i < data.size(); i++){
                    result_buffer[i] = result_cache[i];
                }
                return false;
            }

            if(result_cache.size()!=data.size()){
                result_cache.clear();
                for(auto i = 0; i < data.size(); i++){
                    result_cache.push_back(0.0);
                }
            }
            
            if(prediction_mode.size() != data.size()){
                prediction_mode.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_inst_freq.size() != data.size()){
                //prediction_mode_inst_freq.base->vec->clear();
                //for (auto i = 0; i < data.size(); i++){
                //    prediction_mode_inst_freq.base->vec->push_back(0.0);
                //}
                /*std::unreachable();*/
            }
            if(prediction_mode_ampl.size() != data.size()){
                prediction_mode_ampl.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_ampl.base->vec->push_back(0.0);
                }
            }
            //for (auto i = 0; i < data.size(); i++){
                //prediction_mode_inst_freq[i] = 0.0;
            //}
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_ampl[i] = 0.0;
            }
            for(auto i = 0; i < data.size(); i++){
                prediction_mode[i] = data[i] - result_buffer[i];
            }

            //for(int i = 0; i < data.size(); i++){
                //prediction_mode_inst_freq[i] = (*computer_buffer)[i];
            //}
            //auto label = 0;
            //IC(++label);
            inst_ampl_normalizer->inst_freq = &prediction_mode_inst_freq;
            inst_ampl_normalizer->compute(prediction_mode, *computer_buffer, prediction_mode_ampl);

            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = result_buffer[i] + (*computer_buffer)[i];
            }
            //IC(++label);
            filter->compute(*computer_buffer, prediction_mode, &prediction_mode_inst_freq);
            //IC(++label);
            error = UTILITY_MATH::signalsL2NormedDistance
                    <double, decltype(result_buffer), decltype(prediction_mode)>
                    (result_buffer, prediction_mode);
            //IC(error);

            if (error_old <= error){
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    for(int i = 0; i < data.size(); i++){
                        //(*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return false;
                }
                if(iter_number - good_iter_number > 9){
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    iter_number = good_iter_number;
                    for(int i = 0; i < data.size(); i++){
                        //(*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return computeIter(data, result_buffer, computer_buffer);
                }
            }
            else if(error < error_old){
                //IC(error_old, error, error_treshhold);
                error_old = error;
                good_iter_number = iter_number;
                
                for (int i = 0; i < data.size(); i++){
                    result_cache[i] = result_buffer[i];
                }
                if(static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        //(*computer_buffer)[i] = prediction_mode_inst_freq[i];
                    }
                    return false;
                }

            }

            for(int i = 0; i < data.size(); i++){
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                result_buffer[i] = (result_buffer[i] * coeff1 + prediction_mode[i] * coeff2) / (coeff1 + coeff2);
            }
            for(int i = 0; i < data.size(); i++){
                //(*computer_buffer)[i] = prediction_mode_inst_freq[i];
            }
            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqT>
        void compute(const DataType& data, OutType& result_buffer, InstFreqT* inst_freq){
            iter_number = 0;
            good_iter_number = 0;
            true_iter_number = 0;
            
            //this->inst_freq_computer->compute(data, *inst_freq, &result_buffer);
            this->filter->compute(data, result_buffer, inst_freq);

            new_result.base->vec->clear();
            for (int i = 0; i < data.size(); i++){
                new_result.base->vec->push_back(0.0);
            }
            prediction_mode_inst_freq.base->vec->clear();
            for (int i = 0; i < data.size(); i++){
                prediction_mode_inst_freq.base->vec->push_back((*inst_freq)[i]);
            }

            double der_avg = 0.0;
            auto sum = 0.0;
            for (int i = 1; i < data.size(); i++){
                der_avg += std::abs((data[i] - result_buffer[i])-(data[i-1] - result_buffer[i-1])) / data.size();
                sum += std::abs(data[i]);
            }
            error_old = sum * data.size();
            der_avg = der_avg / data.size();

            error_treshhold = der_avg * error_treshhold_muller;
            
            computeIter(data, result_buffer, &new_result);

            while(error > error_treshhold){
                if(!computeIter(data, result_buffer, &new_result)) break;
            }
        }   
    };



    
    template<typename U, Integrator<U> IntegratorT, Derivator<U> DerivatorT, Filter<U> FilterT, 
            PhaseComputer<U> PhaseComputerT, InstAmplComputer<U> InstAmplComputerT,
                InstFreqComputer<U> InstFreqComputerForModeT>
    struct RecursiveFilterInstAmplChangesDouble{
        using AdditionalDataType = SignalPrototype<U>;

        using BufferT = GenericSignal<SimpleVecWrapper<U>, true>;

        BufferT prediction_mode;
        GenericSignal<SimpleVecWrapper<std::pair<U, U>>, true> prediction_mode_inst_freq;
        BufferT prediction_mode_ampl;
        BufferT new_result;

        static BufferT & prediction_phase(){
            static BufferT * prediction_mode_phase = new BufferT;
            return *prediction_mode_phase;
        }

        //BufferT computer_buffer;

        std::vector<double> result_cache;
        
        constexpr static bool is_filter = true;

        double error;
        double error_treshhold = 0.1;
        double error_treshhold_muller = 0.01;

        size_t good_iter_number = 0;
        size_t true_iter_number = 0;
        size_t iter_number = 0;
        size_t max_iter_number = 1000;

        double error_old = 0.0;


        FilterT * filter;
        //InstFreqComputerT * inst_freq_computer;
        InstFreqComputerForModeT * inst_freq_computer_for_mode;
        PhaseComputerT * phase_computer;
        //PhaseComputerForModeT * phase_computer_for_mode;
        InstAmplComputerT * inst_ampl_computer;
        InstAmplNormalizator<double, IntegratorT, DerivatorT, InstAmplComputerT> 
            * inst_ampl_normalizer;

        RecursiveFilterInstAmplChangesDouble(IntegratorT & integrator, DerivatorT & derivator, 
            FilterT & filter,
                PhaseComputerT & phase_computer, InstAmplComputerT & inst_ampl_computer, 
                    InstFreqComputerForModeT & inst_freq_computer_for_mode){
            this->filter = &filter;
            this->inst_freq_computer_for_mode = &inst_freq_computer_for_mode;
            this->phase_computer = &phase_computer;
            inst_ampl_normalizer = new InstAmplNormalizator
                <double, IntegratorT, DerivatorT, InstAmplComputerT> (integrator, derivator, inst_ampl_computer);
        }

        ~RecursiveFilterInstAmplChangesDouble(){
            delete inst_ampl_normalizer;
        }

        template<Signal DataType, Signal OutType, Signal OldResultType>
        bool computeIter(const DataType& data, OutType& result_buffer, OldResultType* computer_buffer){
            iter_number++;
            true_iter_number++;

            if (true_iter_number > max_iter_number){
                for(int i = 0; i < data.size(); i++){
                    result_buffer[i] = result_cache[i];
                }
                return false;
            }

            using T = typename OutType::SampleType;

            if(result_cache.size()!=data.size()){
                result_cache.clear();
                for(auto i = 0; i < data.size(); i++){
                    result_cache.push_back(0.0);
                }
            }
            
            if(prediction_mode.size() != data.size()){
                prediction_mode.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode.base->vec->push_back(0.0);
                }
            }
            if(prediction_mode_inst_freq.size() != data.size()){
                prediction_mode_inst_freq.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_inst_freq.base->vec->push_back({0.0, 0.0});
                }
                ///*std::unreachable();*/
            }
            if(prediction_mode_ampl.size() != data.size()){
                prediction_mode_ampl.base->vec->clear();
                for (auto i = 0; i < data.size(); i++){
                    prediction_mode_ampl.base->vec->push_back(0.0);
                }
            }
            for (auto i = 0; i < data.size(); i++){
                prediction_mode_ampl[i] = 0.0;
            }
            for (auto i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = 0.0;
            }
            for(auto i = 0; i < data.size(); i++){
                prediction_mode[i] = data[i] - result_buffer[i];
            }

            if constexpr (InstFreqComputerForModeT::is_phase_based()){
                if (prediction_phase().size()!= data.size){
                    prediction_phase().base->vec->clear();
                    for (auto i = 0; i < data.size(); i++){
                        prediction_phase().base->vec->push_back(0.0);
                    }
                }
                for (auto i = 0; i < data.size(); i++){
                    prediction_phase()[i] = 0.0;
                }
                phase_computer->compute(prediction_mode, prediction_phase(), computer_buffer);
                inst_freq_computer_for_mode->compute(prediction_phase(), prediction_mode_inst_freq, computer_buffer);
            }
            else{
                inst_freq_computer_for_mode->compute(prediction_mode, prediction_mode_inst_freq, computer_buffer);
            }
            
            
            inst_ampl_normalizer->compute(prediction_mode, *computer_buffer, prediction_mode_ampl);

            for(int i = 0; i < data.size(); i++){
                (*computer_buffer)[i] = result_buffer[i] + (*computer_buffer)[i];
            }

            filter->compute(*computer_buffer, prediction_mode, &prediction_mode_inst_freq);
            error = UTILITY_MATH::signalsL2NormedDistance
                    <double, decltype(result_buffer), decltype(prediction_mode)>
                    (result_buffer, prediction_mode);

            if (error_old <= error){
                if (static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    return false;
                }
                if(iter_number - good_iter_number > 9){
                    for(int i = 0; i < data.size(); i++){
                        result_buffer[i] = result_cache[i];
                    }
                    iter_number = good_iter_number;
                    return computeIter(data, result_buffer, computer_buffer);
                }
            }
            else if(error < error_old){
                //IC(error_old, error, error_treshhold);
                error_old = error;
                good_iter_number = iter_number;
                
                for (int i = 0; i < data.size(); i++){
                    result_cache[i] = result_buffer[i];
                }
                if(static_cast<double>(true_iter_number) / static_cast<double>(iter_number)
                    > 10 && iter_number > 50) {
                    return false;
                }

            }

            for(int i = 0; i < data.size(); i += std::rand() % 3){
                double coeff1 = (std::rand() % 300 + 1) / 10.0;
                double coeff2 = (std::rand() % 300 + 1) / 10.0;
                result_buffer[i] = (result_buffer[i] * coeff1 + prediction_mode[i] * coeff2) / (coeff1 + coeff2);
            }
            return true;
        }

        template<Signal DataType, Signal OutType, Signal InstFreqT>
        void compute(const DataType& data, OutType& result_buffer, InstFreqT* inst_freq){
            iter_number = 0;
            good_iter_number = 0;
            true_iter_number = 0;
            
            //this->inst_freq_computer->compute(data, *inst_freq, &result_buffer);
            this->filter->compute(data, result_buffer, inst_freq);

            double der_avg = 0.0;
            auto sum = 0.0;
            for (int i = 1; i < data.size(); i++){
                der_avg += std::abs((data[i] - result_buffer[i])-(data[i-1] - result_buffer[i-1])) / data.size();
                sum += std::abs(data[i]);
            }
            new_result.base->vec->clear();
            for (int i = 0; i < data.size(); i++){
                new_result.base->vec->push_back(0.0);
            }
            prediction_mode_inst_freq.base->vec->clear();
            for (int i = 0; i < data.size(); i++){
                prediction_mode_inst_freq.base->vec->push_back((*inst_freq)[i]);
            }

            error_old = sum * data.size();
            der_avg = der_avg / data.size();

            error_treshhold = der_avg * error_treshhold_muller;
            
            computeIter(data, result_buffer, &new_result);

            while(error > error_treshhold){
                if(!computeIter(data, result_buffer, &new_result)) break;
            }
        }   
    };

    template<typename U, PhaseComputer<U> PhaseComputerT>
    struct SincResLocalFilterWithRes{
        constexpr static bool is_filter = true;
        GenericSignal<SimpleVecWrapper<U>, true> buffer1;
        GenericSignal<SimpleVecWrapper<U>, true> buffer2;
        using SignalT = decltype(buffer1);

        double freq;
        bool is_low_pass;
        double locality_coeff = 5.0;

        double period_muller = 1.0;

        PhaseComputerT * phase_computer;
        SincResLocalFilter<U> filter;

        template<Signal DataT, Signal OutT, Signal InstFreqT>
        void compute(const DataT & data, OutT & out, InstFreqT * inst_freq){
            std::vector<double> freq_conv;
            std::vector<double> freq_conv_image;
            for (int i = 0; i < data.size(); i++){
                freq_conv.push_back(1.0);
                freq_conv_image.push_back(1.0);
            }
            //std::cout << "created freq_conv and freq_conv_image" << std::endl;

            SignalT compute_buffer;
            for (auto i = 0; i < data.size(); i++){
                compute_buffer.base->vec->push_back(0.0);
            }

            compute_buffer.has_ovnership = true;
            //double base_inst_freq = INST_FREQ_COMPUTERS::instFreqNorm(data, out, *inst_freq, freq_conv, freq_conv_image);
            
            double base_inst_freq = 
                INST_FREQ_COMPUTERS::instFreqNormOnce(data, out, *inst_freq, freq_conv);

            //std::cout << "computed inst freq norm" << std::endl;

            phase_computer->compute(out, compute_buffer, nullptr);

            base_inst_freq = 1.0 /
                (static_cast<double>(data.size()) /
                (compute_buffer[data.size() - 1] / 
                2.0 / std::numbers::pi)) / period_muller;

            filter.freq = base_inst_freq;
            filter.is_low_pass = is_low_pass;
            filter.locality_coeff = locality_coeff;
            /*std::cout << base_inst_freq << " "  
                << is_low_pass << " "
                << locality_coeff << " "
                << period_muller << std::endl;*/
            filter.compute(out, compute_buffer, nullptr);

            /*std::cout << "resampled_filtering" << std::endl;
            matplot::plot(*out.base->vec);
            matplot::hold(true);
            matplot::plot(*compute_buffer.base->vec);
            matplot::hold(false);
            matplot::show();*/


            INST_FREQ_COMPUTERS::backInstFreqNormOnce(compute_buffer, out, freq_conv);

            /*std::cout << "after resampled_filtering" << std::endl;
            matplot::plot(*data.base->vec);
            matplot::hold(true);
            matplot::plot(*out.base->vec);
            matplot::hold(false);
            matplot::show();*/
        }
    };


    template<typename U, PhaseComputer<U> PhaseComputerT, InstFreqComputer<U> InstFreqComputerT>
    struct SincResLocalFilterWithResReq{
        constexpr static bool is_filter = true;
        GenericSignal<SimpleVecWrapper<U>, true> buffer1;
        GenericSignal<SimpleVecWrapper<U>, true> buffer2;
        using SignalT = decltype(buffer1);

        double locality_coeff = 5.0;

        double period_muller = 1.0;

        PhaseComputerT * phase_computer;
        InstFreqComputerT * inst_freq_computer;

        SincResLocalFilterWithRes<U, PhaseComputerT> filter;

        size_t max_iters = 3;

        bool debug = true;

        template<Signal DataT, Signal OutT, Signal ComputeBufferT>
        void compute(const DataT & data, OutT & out, ComputeBufferT * compute_buffer){
            bool flag = true;
            SignalT low_freq_mode;
            SignalT low_freq_mode_phase;
            SignalT high_freq_mode;
            SignalT compute_buffer2;

            for (auto i = 0; i < data.size(); i++){
                low_freq_mode.base->vec->push_back(data[i]);
                low_freq_mode_phase.base->vec->push_back(data[i]);;
                high_freq_mode.base->vec->push_back(data[i]);
                compute_buffer2.base->vec->push_back(data[i]);
            }

            low_freq_mode.has_ovnership = true;
            low_freq_mode_phase.has_ovnership = true;
            high_freq_mode.has_ovnership = true;
            compute_buffer2.has_ovnership = true;
            
            if constexpr (inst_freq_computer->is_phase_based()){
                phase_computer->compute(data, out, compute_buffer);
                if (debug){
                    //std::cout << "phase" << std::endl;
                    //out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                inst_freq_computer->compute(out, *compute_buffer, &compute_buffer2);
                if (debug){
                    //std::cout << "inst_freq" << std::endl;
                    //compute_buffer->show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
            }
            else{
                inst_freq_computer->compute(data, *compute_buffer, &compute_buffer2);
                if (debug){
                    //std::cout << "inst_freq" << std::endl;
                    //compute_buffer->show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
            }
            
            filter.phase_computer = phase_computer;
            filter.locality_coeff = locality_coeff;
            filter.period_muller = period_muller;
            filter.is_low_pass = true;
            filter.compute(data, out, compute_buffer);
            if (debug){
                std::cout << "filtered_data" << std::endl;
                matplot::plot(*data.base->vec);
                matplot::hold(true);
                matplot::plot(*out.base->vec);
                matplot::hold(false);
                matplot::show();
                //out.show(NP_DSP::ONE_D::PlottingKind::Simple);
            }

            for (auto i = 0; i < data.size(); i++){
                high_freq_mode[i] -= out[i];
            } 
            if (debug){
                std::cout << "first iter high_freq_mode of signal" << std::endl;
                high_freq_mode.show(NP_DSP::ONE_D::PlottingKind::Simple);
            }
            size_t iter_number = 1;
            while (flag){
                if (iter_number >= max_iters){
                    flag = false;
                    break;
                }
                iter_number++;
                if constexpr (inst_freq_computer->is_phase_based()){
                    phase_computer->compute(high_freq_mode, compute_buffer2, compute_buffer);
                    if (debug){
                        std::cout << "phase of high_freq_mode" << std::endl;
                        compute_buffer2.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    }
                    inst_freq_computer->compute(compute_buffer2, *compute_buffer, &low_freq_mode);
                    if (debug){
                        std::cout << "inst_freq of high_freq_mode" << std::endl;
                        compute_buffer->show(NP_DSP::ONE_D::PlottingKind::Simple);
                    }
                }
                else{
                    inst_freq_computer->compute(high_freq_mode, *compute_buffer, &low_freq_mode);
                    if (debug){
                        std::cout << "inst_freq of high_freq_mode" << std::endl;
                        compute_buffer->show(NP_DSP::ONE_D::PlottingKind::Simple);
                    }
                }
                filter.compute(high_freq_mode, low_freq_mode, compute_buffer);
                if (debug){
                    std::cout << "low_freq_mode of high_freq_mode" << std::endl;
                    //low_freq_mode.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    if (debug){
                        std::cout << "filtered_data" << std::endl;
                        matplot::plot(*high_freq_mode.base->vec);
                        matplot::hold(true);
                        matplot::plot(*low_freq_mode.base->vec);
                        matplot::hold(false);
                        matplot::show();
                        //out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    }
                }
                for (int i = 0; i < data.size(); i++){
                    high_freq_mode[i] -= low_freq_mode[i];
                }
                if (debug){
                    std::cout << "high_freq_mode of high_freq_mode" << std::endl;
                    high_freq_mode.show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                phase_computer->compute(low_freq_mode, low_freq_mode_phase, compute_buffer);
                if (debug){
                    std::cout << "phase of low_freq_mode of high_freq_mode" << std::endl;
                    low_freq_mode_phase.show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                for (int i = 0; i < data.size(); i++){
                    out[i] = out[i] + low_freq_mode[i];
                }
                if (debug){
                    std::cout << "low_freq_mode_of Signal" << std::endl;
                    out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                if (low_freq_mode_phase[data.size() - 1] < 6.28) {
                    flag = false;
                }
            }
        }
    };

    template<typename U, PhaseComputer<U> PhaseComputerT>
    struct SincResLocalFilterWithResReqV2{
        constexpr static bool is_filter = true;
        GenericSignal<SimpleVecWrapper<U>, true> buffer1;
        GenericSignal<SimpleVecWrapper<U>, true> buffer2;
        using SignalT = decltype(buffer1);

        double locality_coeff = 5.0;

        double period_muller = 1.0;

        PhaseComputerT * phase_computer;

        SincResLocalFilterWithRes<U, PhaseComputerT> filter;

        size_t max_iters = 3;

        bool debug = true;

        template<Signal DataT, Signal OutT, Signal ComputeBufferT>
        void compute(const DataT & data, OutT & out, ComputeBufferT * compute_buffer){
            bool flag = true;
            SignalT low_freq_mode;
            SignalT low_freq_mode_phase;
            SignalT high_freq_mode;
            SignalT compute_buffer2;

            for (auto i = 0; i < data.size(); i++){
                low_freq_mode.base->vec->push_back(data[i]);
                low_freq_mode_phase.base->vec->push_back(data[i]);;
                high_freq_mode.base->vec->push_back(data[i]);
                compute_buffer2.base->vec->push_back(data[i]);
            }

            low_freq_mode.has_ovnership = true;
            low_freq_mode_phase.has_ovnership = true;
            high_freq_mode.has_ovnership = true;
            compute_buffer2.has_ovnership = true;
            
            phase_computer->compute(data, out, compute_buffer);
            if (debug){
                //std::cout << "phase" << std::endl;
                //out.show(NP_DSP::ONE_D::PlottingKind::Simple);
            }
            for (size_t i = 0; i < out.size(); i++){
                (*compute_buffer)[i] = phase_computer->derive(i) / std::numbers::pi / 2.0;
            }
            //inst_freq_computer->compute(out, *compute_buffer, &compute_buffer2);
            if (debug){
                //std::cout << "inst_freq" << std::endl;
                //compute_buffer->show(NP_DSP::ONE_D::PlottingKind::Simple);
            }
            
            
            filter.phase_computer = phase_computer;
            filter.locality_coeff = locality_coeff;
            filter.period_muller = period_muller;
            filter.is_low_pass = true;
            filter.compute(data, out, compute_buffer);
            if (debug){
                std::cout << "filtered_data" << std::endl;
                matplot::plot(*data.base->vec);
                matplot::hold(true);
                matplot::plot(*out.base->vec);
                matplot::hold(false);
                matplot::show();
                //out.show(NP_DSP::ONE_D::PlottingKind::Simple);
            }

            for (auto i = 0; i < data.size(); i++){
                high_freq_mode[i] -= out[i];
            } 
            if (debug){
                std::cout << "first iter high_freq_mode of signal" << std::endl;
                high_freq_mode.show(NP_DSP::ONE_D::PlottingKind::Simple);
            }
            size_t iter_number = 1;
            while (flag){
                if (iter_number >= max_iters){
                    flag = false;
                    break;
                }
                iter_number++;

                phase_computer->compute(high_freq_mode, compute_buffer2, compute_buffer);
                if (debug){
                    //std::cout << "phase" << std::endl;
                    //out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                for (size_t i = 0; i < out.size(); i++){
                    (*compute_buffer)[i] = phase_computer->derive(i) / std::numbers::pi / 2.0;
                }
                //inst_freq_computer->compute(out, *compute_buffer, &compute_buffer2);
                if (debug){
                    //std::cout << "inst_freq" << std::endl;
                    //compute_buffer->show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                
                filter.compute(high_freq_mode, low_freq_mode, compute_buffer);
                if (debug){
                    std::cout << "low_freq_mode of high_freq_mode" << std::endl;
                    //low_freq_mode.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    if (debug){
                        std::cout << "filtered_data" << std::endl;
                        matplot::plot(*high_freq_mode.base->vec);
                        matplot::hold(true);
                        matplot::plot(*low_freq_mode.base->vec);
                        matplot::hold(false);
                        matplot::show();
                        //out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    }
                }
                for (int i = 0; i < data.size(); i++){
                    high_freq_mode[i] -= low_freq_mode[i];
                }
                if (debug){
                    std::cout << "high_freq_mode of high_freq_mode" << std::endl;
                    high_freq_mode.show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                phase_computer->compute(low_freq_mode, low_freq_mode_phase, compute_buffer);
                if (debug){
                    std::cout << "phase of low_freq_mode of high_freq_mode" << std::endl;
                    low_freq_mode_phase.show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                for (int i = 0; i < data.size(); i++){
                    out[i] = out[i] + low_freq_mode[i];
                }
                if (debug){
                    std::cout << "low_freq_mode_of Signal" << std::endl;
                    out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                }
                if (low_freq_mode_phase[data.size() - 1] < 6.28) {
                    flag = false;
                }
            }
        }
    };
}
