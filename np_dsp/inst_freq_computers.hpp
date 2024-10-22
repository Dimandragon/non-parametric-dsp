#pragma once

#include <cstddef>
#include <phase_computers.hpp>
#include <npdsp_concepts.hpp>
#include <signals.hpp>
#include <cmath>
#include <numbers>
#include <vector>
#include <utility>
#include <utility_math.hpp>
#include <approximators.hpp>
#include <npdsp_config.hpp>
#include <complex>
#include <math.h>

namespace NP_DSP::ONE_D::INST_FREQ_COMPUTERS {
    
    using InstFreqDerivativeBasedKind = PHASE_COMPUTERS::InstFreqDerivativeBasedKind;

    
    template<typename U, Integrator<U> IntegratorT,
        Derivator<U> DerivatorT, InstFreqDerivativeBasedKind kind>
    struct DerivativeBased {
        using AdditionalDataType = SignalPrototype<U>;

        constexpr static bool is_inst_freq_computer = true;
        constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

        using IntegratorType = IntegratorT;
        using DerivatorType = DerivatorT;


        IntegratorType integrator;
        DerivatorType derivator;

        double variability = 1.0;

        DerivativeBased(IntegratorT integrator_o,
                        DerivatorT derivator_o) {
            integrator = integrator_o;
            derivator = derivator_o;
        }

        constexpr static bool is_phase_based() {
            return false;
        }

        template <Signal DataType, Signal OutType, Signal ComputerNufferType>
        void compute(const DataType& data, OutType& out, ComputerNufferType* computer_buffer) {
            using T = typename OutType::SampleType;
            if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental) {
                derivator.compute(data, *computer_buffer, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    (*computer_buffer)[i] = std::atan((*computer_buffer)[i]);
                }
                derivator.compute(*computer_buffer, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    out[i] = std::abs(out[i]) / (std::numbers::pi * 2.0);
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage) {
                derivator.compute(data, *computer_buffer, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    (*computer_buffer)[i] = std::atan((*computer_buffer)[i]);
                }
                derivator.compute(*computer_buffer, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    out[i] = std::abs(out[i]);
                }
                integrator.compute(out, *computer_buffer, nullptr);

                for (auto i = 0; i < data.size(); i++) {
                    auto approx_answer = static_cast<T>(0.0);
                    auto old_approx_answer = approx_answer;
                    auto counter = 0;
                    while (approx_answer < 2.0 * std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer = approx_answer;
                        approx_answer = computer_buffer->interpolate(i + counter, SignalKind::Monotone) 
                            - computer_buffer->interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss = std::numbers::pi * 2.0 * variability - old_approx_answer;
                    auto right_loss = approx_answer - std::numbers::pi * 2.0 * variability;
                    auto sum_loss = left_loss + right_loss;
                    auto period = (static_cast<T>(counter) - right_loss / sum_loss) * 2;
                    out[i] = static_cast<T>(1.0 / period) * variability;
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage) {
                derivator.compute(data, *computer_buffer, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    (*computer_buffer)[i] = std::atan((*computer_buffer)[i]);
                }
                derivator.compute(*computer_buffer, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    out[i] = std::abs(out[i]);
                }
                integrator.compute(out, *computer_buffer, nullptr);
                for (auto i = 0; i < data.size(); i++) {
                    auto approx_answer_right = static_cast<T>(0.0);
                    auto old_approx_answer_right = approx_answer_right;
                    auto counter = 0;
                    while (approx_answer_right < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_right = approx_answer_right;
                        approx_answer_right = computer_buffer->interpolate(i + counter, SignalKind::Monotone) - computer_buffer->
                                              interpolate(i, SignalKind::Monotone);
                    }
                    auto left_loss_right_edge = std::numbers::pi * variability - old_approx_answer_right;
                    auto right_loss_right_edge = approx_answer_right - std::numbers::pi * variability;
                    auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                    auto right_edge = static_cast<T>(counter) - right_loss_right_edge / sum_loss_right_edge;

                    auto approx_answer_left = static_cast<T>(0.0);
                    auto old_approx_answer_left = approx_answer_left;
                    counter = 0;
                    while (approx_answer_left < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_left = approx_answer_left;
                        approx_answer_left = computer_buffer->interpolate(i, SignalKind::Monotone) - computer_buffer->
                                             interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss_left_edge = std::numbers::pi * variability - old_approx_answer_left;
                    auto right_loss_left_edge = approx_answer_left - std::numbers::pi * variability;
                    auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                    auto left_edge = static_cast<T>(counter) - right_loss_left_edge / sum_loss_left_edge;

                    auto period = right_edge + left_edge;
                    out[i] = static_cast<T>(1.0 / period) * variability;
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveDouble) {
                //todo
            }
        }
    };

    template<typename PhaseT>
    double findPeriodDistancDouble(const PhaseT & phase, double idx, double pad){
        return std::abs(phase.interpolate(idx, SignalKind::Monotone) - 
            phase.interpolate(idx + pad, SignalKind::Monotone));
    }

    template<typename PhaseT>
    double findPeriodDistanceAverage(const PhaseT & phase, double idx, double pad){
        return phase.interpolate(idx + pad, SignalKind::Monotone) - 
            phase.interpolate(idx - pad, SignalKind::Monotone);
    }

    
    template<typename U, Integrator<U> IntegratorT,
        Derivator<U> DerivatorT, InstFreqDerivativeBasedKind kind>
    struct ComputedOnPhase {
        using AdditionalDataType = GENERAL::Nil;

        constexpr static bool is_inst_freq_computer = true;
        constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

        //using OptFunction = OptFn;
        //OptFunction opt_function = opt_fn;

        using IntegratorType = IntegratorT;
        using DerivatorType = DerivatorT;


        IntegratorType integrator;
        DerivatorType derivator;
        double variability = 1.0;

        constexpr static bool is_phase_based() {
            return true;
        }

        //static_assert(OutType::is_writable == true);
        ComputedOnPhase(){}

        ComputedOnPhase(IntegratorT integrator_o,
                        DerivatorT derivator_o) {
            integrator = integrator_o;
            derivator = derivator_o;
        }

        template<Signal DataType, Signal OutType>
        void compute(const DataType& phase, OutType& out, auto * nil) {
            compute(phase, out, nullptr);
        }

        template<Signal DataType, Signal OutType>
        void compute(const DataType& phase, OutType& out, std::nullptr_t nil) {
            using T = typename OutType::SampleType;
            //nil may bee nullptr
            auto delta = 0.01;
            if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental) {
                derivator.compute(phase, out, nullptr);
                for (int i = 0; i < phase.size(); i++) {
                    out[i] = out[i] / (std::numbers::pi * 2.0);
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage) {
                for (auto i = 0; i < phase.size(); i++) {
                    #if 0
                    auto approx_answer_x = 0.0;
                    auto approx_answer_y = 0.0;
                    auto counter = 0;
                    double x_left = 0.0;
                    double x_right = phase.size();
                    double y_left = 0.0;
                    double y_right = findPeriodDistanceAverage(phase, i, x_right);

                    while (x_right - x_left > 1.0) {
                        approx_answer_x = UTILITY_MATH::backLinearInterpolate<double, double>
                            ({x_left, y_left}, {x_right, y_right}, 
                                2.0 * std::numbers::pi * variability);
                        approx_answer_y = findPeriodDistanceAverage(phase, i, approx_answer_x);

                        if (approx_answer_y > y_right){
                            y_left = y_right;
                            y_right = approx_answer_y;
                            x_left = x_right;
                            x_right = approx_answer_x;
                        }
                        else if (approx_answer_y < y_left){
                            y_right = y_left;
                            x_right = x_left;
                            y_left = approx_answer_y;
                            x_left = approx_answer_x;
                        }
                        else if (approx_answer_y > 2.0 * std::numbers::pi * variability){
                            y_right = approx_answer_y;
                            x_right = approx_answer_x;
                        }
                        else if (approx_answer_y < 2.0 * std::numbers::pi * variability){
                            y_left = approx_answer_y;
                            x_left = approx_answer_x;
                        }
                        if (std::abs(approx_answer_y - 2.0 * std::numbers::pi * variability) < delta / approx_answer_x){
                            /*x_left = static_cast<int64_t>(approx_answer_x);
                            x_right = x_left + 1;
                            y_left = findPeriodDistanceAverage(phase, i, x_left);
                            y_right = findPeriodDistanceAverage(phase, i, x_right);*/
                            break;
                        }
                        //IC(i, approx_answer_y, approx_answer_x, x_left, x_right, y_left, y_right, 2.0 * std::numbers::pi * variability);
                    }
                    #endif
                    auto val_expr = [&](double idx){
                        return std::abs(phase.interpolate(i + idx, SignalKind::Monotone) - phase.interpolate(i - idx, SignalKind::Monotone));
                    };
                    auto size_expr = [&](){
                        return phase.size();
                    };
                    ExpressionWrapper<double, double, decltype(val_expr), GENERAL::Nil,
                        decltype(size_expr), false> expr_wrapper(val_expr, size_expr);
                    GenericSignal<decltype(expr_wrapper), false> phase_(expr_wrapper);
                    double approx_answer_x = phase_.findMonotone
                        (2.0 * std::numbers::pi * variability, {}, {}, {}, 0.01);

                    auto period = approx_answer_x * 2.0;
                    out[i] = static_cast<T>(1.0 / period) * variability;
                    //IC(out[i], i, period, approx_answer_x);
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage) {
                for (auto i = 0; i < phase.size(); i++) {
                    auto approx_answer_right = static_cast<T>(0.0);
                    auto old_approx_answer_right = approx_answer_right;
                    auto counter = 0;
                    while (approx_answer_right < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_right = approx_answer_right;
                        approx_answer_right = phase.interpolate(i + counter, SignalKind::Monotone) - phase.interpolate(i, SignalKind::Monotone);
                    }
                    auto left_loss_right_edge = std::numbers::pi * variability - old_approx_answer_right;
                    auto right_loss_right_edge = approx_answer_right - std::numbers::pi * variability;
                    auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                    auto right_edge = static_cast<T>(counter) - right_loss_right_edge / sum_loss_right_edge;

                    auto approx_answer_left = static_cast<T>(0.0);
                    auto old_approx_answer_left = approx_answer_left;
                    counter = 0;
                    while (approx_answer_left < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_left = approx_answer_left;
                        approx_answer_left = phase.interpolate(i, SignalKind::Monotone) - phase.interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss_left_edge = std::numbers::pi * variability - old_approx_answer_left;
                    auto right_loss_left_edge = approx_answer_left - std::numbers::pi * variability;
                    auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                    auto left_edge = static_cast<T>(counter) - right_loss_left_edge / sum_loss_left_edge;

                    auto period = right_edge + left_edge;
                    out[i] = static_cast<T>(1.0 / period) * variability;
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveDouble) {
                for (auto i = 0; i < phase.size(); i++) {
                    auto approx_answer_right = 0.0;
                    auto old_approx_answer_right = approx_answer_right;
                    auto counter = 0;
                    while (approx_answer_right < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_right = approx_answer_right;
                        approx_answer_right = phase.interpolate(i + counter, SignalKind::Monotone) - phase.interpolate(i, SignalKind::Monotone);
                    }
                    auto left_loss_right_edge = std::numbers::pi * variability - old_approx_answer_right;
                    auto right_loss_right_edge = approx_answer_right - std::numbers::pi * variability;
                    auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                    auto right_edge = static_cast<double>(counter) - right_loss_right_edge / sum_loss_right_edge;

                    auto approx_answer_left = 0.0;
                    auto old_approx_answer_left = approx_answer_left;
                    counter = 0;
                    while (approx_answer_left < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_left = approx_answer_left;
                        approx_answer_left = phase.interpolate(i, SignalKind::Monotone) - phase.interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss_left_edge = std::numbers::pi * variability - old_approx_answer_left;
                    auto right_loss_left_edge = approx_answer_left - std::numbers::pi * variability;
                    auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                    auto left_edge = static_cast<double>(counter) - right_loss_left_edge / sum_loss_left_edge;

                    out[i].first = 0.5 / left_edge * variability;
                    out[i].second = 0.5 / right_edge * variability;
                }
            }
        }
    };

    
    template<typename U, Integrator<U> IntegratorT,
        Derivator<U> DerivatorT, InstFreqDerivativeBasedKind kind, PhaseComputer<U> PhaseComputerT>
    struct PhaseBased {
        using AdditionalDataType = SignalPrototype<U>;

        constexpr static bool is_inst_freq_computer = true;
        constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

        //using OptFunction = OptFn;
        //OptFunction opt_function = opt_fn;

        using IntegratorType = IntegratorT;
        using DerivatorType = DerivatorT;
        using PhaseComputerType = PhaseComputerT;

        PhaseComputerType* phase_computer;
        IntegratorType integrator;
        DerivatorType derivator;

        constexpr static bool is_phase_based() {
            return false;
        }

        double variability = 1.0;

        //static_assert(OutType::is_writable == true);

        PhaseBased(IntegratorT integrator_o,
                   DerivatorT derivator_o, PhaseComputerT& phase_computer_o) {
            integrator = integrator_o;
            derivator = derivator_o;
            phase_computer = &phase_computer_o;
        }

        template<Signal DataType, Signal OutType, Signal ComputerBufferType>
        void compute(const DataType& data, OutType& out, ComputerBufferType* computer_buffer) {
            using T = typename OutType::SampleType;
            //auto phase_computer =
            //        ONE_D::PHASE_COMPUTERS::ArctgScaledToExtremums<DataType, AdditionalDataType, OutType,
            //            IntegratorType, DerivatorType>(integrator, derivator);
            if constexpr (std::convertible_to<typename PhaseComputerT::AdditionalDataType, GENERAL::Nil>) {
                phase_computer->compute(data, *computer_buffer, nullptr);
            } else {
                phase_computer->compute(data, *computer_buffer, &out);
            }

            if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental) {
                derivator.compute(*computer_buffer, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    out[i] = out[i] / (std::numbers::pi * 2.0);
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage) {
                for (auto i = 0; i < data.size(); i++) {
                    auto approx_answer = static_cast<T>(0.0);
                    auto old_approx_answer = approx_answer;
                    auto counter = 0;
                    while (approx_answer < 2.0 * std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer = approx_answer;
                        approx_answer = computer_buffer->interpolate(i + counter, SignalKind::Universal) - 
                            computer_buffer->interpolate(i - counter, SignalKind::Universal);
                    }
                    auto left_loss = std::numbers::pi * 2.0 * variability - old_approx_answer;
                    auto right_loss = approx_answer - std::numbers::pi * 2.0 * variability;
                    auto sum_loss = left_loss + right_loss;
                    auto period = (static_cast<T>(counter) - right_loss / sum_loss) * 2;
                    out[i] = static_cast<T>(1.0 / period) * variability;
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage) {
                for (auto i = 0; i < data.size(); i++) {
                    auto approx_answer_right = static_cast<T>(0.0);
                    auto old_approx_answer_right = approx_answer_right;
                    auto counter = 0;
                    while (approx_answer_right < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_right = approx_answer_right;
                        approx_answer_right = computer_buffer->interpolate(i + counter, SignalKind::Monotone) - computer_buffer->
                                              interpolate(i, SignalKind::Monotone);
                    }
                    auto left_loss_right_edge = std::numbers::pi * variability - old_approx_answer_right;
                    auto right_loss_right_edge = approx_answer_right - std::numbers::pi * variability;
                    auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                    auto right_edge = static_cast<T>(counter) - right_loss_right_edge / sum_loss_right_edge;

                    auto approx_answer_left = static_cast<T>(0.0);
                    auto old_approx_answer_left = approx_answer_left;
                    counter = 0;
                    while (approx_answer_left < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_left = approx_answer_left;
                        approx_answer_left = computer_buffer->interpolate(i, SignalKind::Monotone) - computer_buffer->
                                             interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss_left_edge = std::numbers::pi * variability - old_approx_answer_left;
                    auto right_loss_left_edge = approx_answer_left - std::numbers::pi * variability;
                    auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                    auto left_edge = static_cast<T>(counter) - right_loss_left_edge / sum_loss_left_edge;

                    auto period = right_edge + left_edge;
                    out[i] = static_cast<T>(1.0 / period) * variability;
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveDouble) {
                for (auto i = 0; i < data.size(); i++) {
                    auto approx_answer_right = 0.0;
                    auto old_approx_answer_right = approx_answer_right;
                    auto counter = 0;
                    while (approx_answer_right < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_right = approx_answer_right;
                        approx_answer_right = computer_buffer->interpolate(i + counter, SignalKind::Monotone) - computer_buffer->
                                              interpolate(i, SignalKind::Monotone);
                    }
                    auto left_loss_right_edge = std::numbers::pi * variability - old_approx_answer_right;
                    auto right_loss_right_edge = approx_answer_right - std::numbers::pi * variability;
                    auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                    auto right_edge = static_cast<double>(counter) - right_loss_right_edge / sum_loss_right_edge;

                    auto approx_answer_left = static_cast<double>(0.0);
                    auto old_approx_answer_left = approx_answer_left;
                    counter = 0;
                    while (approx_answer_left < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_left = approx_answer_left;
                        approx_answer_left = computer_buffer->interpolate(i, SignalKind::Monotone) - computer_buffer->
                                             interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss_left_edge = std::numbers::pi * variability - old_approx_answer_left;
                    auto right_loss_left_edge = approx_answer_left - std::numbers::pi * variability;
                    auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                    auto left_edge = static_cast<double>(counter) - right_loss_left_edge / sum_loss_left_edge;

                    out[i].first = static_cast<double>(0.5 / left_edge) * variability;
                    out[i].second = static_cast<double>(0.5 / right_edge) * variability;
                }
            }
        }
    };

    
    template<typename U, Integrator<U> IntegratorT,
        Derivator<U> DerivatorT, InstFreqDerivativeBasedKind kind>
    struct DerivativeBasedWithExternalOptParametr {
        using AdditionalDataType = SignalPrototype<U>;

        constexpr static bool is_inst_freq_computer = true;
        constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

        using IntegratorType = IntegratorT;
        using DerivatorType = DerivatorT;


        IntegratorType integrator;
        DerivatorType derivator;

        double variability = 1.0;

        DerivativeBasedWithExternalOptParametr(IntegratorT integrator_o,
                                               DerivatorT derivator_o) {
            integrator = integrator_o;
            derivator = derivator_o;
        }

        constexpr static bool is_phase_based() {
            return false;
        }

        template<Signal DataType, Signal OutType, Signal ComputerBufferType>
        void compute(const DataType& data, OutType& out, ComputerBufferType* computer_buffer) {
            using T = typename OutType::SampleType;
            if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental) {
                derivator.compute(data, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    out[i] = std::atan(out[i]) * (*computer_buffer)[i];
                }
                derivator.compute(out, *computer_buffer, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    out[i] = std::abs((*computer_buffer[i])) / (std::numbers::pi * 2.0);
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage) {
                derivator.compute(data, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    (*computer_buffer)[i] = std::atan(out[i]) * (*computer_buffer)[i];
                }
                derivator.compute(*computer_buffer, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    out[i] = std::abs(out[i]);
                }
                integrator.compute(out, *computer_buffer, nullptr);

                for (auto i = 0; i < data.size(); i++) {
                    auto approx_answer = static_cast<T>(0.0);
                    auto old_approx_answer = approx_answer;
                    auto counter = 0;
                    while (approx_answer < 2.0 * std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer = approx_answer;
                        approx_answer = computer_buffer->interpolate(i + counter, SignalKind::Monotone) - computer_buffer->interpolate(
                                            i - counter, SignalKind::Monotone);
                    }
                    auto left_loss = std::numbers::pi * 2.0 * variability - old_approx_answer;
                    auto right_loss = approx_answer - std::numbers::pi * 2.0 * variability;
                    auto sum_loss = left_loss + right_loss;
                    auto period = (static_cast<T>(counter) - right_loss / sum_loss) * 2;
                    out[i] = static_cast<T>(1.0 / period) * variability;
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage) {
                derivator.compute(data, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    (*computer_buffer)[i] = std::atan(out[i]) * (*computer_buffer)[i];
                }
                derivator.compute(*computer_buffer, out, nullptr);
                for (int i = 0; i < data.size(); i++) {
                    out[i] = std::abs(out[i]);
                }
                integrator.compute(out, *computer_buffer, nullptr);
                if constexpr (NP_DSP::CONFIG::debug) {
                    //computer_buffer.show(PlottingKind::Interpolate);
                }
                for (auto i = 0; i < data.size(); i++) {
                    auto approx_answer_right = static_cast<T>(0.0);
                    auto old_approx_answer_right = approx_answer_right;
                    auto counter = 0;
                    while (approx_answer_right < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_right = approx_answer_right;
                        approx_answer_right = computer_buffer->interpolate(i + counter, SignalKind::Monotone) - computer_buffer->
                                              interpolate(i, SignalKind::Monotone);
                    }
                    auto left_loss_right_edge = std::numbers::pi * variability - old_approx_answer_right;
                    auto right_loss_right_edge = approx_answer_right - std::numbers::pi * variability;
                    auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                    auto right_edge = static_cast<T>(counter) - right_loss_right_edge / sum_loss_right_edge;

                    auto approx_answer_left = static_cast<T>(0.0);
                    auto old_approx_answer_left = approx_answer_left;
                    counter = 0;
                    while (approx_answer_left < std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer_left = approx_answer_left;
                        approx_answer_left = computer_buffer->interpolate(i, SignalKind::Monotone) - computer_buffer->
                                             interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss_left_edge = std::numbers::pi * variability - old_approx_answer_left;
                    auto right_loss_left_edge = approx_answer_left - std::numbers::pi * variability;
                    auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                    auto left_edge = static_cast<T>(counter) - right_loss_left_edge / sum_loss_left_edge;

                    auto period = right_edge + left_edge;
                    out[i] = static_cast<T>(1.0 / period) * variability;
                }
            }
        }
    };

     enum class ExtremumsBasedComputeInstFreqKind { Simple, Linear };

    
    template<ExtremumsBasedComputeInstFreqKind compute_kind>
    struct ExtremumsBased {
        using AdditionalDataType = GENERAL::Nil;

        constexpr static ExtremumsBasedComputeInstFreqKind kind = compute_kind;
        constexpr static bool is_inst_freq_computer = true;

        constexpr static bool is_phase_based() {
            return false;
        }

        template<Signal DataType, Signal OutType>
        void compute(const DataType& data, OutType& out, auto * nil) {
            using T = typename OutType::SampleType;
            std::vector<int> extremums;
            extremums.push_back(0);
            for (int i = 1; i < data.size() - 1; i++) {
                if ((data[i] >= data[i - 1] &&
                     data[i] > data[i + 1]) ||
                    (data[i] > data[i - 1] &&
                     data[i] >= data[i + 1]) ||
                    (data[i] <= data[i - 1] &&
                     data[i] < data[i + 1]) ||
                    (data[i] < data[i - 1] &&
                     data[i] <= data[i + 1])) {
                    extremums.push_back(i);
                }
            }
            extremums.push_back(static_cast<int>(data.size() - 1));

            if (extremums.size() > 2) {
                extremums[0] = extremums[1] * 2 - extremums[2];
                auto last = extremums.size() - 1;
                extremums[last] = extremums[last - 1] * 2 - extremums[last - 2];
                if (extremums[0] > 0) {
                    extremums[0] = -extremums[1];
                }
                if (extremums[last] < extremums.size() - 1) {
                    extremums[last] = data.size() + data.size() - extremums[last - 1];
                }
            } else {
                auto last = extremums.size() - 1;
                extremums[0] = -extremums[1];
                extremums[last] = data.size() + data.size() - extremums[last - 1];
            }


            if constexpr (kind == ExtremumsBasedComputeInstFreqKind::Simple) {
                auto left = extremums[0];
                auto right = extremums[1];
                auto counter = 0;
                for (auto i = 0; i < out.size(); i++) {
                    if (i > right) {
                        if (counter + 2 < extremums.size()) {
                            counter++;
                        }
                        left = extremums[counter];
                        right = extremums[counter + 1];
                    }
                    out[i] = static_cast<T>(0.5 / (right - left));
                }
            } else if constexpr (kind == ExtremumsBasedComputeInstFreqKind::Linear) {
                std::vector<std::pair<T, T>> points;
                points.push_back({static_cast<T>(0), static_cast<T>(0.5 / (extremums[1] - extremums[0]))});
                for (auto i = 0; i < extremums.size() - 1; i++) {
                    points.push_back({
                        static_cast<T>((extremums[i + 1] + extremums[i]) / 2.0),
                        static_cast<T>(0.5 / (extremums[i + 1] - extremums[i]))
                    });
                }
                points.push_back({
                    static_cast<T>(data.size() - 1), static_cast<T>
                    (0.5 / (extremums[extremums.size() - 1] - extremums[extremums.size() - 2]))
                });

                auto left = points[0];
                auto right = points[1];
                auto counter = 0;

                for (auto i = 0; i < out.size(); i++) {
                    if (i > right.first) {
                        counter++;
                        left = points[counter];
                        right = points[counter + 1];
                    }
                    out[i] = UTILITY_MATH::linearInterpolate(left, right, static_cast<T>(i));
                }
            } else {
                /*std::unreachable();*/
            }
        }

        template<Signal DataType, Signal OutType>
        void compute(const DataType& data, OutType& out, std::nullptr_t nil) {
            using T = typename OutType::SampleType;
            std::vector<int> extremums;
            extremums.push_back(0);
            for (int i = 1; i < data.size() - 1; i++) {
                if ((data[i] >= data[i - 1] &&
                     data[i] > data[i + 1]) ||
                    (data[i] > data[i - 1] &&
                     data[i] >= data[i + 1]) ||
                    (data[i] <= data[i - 1] &&
                     data[i] < data[i + 1]) ||
                    (data[i] < data[i - 1] &&
                     data[i] <= data[i + 1])) {
                    extremums.push_back(i);
                }
            }
            extremums.push_back(static_cast<int>(data.size() - 1));

            if (extremums.size() > 2) {
                extremums[0] = extremums[1] * 2 - extremums[2];
                auto last = extremums.size() - 1;
                extremums[last] = extremums[last - 1] * 2 - extremums[last - 2];
                if (extremums[0] > 0) {
                    extremums[0] = -extremums[1];
                }
                if (extremums[last] < extremums.size() - 1) {
                    extremums[last] = data.size() + data.size() - extremums[last - 1];
                }
            } else {
                auto last = extremums.size() - 1;
                extremums[0] = -extremums[1];
                extremums[last] = data.size() + data.size() - extremums[last - 1];
            }


            if constexpr (kind == ExtremumsBasedComputeInstFreqKind::Simple) {
                auto left = extremums[0];
                auto right = extremums[1];
                auto counter = 0;
                for (auto i = 0; i < out.size(); i++) {
                    if (i > right) {
                        if (counter + 2 < extremums.size()) {
                            counter++;
                        }
                        left = extremums[counter];
                        right = extremums[counter + 1];
                    }
                    out[i] = static_cast<T>(0.5 / (right - left));
                }
            } else if constexpr (kind == ExtremumsBasedComputeInstFreqKind::Linear) {
                std::vector<std::pair<T, T>> points;
                points.push_back({static_cast<T>(0), static_cast<T>(0.5 / (extremums[1] - extremums[0]))});
                for (auto i = 0; i < extremums.size() - 1; i++) {
                    points.push_back({
                        static_cast<T>((extremums[i + 1] + extremums[i]) / 2.0),
                        static_cast<T>(0.5 / (extremums[i + 1] - extremums[i]))
                    });
                }
                points.push_back({
                    static_cast<T>(data.size() - 1), static_cast<T>
                    (0.5 / (extremums[extremums.size() - 1] - extremums[extremums.size() - 2]))
                });

                auto left = points[0];
                auto right = points[1];
                auto counter = 0;

                for (auto i = 0; i < out.size(); i++) {
                    if (i > right.first) {
                        counter++;
                        left = points[counter];
                        right = points[counter + 1];
                    }
                    out[i] = UTILITY_MATH::linearInterpolate(left, right, static_cast<T>(i));
                }
            } else {
                /*std::unreachable();*/
            }
        }
    };

    
    template<typename U, Integrator<U> IntegratorT,
        Derivator<U> DerivatorT, InstFreqDerivativeBasedKind kind>
    struct PeriodAndExtremumsBased {
        using AdditionalDataType = GENERAL::Nil;

        constexpr static bool is_inst_freq_computer = true;
        constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

        using IntegratorType = IntegratorT;
        using DerivatorType = DerivatorT;

        IntegratorType integrator;
        DerivatorType derivator;

        double approx_order_coeff = 1.0;
        int tile_size = 128;

        double variability = 0.0;

        double max_error = 1000000.0;

        PeriodAndExtremumsBased(IntegratorType integrator_in, DerivatorType derivator_in) {
            integrator = integrator_in;
            derivator = derivator_in;
        }

        constexpr static bool is_phase_based() {
            return false;
        }

        template<Signal DataType, Signal OutType>
        void compute(const DataType& data, OutType& out, auto * nil) {
            // compute extremums inst freq ->
            // compute period based with external opt parameter equal const 1
            // use fourier based approximator
            // train it to approximate external opt parameter to minimize loss with extremums
            using T = typename OutType::SampleType;
            using SampleType = T;

            std::vector<T> extremums_freq_vec(data.size());
            SimpleVecWrapper<T> extremums_freq_wrapper(extremums_freq_vec);
            GenericSignal<SimpleVecWrapper<T>, true> extremums_freq(extremums_freq_wrapper);

            ExtremumsBased<ExtremumsBasedComputeInstFreqKind::Linear> extremums_based;
            extremums_based.compute(data, extremums_freq, nullptr);

            std::vector<T> external_opt_parametr_vector(data.size());
            SimpleVecWrapper<T> external_opt_parametr_wrapper(external_opt_parametr_vector);
            GenericSignal<T, true> external_opt_parametr(external_opt_parametr_wrapper);

            DerivativeBasedWithExternalOptParametr<T, IntegratorType, DerivatorType, kind>
                    inst_freq_computer(integrator, derivator);

            inst_freq_computer.variability = variability;

            for (auto i = 0; i < external_opt_parametr.size(); i++) {
                external_opt_parametr[i] = static_cast<SampleType>(1.);
            }
            inst_freq_computer.compute(data, out, &external_opt_parametr);

            auto loss = [&](auto& approximator) {
                approximator.is_actual = false;
                std::vector<double> add_error = {};
                for (auto i = 0; i < data.size(); i++) {
                    external_opt_parametr[i] = approximator.compute(i);
                    if (external_opt_parametr[i] < 0.05) {
                        add_error.push_back(0.05 - external_opt_parametr[i]);
                        external_opt_parametr[i] = 0.05;
                    } else {
                        add_error.push_back(0.0);
                    }
                }
                double accum = 0.0;

                inst_freq_computer.compute(data, out, &external_opt_parametr);

                for (auto i = 0; i < extremums_freq.size(); i++) {
                    accum += add_error[i] + (extremums_freq[i] - out[i]) * (extremums_freq[i] - out[i]) / 1000;
                    /// data[i];
                }
                return accum;
            };

            auto bySampleError = [&](auto& approximator, auto i) {
                double add_error = 0.0;
                external_opt_parametr[i] = approximator.approximated_data[i].real();
                if (external_opt_parametr[i] < 0.05) {
                    add_error += 0.05 - external_opt_parametr[i];
                    external_opt_parametr[i] = 0.05;
                }
                inst_freq_computer.compute(data, out, &external_opt_parametr);
                return add_error + (extremums_freq[i] - out[i]) * (extremums_freq[i] - out[i]) / approximator.tile_size;
            };

            auto stopPoint = [](auto losses_different, auto& approximator) {
                if (losses_different > 0.00001) {
                    return false;
                } else {
                    return true;
                }
            };
            auto approximator = APPROX::FourierSeriesBased<decltype(loss), decltype(stopPoint),
                        APPROX::FSApproxKind::Positive, decltype(bySampleError)>
                    (loss, external_opt_parametr, stopPoint);
            approximator.tile_size = tile_size;
            approximator.bySampleLoss = &bySampleError;
            approximator.is_actual = false;
            approximator.setApproxOrderRatio(approx_order_coeff);
            approximator.max_value = 10;
            approximator.train();
        }

        template<Signal DataType, Signal OutType>
        void compute(const DataType& data, OutType& out, std::nullptr_t * nil) {
            // compute extremums inst freq ->
            // compute period based with external opt parameter equal const 1
            // use fourier based approximator
            // train it to approximate external opt parameter to minimize loss with extremums
            using T = typename OutType::SampleType;
            using SampleType = T;

            std::vector<T> extremums_freq_vec(data.size());
            SimpleVecWrapper<T> extremums_freq_wrapper(extremums_freq_vec);
            GenericSignal<SimpleVecWrapper<T>, true> extremums_freq(extremums_freq_wrapper);

            ExtremumsBased<ExtremumsBasedComputeInstFreqKind::Linear> extremums_based;
            extremums_based.compute(data, extremums_freq, nullptr);

            std::vector<T> external_opt_parametr_vector(data.size());
            SimpleVecWrapper<T> external_opt_parametr_wrapper(external_opt_parametr_vector);
            GenericSignal<T, true> external_opt_parametr(external_opt_parametr_wrapper);

            DerivativeBasedWithExternalOptParametr<T, IntegratorType, DerivatorType, kind>
                    inst_freq_computer(integrator, derivator);

            inst_freq_computer.variability = variability;

            for (auto i = 0; i < external_opt_parametr.size(); i++) {
                external_opt_parametr[i] = static_cast<SampleType>(1.);
            }
            inst_freq_computer.compute(data, out, &external_opt_parametr);

            auto loss = [&](auto& approximator) {
                approximator.is_actual = false;
                std::vector<double> add_error = {};
                for (auto i = 0; i < data.size(); i++) {
                    external_opt_parametr[i] = approximator.compute(i);
                    if (external_opt_parametr[i] < 0.05) {
                        add_error.push_back(0.05 - external_opt_parametr[i]);
                        external_opt_parametr[i] = 0.05;
                    } else {
                        add_error.push_back(0.0);
                    }
                }
                double accum = 0.0;

                inst_freq_computer.compute(data, out, &external_opt_parametr);

                for (auto i = 0; i < extremums_freq.size(); i++) {
                    accum += add_error[i] + (extremums_freq[i] - out[i]) * (extremums_freq[i] - out[i]) / 1000;
                    /// data[i];
                }
                return accum;
            };

            auto bySampleError = [&](auto& approximator, auto i) {
                double add_error = 0.0;
                external_opt_parametr[i] = approximator.approximated_data[i].real();
                if (external_opt_parametr[i] < 0.05) {
                    add_error += 0.05 - external_opt_parametr[i];
                    external_opt_parametr[i] = 0.05;
                }
                inst_freq_computer.compute(data, out, &external_opt_parametr);
                return add_error + (extremums_freq[i] - out[i]) * (extremums_freq[i] - out[i]) / approximator.tile_size;
            };

            auto stopPoint = [](auto losses_different, auto& approximator) {
                if (losses_different > 0.00001) {
                    return false;
                } else {
                    return true;
                }
            };
            auto approximator = APPROX::FourierSeriesBased<decltype(loss), decltype(stopPoint),
                        APPROX::FSApproxKind::Positive, decltype(bySampleError)>
                    (loss, external_opt_parametr, stopPoint);
            approximator.tile_size = tile_size;
            approximator.bySampleLoss = &bySampleError;
            approximator.is_actual = false;
            approximator.setApproxOrderRatio(approx_order_coeff);
            approximator.max_value = 10;
            approximator.train();
        }
    };

    
    template<typename U, Integrator<U> IntegratorT,
        Derivator<U> DerivatorT, InstFreqDerivativeBasedKind kind>
    struct PeriodAndExtremumsBasedExternal {
        using AdditionalDataType = SignalPrototype<U>;

        constexpr static bool is_inst_freq_computer = true;
        constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

        using IntegratorType = IntegratorT;
        using DerivatorType = DerivatorT;


        IntegratorType integrator;
        DerivatorType derivator;

        double approx_order_coeff = 1.0;
        int tile_size = 128;
        double variability = 0.0;

        double max_error = 10000000.;

        PeriodAndExtremumsBasedExternal(IntegratorType integrator_in, DerivatorType derivator_in) {
            integrator = integrator_in;
            derivator = derivator_in;
        }

        constexpr static bool is_phase_based() {
            return false;
        }

        template<Signal DataType, Signal OutType, Signal ComputerBufferType>
        void compute(const DataType& data, OutType& out, ComputerBufferType* computer_buffer) {
            // compute extremums inst freq ->
            // compute period based with external opt parameter equal const 1
            // use fourier based approximator
            // train it to approximate external opt parameter to minimize loss with extremums
            using T = typename OutType::SampleType;

            GenericSignal<SimpleVecWrapper<T>, true> extremums_freq;

            for (int i = 0; i < data.size(); i++) {
                static_cast<SimpleVecWrapper<double> *>(extremums_freq.base)->vec->push_back(0.);
            }

            ExtremumsBased<ExtremumsBasedComputeInstFreqKind::Linear> extremums_based;
            extremums_based.compute(data, extremums_freq, nullptr);

            GenericSignal<SimpleVecWrapper<T>, true> external_opt_parametr;
            for (int i = 0; i < data.size(); i++) {
                static_cast<SimpleVecWrapper<double> *>(external_opt_parametr.base)->vec->push_back(1.);
            }

            DerivativeBased<T, IntegratorType, DerivatorType, kind>
                    inst_freq_computer(integrator, derivator);

            inst_freq_computer.variability = variability;

            inst_freq_computer.compute(data, out, computer_buffer);

            auto loss = [&](auto& approximator) {
                approximator.is_actual = false;
                for (auto i = 0; i < data.size(); i++) {
                    external_opt_parametr[i] = approximator.template compute<double, size_t>(i);
                }
                double accum = 0.0;

                for (auto i = 0; i < extremums_freq.size(); i++) {
                    if (out[i] != 0.0) {
                        accum += (extremums_freq[i] - out[i] * external_opt_parametr[i])
                                * (extremums_freq[i] - out[i] * external_opt_parametr[i]) / 1000 / out[i];
                    }
                }
                return accum;
            };

            auto stopPoint = [](auto losses_different, auto& approximator) {
                if (losses_different > 0.000001) {
                    //todo move precision to external parameter
                    return false;
                } else {
                    return true;
                }
            };

            auto bySampleError = [&](auto& approximator, auto i) {
                external_opt_parametr[i] = approximator.approximated_data[i].real();
                return (extremums_freq[i] - out[i] * external_opt_parametr[i])
                       * (extremums_freq[i] - out[i] * external_opt_parametr[i]) / approximator.tile_size / out[i];
            };
            auto approximator = APPROX::FourierSeriesBased<decltype(loss), decltype(stopPoint),
                        APPROX::FSApproxKind::Simple, decltype(bySampleError)>
                    (loss, external_opt_parametr, stopPoint);
            approximator.is_actual = false;

            approximator.setApproxOrderRatio(approx_order_coeff);
            approximator.max_value = 10;
            approximator.tile_size = tile_size;
            approximator.bySampleLoss = &bySampleError;
            approximator.train();
            for (auto i = 0; i < data.size(); i++) {
                out[i] = out[i] * approximator.template compute<double, size_t>(i);
            }
        }
    };

    

    template<Signal DataT, Signal OutT, Signal InstFreqT>
    double instFreqNorm(const DataT & data, OutT & out, const InstFreqT & inst_freq, 
        std::vector<double> & freq_conv, std::vector<double> & freq_conv_image){
        freq_conv_image.clear();
        auto delta = 0.5;
        double freq_avg = 0.0;
        for (int i = 0; i < data.size(); i++){
            freq_avg += inst_freq[i];
        }
        freq_avg = freq_avg / data.size();
        if(freq_avg == 0){
            return 0.;
        }

        double iter_predict_int = 0.0;
        double iter_predict = 0.0;
        double temp = 0.0;
        double temp1;

        std::optional<std::pair<double, double>> max;
        std::optional<std::pair<double, double>> min;

        const double lim = static_cast<double>(out.size()) + 0.1 * delta - 1.0;
        int iter = 0;
        

        auto findIterLambda = [&]() ->bool{
            temp = 0.0;
            for (int i = 0; i < out.size(); i++){
                double i_f = inst_freq.interpolate(temp, SignalKind::Universal);
                double temp1;
                if (i_f == 0.){
                    temp1 = freq_avg;
                }
                else{
                    temp1 = freq_avg / i_f;
                }
                temp += std::abs(temp1);
            }
            
            if (temp >= (double)out.size() + delta - 1.){
                if (!max){
                    max = {freq_avg, temp};
                }
                else if (max){
                    if (freq_avg < max->first){
                        max = {freq_avg, temp};
                    }
                }
            }
            else if (temp <= (double)out.size() - 1.0){
                if (!min){
                    min = {freq_avg, temp};
                }
                else{
                    if (freq_avg > min->first){
                        min = {freq_avg, temp};
                    }
                }
            }
            else{
                return false;
            }
            if (!max && !min){
                freq_avg = std::sqrt(freq_avg*freq_avg * lim * lim / temp / temp);
            }
            else if (max && !min){
                freq_avg = std::sqrt(freq_avg*freq_avg * (lim - delta) * (lim - delta) / max->second / max->second);
            }
            else if (min && !max){
                freq_avg = freq_avg * (lim + lim) / min->second;
                if (freq_avg == 0.0){
                    int i = *reinterpret_cast<int *>(0);
                }
            }
            else{
                if (max->first == min->first){
                    freq_avg = max->first;
                    return false;
                }
                if (iter % 2){
                    freq_avg = UTILITY_MATH::backLinearInterpolate<double, double>
                        (*min, *max, lim);
                }
                else{
                    freq_avg = (min->first / 2.0 + max->first / 2.0);
                }
                iter++;
            }
            return true;
        };
        while (findIterLambda()){}

        double temp_res = temp;
        temp = 0.0;
        auto counter = 0;
        while (temp<=temp_res){
            double i_f = inst_freq.interpolate(temp, SignalKind::Universal);
            double temp1;
            if (i_f == 0.){
                temp1 = 1;
            }
            else{
                temp1 = freq_avg / i_f;
            }
            //out_signal.push(linear_interpolate(signal_in, temp));
            if(counter == out.size()){
                break;
            }
            out[counter] = data.interpolate(temp, SignalKind::Universal);
            //freq_conv_image.push(temp1 * linear_interpolate(freq_conv, temp));
            if (static_cast<int>(temp) >= freq_conv.size() - 1){
                auto idx_ = freq_conv.size() - 1;
                freq_conv_image.push_back(temp1 * UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[idx_-1]},
                    {static_cast<int>(temp) + 1, freq_conv[idx_]}, temp));
            }
            else{
                freq_conv_image.push_back(temp1 * UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp + 1)]}, temp));
            }
            temp = temp + temp1;
            counter++;
        }

        std::swap(freq_conv, freq_conv_image);
        freq_conv_image.clear();

        double freq_conv_red = 0.0;
        for (int i = 0; i < freq_conv.size(); i++){
            freq_conv_red += freq_conv[i];
        }

        return freq_avg;
    }

    template<Signal DataT, Signal OutT, Signal InstFreqT>
    double instFreqNormOnce(const DataT & data, OutT & out, const InstFreqT & inst_freq, 
        std::vector<double> & freq_conv)
    {   
        APPROX::ModifiedAkimaBasedWithNoTrain<DataT> approximator;
        approximator.loadData(data);
        double d_avg = 0.0;
        freq_conv.clear();

        for (int i = 0; i < data.size() - 1; i++){
            d_avg += 1.0 / (inst_freq[i] + inst_freq[i + 1]) * 2.0;
        }

        d_avg = d_avg / (data.size() - 1);

        if (d_avg == 0.0){
            return 0.0;
        }

        double iter_predict_int = 0.0;
        double iter_predict = 0.0;
        double temp = 0.0;

        for (size_t i = 0; i < data.size(); i++){
            freq_conv.push_back(temp);
            //out[i] = data.interpolate(temp, SignalKind::Universal);
            out[i] = approximator.compute(temp);
            
            if (i == data.size() - 1){
                break;
            }
            temp += 1.0 / (inst_freq[i] + inst_freq[i + 1]) * 2.0 / d_avg;
        }
        //IC(temp);
        return d_avg;
        //std::optional<std::pair<>>
    }

    template<Signal DataT, Signal OutT>
    void backInstFreqNormOnce(DataT const & data, OutT & out, std::vector<double> & freq_conv){
        size_t counter = 0;
        APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> approximator;
        std::vector<double> y_data(data.size());
        std::vector<double> x_data(data.size());
        for (int i = 0; i < data.size(); i++){
            y_data[i] = data[i];
            x_data[i] = freq_conv[i];
        }
        approximator.loadData(x_data, y_data);

        for (size_t i = 0; i < data.size(); i++){
            //out[i] = data.interpolate(freq_conv[i], SignalKind::Universal);
            //IC(freq_conv[i], i);
            while (freq_conv[counter] < i){
                counter++;
            }
            double idx = 0.0;
            if (counter != 0){
                idx = UTILITY_MATH::linearInterpolate<double, double>
                    ({freq_conv[counter-1], double(counter-1)}, {freq_conv[counter], double(counter-1)}, i);
            }
            //out[i] = data.interpolate(idx, SignalKind::Universal);
            out[i] = approximator.compute(i);
        }
    }

    template<Signal DataT, Signal OutT>
    double instFreqNormExtrBased(DataT const & data, OutT & out,
        APPROX::PiecewiseCubicHermitePolynomialBasedWithNoTrain<std::vector<double>> & idx_approx, 
        std::vector<double> & extremums)
    {
        //std::vector<double> extremums;
        UTILITY_MATH::computeExtremums<DataT, double>
            (data, extremums, UTILITY_MATH::ExtremumsKind::Simple);
        
        double period = static_cast<double>(data.size() - 1) / static_cast<double>(extremums.size() - 1);

        std::vector<double> extremums_old_idx;
        std::vector<double> extremums_new_idx;
        extremums_old_idx.push_back(-extremums[1]);
        extremums_new_idx.push_back(-period);
        for (int i = 0; i < extremums.size(); i++){
            extremums_old_idx.push_back(extremums[i]);
            extremums_new_idx.push_back(i * period);
        }
        extremums_old_idx.push_back(extremums[extremums.size() - 1] * 2 - extremums[extremums.size() - 2]);
        extremums_new_idx.push_back(extremums.size() * period);

        idx_approx.loadData(extremums_new_idx, extremums_old_idx);

        APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> signal_approx;
        //signal_approx.loadData(*(data.base->vec));
        std::vector<double> data_vec;
        for (int i = 0; i < data.size(); i++){
            data_vec.push_back(data[i]);
        }
        signal_approx.loadData(data_vec);

        for (int i = 0; i < data.size(); i++){
            out[i] = signal_approx.compute(idx_approx.compute(i));
            //IC(i, idx_approx.compute(i));
        }

        return period;
    }

    template<Signal DataT, Signal OutT>
    double instFreqNormComputedOnExtremums(DataT const & data, OutT & out,
        APPROX::PiecewiseCubicHermitePolynomialBasedWithNoTrain<std::vector<double>> & idx_approx, 
        const std::vector<double> & extremums)
    {
        //std::vector<double> extremums;
        //UTILITY_MATH::computeExtremums<DataT, double>
        //    (data, extremums, UTILITY_MATH::ExtremumsKind::Simple);
        
        double period = static_cast<double>(data.size() - 1) / static_cast<double>(extremums.size() - 1);

        std::vector<double> extremums_old_idx;
        std::vector<double> extremums_new_idx;
        extremums_old_idx.push_back(-extremums[1]);
        extremums_new_idx.push_back(-period);
        for (int i = 0; i < extremums.size(); i++){
            extremums_old_idx.push_back(extremums[i]);
            extremums_new_idx.push_back(i * period);
        }
        extremums_old_idx.push_back(extremums[extremums.size() - 1] * 2 - extremums[extremums.size() - 2]);
        extremums_new_idx.push_back(extremums.size() * period);

        idx_approx.loadData(extremums_new_idx, extremums_old_idx);

        APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> signal_approx;
        //signal_approx.loadData(*(data.base->vec));
        std::vector<double> data_vec;
        for (int i = 0; i < data.size(); i++){
            data_vec.push_back(data[i]);
        }
        signal_approx.loadData(data_vec);

        for (int i = 0; i < data.size(); i++){
            out[i] = signal_approx.compute(idx_approx.compute(i));
            //IC(i, idx_approx.compute(i));
        }

        return period;
    }

    template<Signal DataT, Signal OutT>
    void backInstFreqNormExtrBased(DataT const & data, OutT & out,
        APPROX::PiecewiseCubicHermitePolynomialBasedWithNoTrain<std::vector<double>> & idx_approx){
        APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>> signal_res_approx;
        std::vector<double> data_vec;

        auto val_expr = [&](double idx){
            return idx_approx.compute(idx);
        };
        auto size_expr = [&](){
            return data.size();
        };
        ExpressionWrapper<double, double, decltype(val_expr), GENERAL::Nil, decltype(size_expr), false> expr_wrapper(val_expr, size_expr);
        GenericSignalRExpr<decltype(expr_wrapper)> resampled_idx_signal(expr_wrapper);

        for (int i = 0; i < data.size(); i++){
            data_vec.push_back(data[i]);
        }
        signal_res_approx.loadData(data_vec);

        for(int i = 0; i < data.size(); i++){
            //IC(i, resampled_idx_signal.findMonotone(i, {}, {}, {i}, 0.01));
            out[i] = signal_res_approx.compute(resampled_idx_signal.findMonotone(i, {}, {}, {}, 0.01));
        }
    }
    
    template<Signal DataT, Signal OutT, Signal InstFreqT>
    double instFreqNormDouble(const DataT & data, const OutT & out, const InstFreqT & inst_freq, 
        std::vector<double> & freq_conv, std::vector<double> & freq_conv_image){
        
        auto val_expr = [&](size_t idx) {
            return ((*inst_freq)[idx].first + (*inst_freq)[idx].second)/ 2.0;;
        };
        auto size_expr = [&]() {
            return data.size();
        };
        
        ExpressionWrapper<double, size_t,
            decltype(val_expr), GENERAL::Nil, decltype(size_expr), false>
                inst_freq_expr(val_expr, size_expr);

        const GenericSignal<decltype(inst_freq_expr), false> inst_freq_e(inst_freq_expr);
        
        freq_conv_image.clear();
        double freq_avg = 0.0;
        for (int i = 0; i < data.size(); i++){
            freq_avg += inst_freq_e[i];
        }
        freq_avg = freq_avg / data.size();
        
        

        size_t iter_predict = 0;
        double temp = 0.0;
        double temp1;
        while (iter_predict != out.size()){
            freq_avg = std::sqrt(freq_avg*freq_avg / data.size() / data.size() * iter_predict * iter_predict);
            iter_predict = 0;
            double temp = 0.0;
            while (temp <= (static_cast<double>(data.size()) - 1.0)){
                temp1 = freq_avg / inst_freq_e.interpolate(temp, SignalKind::Universal);//linear_interpolate(freq_arr, temp);
                temp = temp + temp1;
                iter_predict++;
            }
        }

        temp = 0.0;
        auto counter = 0;
        while (temp<=data.size()-1.0){
            temp1 = freq_avg / inst_freq_e.interpolate(temp, SignalKind::Universal);
            //out_signal.push(linear_interpolate(signal_in, temp));
            if(counter == out.size()){
                break;
            }
            out[counter] = data.interpolate(temp, SignalKind::Universal);
            //freq_conv_image.push(temp1 * linear_interpolate(freq_conv, temp));
            freq_conv_image.push_back(UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp + 1)]}, temp));
            temp = temp + temp1;
            counter++;
        }

        std::swap(freq_conv, freq_conv_image);
        freq_conv_image.clear();
        return freq_avg;
    }

    
    template<Signal DataT, Signal OutT>
    void backInstFreqNorm(DataT const & data, OutT & out, std::vector<double> & freq_conv){
        double temp = 0.0;
        double temp1 = 0.0;
        double sum = 0.0;

        for (int i = 0; i < freq_conv.size(); i++){
            sum += freq_conv[i];
        }
        
        size_t iter_predict = 0;

        auto first_find_iter_loop_body = [&](int flag){
            double fc = UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp) + 1]}, temp);
            
            if (0.0 == fc)
            {
                temp = temp + 1.0;
                iter_predict++;
            }
            else{
                temp = temp + 1.0 / fc;
                iter_predict++;
            }
        };

        while (temp <= (data.size() - 1)){
            first_find_iter_loop_body(1);
        }

        double one = 1.0;
        one = one / data.size() * iter_predict;
        temp = 0.0;
        size_t counter = 0;
        auto find_iter = [&](bool arg){
            double fc = UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp) + 1]}, temp);
            if (0.0 == fc)
            {
                temp = temp + 1.0;
                counter++;
            }
            else{
                //arrout.push(linear_interpolate::<T>(arrin, temp));
                if (counter < data.size()){
                    out[counter] = data.interpolate(temp, SignalKind::Universal);
                    //IC(counter, temp, out[counter]);
                    counter++;
                    temp = temp + one / fc;
                }
                else{
                    return false;
                }
            }
            return true;
        };
        while (temp <= static_cast<double>(data.size()) - 1.0){
            if (!find_iter(true)){
                break;
            }
        }
        //IC(sum, iter_predict, temp);
    }

    template<Signal DataT, Signal OutT>
    void backInstFreqNormNew(DataT const & data, OutT & out, std::vector<double> & freq_conv){
        //todo
        double temp = 0.0;
        double temp1 = 0.0;

        double one_ = 1.0;
        auto delta = 0.5;

        const double lim = static_cast<double>(out.size()) + 0.1 * delta - 1.0;
        
        std::optional<std::pair<double, double>> max;
        std::optional<std::pair<double, double>> min;

        size_t iter;

        auto find_iter_loop_body_ext = [&](){
            temp = 0.0;
            for (int i = 0; i < data.size(); i++){
                double id;
                if (temp >= data.size() - 1){
                    id = data.size() - 1;
                }
                else{
                    id = static_cast<int>(temp + 1);
                }
                double fc = UTILITY_MATH::linearInterpolate<double, double>
                    ({id - 1, freq_conv[id - 1]},
                    {id, freq_conv[id]}, temp);

                if (0.0 < delta / 1000.)
                {
                    temp = temp + one_;
                }
                else{
                    temp = temp + one_ / fc;
                }
            }

            if (temp >= (double)data.size() + delta - 1.0){
                if (!max){
                    max = {one_, temp};
                }
                else if (max){
                    if (one_ < max->first){
                        max = {one_, temp};
                    }
                }
            }
            else if (temp <= (double)out.size() - 1.0){
                if (!min){
                    min = {one_, temp};
                }
                else{
                    if (one_ > min->first){
                        min = {one_, temp};
                    }
                }
            }
            else{
                return false;
            }
            if (!max && !min){
                one_ = std::sqrt(one_*one_ * lim * lim / temp / temp);
            }
            else if (max && !min){
                one_ = std::sqrt(one_*one_ * (lim - delta) * (lim - delta) / max->second / max->second);
            }
            else if (min && !max){
                one_ = one_ * (lim + lim) / min->second;
            }
            else{
                if (max->first == min->first){
                    return false;
                }
                if (iter % 2){
                    one_ = UTILITY_MATH::backLinearInterpolate<double, double>
                        (*min, *max, lim);
                }
                else{
                    one_ = (min->first / 2.0 + max->first / 2.0);
                }
                iter++;
            }
            return true;
        };

        while(find_iter_loop_body_ext()){}

        double temp_res = temp;
        temp = 0.0;
        size_t counter = 0;
        while(temp < temp_res){
            double fc = UTILITY_MATH::linearInterpolate<double, double>
                ({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp) + 1]}, temp);
            if (fc == 0.){
                temp1 = one_;
            }
            else{
                temp1 = one_ / fc;
            }
            if(counter == out.size()){
                break;
            }
            out[counter] = data.interpolate(temp, SignalKind::Universal);
            temp += temp1;
            counter++;
        }

        IC(counter, temp);
    }
}
