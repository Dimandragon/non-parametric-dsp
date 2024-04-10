#pragma once

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

        ComputedOnPhase(IntegratorT integrator_o,
                        DerivatorT derivator_o) {
            integrator = integrator_o;
            derivator = derivator_o;
        }

        template<Signal DataType, Signal OutType>
        void compute(const DataType& phase, OutType& out, auto * nil) {
            using T = typename OutType::SampleType;
            //nil may bee nullptr
            if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental) {
                derivator.compute(phase, out, nullptr);
                for (int i = 0; i < phase.size(); i++) {
                    out[i] = out[i] / (std::numbers::pi * 2.0);
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage) {
                for (auto i = 0; i < phase.size(); i++) {
                    auto approx_answer = static_cast<T>(0.0);
                    auto old_approx_answer = approx_answer;
                    auto counter = 0;
                    while (approx_answer < 2.0 * std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer = approx_answer;
                        approx_answer = phase.interpolate(i + counter, SignalKind::Monotone) - 
                            phase.interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss = std::numbers::pi * 2.0 * variability - old_approx_answer;
                    auto right_loss = approx_answer - std::numbers::pi * 2.0 * variability;
                    auto sum_loss = left_loss + right_loss;
                    auto period = (static_cast<T>(counter) - right_loss / sum_loss) * 2;
                    out[i] = static_cast<T>(1.0 / period) * variability;
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

        template<Signal DataType, Signal OutType>
        void compute(const DataType& phase, OutType& out, std::nullptr_t nil) {
            using T = typename OutType::SampleType;
            //nil may bee nullptr
            if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental) {
                derivator.compute(phase, out, nullptr);
                for (int i = 0; i < phase.size(); i++) {
                    out[i] = out[i] / (std::numbers::pi * 2.0);
                }
            } else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage) {
                for (auto i = 0; i < phase.size(); i++) {
                    auto approx_answer = static_cast<T>(0.0);
                    auto old_approx_answer = approx_answer;
                    auto counter = 0;
                    while (approx_answer < 2.0 * std::numbers::pi * variability) {
                        counter++;
                        old_approx_answer = approx_answer;
                        approx_answer = phase.interpolate(i + counter, SignalKind::Monotone) -
                            phase.interpolate(i - counter, SignalKind::Monotone);
                    }
                    auto left_loss = std::numbers::pi * 2.0 * variability - old_approx_answer;
                    auto right_loss = approx_answer - std::numbers::pi * 2.0 * variability;
                    auto sum_loss = left_loss + right_loss;
                    auto period = (static_cast<T>(counter) - right_loss / sum_loss) * 2;
                    out[i] = static_cast<T>(1.0 / period) * variability;
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
                std::unreachable();
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
                std::unreachable();
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
    double InstFreqNorm(const DataT & data, OutT & out, const InstFreqT & inst_freq, 
        std::vector<double> & freq_conv, std::vector<double> & freq_conv_image){
        freq_conv_image.clear();
        double freq_avg = 0.0;
        for (int i = 0; i < data.size(); i++){
            freq_avg += inst_freq[i];
        }
        freq_avg = freq_avg / data.size();

        size_t iter_predict = 0;
        double temp = 0.0;
        double temp1;
        while (iter_predict != out.size()){
            if (iter_predict != 0){
                freq_avg = std::sqrt(freq_avg*freq_avg / data.size() / data.size() * iter_predict * iter_predict);
            }
            iter_predict = 0;
            double temp = 0.0;
            while (temp <= (static_cast<double>(data.size()) - 1.0)){
                temp1 = freq_avg / inst_freq.interpolate(temp, SignalKind::Universal);//linear_interpolate(freq_arr, temp);
                temp = temp + temp1;
                iter_predict++;
            }
        }

        temp = 0.0;
        auto counter = 0;
        while (temp<=(data.size()-1.0)){
            temp1 = freq_avg / inst_freq.interpolate(temp, SignalKind::Universal);
            //out_signal.push(linear_interpolate(signal_in, temp));
            if(counter == out.size()){
                break;
            }
            out[counter] = data.interpolate(temp, SignalKind::Universal);
            //freq_conv_image.push(temp1 * linear_interpolate(freq_conv, temp));
            freq_conv_image.push_back(temp1 * UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp + 1)]}, temp));
            temp = temp + temp1;
            counter++;
        }

        std::swap(freq_conv, freq_conv_image);
        freq_conv_image.clear();
        return freq_avg;
    }

    
    template<Signal DataT, Signal OutT, Signal InstFreqT>
    double InstFreqNormDouble(const DataT & data, const OutT & out, const InstFreqT & inst_freq, 
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

        while (temp <= (data.size() - 1)){
            temp1 = 1.0 / UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp) + 1]}, temp);
            
            if (0.0 == UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp) + 1]}, temp))
            {
                temp = temp + 1.0;
            }
            else{
                temp = temp + temp1;
                iter_predict = iter_predict+1;
            }
        }

        double one = 1.0;
        one = one / data.size() * iter_predict;
        temp = 0.0;
        size_t counter = 0;
        while (temp <= static_cast<double>(data.size()) - 1.0){
            temp1 = one / UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp) + 1]}, temp);
            if (0.0 == UTILITY_MATH::linearInterpolate<double, double>({static_cast<int>(temp), freq_conv[static_cast<int>(temp)]},
                {static_cast<int>(temp) + 1, freq_conv[static_cast<int>(temp) + 1]}, temp))
            {
                temp = temp + 1.0;
            }
            else{
                //arrout.push(linear_interpolate::<T>(arrin, temp));
                if (counter < data.size()){
                    out[counter] = data.interpolate(temp, SignalKind::Universal);
                    //IC(counter, temp, out[counter]);
                    counter++;
                    temp = temp + temp1;
                }
                else{
                    break;
                }
            }
        }
    }

    
}
