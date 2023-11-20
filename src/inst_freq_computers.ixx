module;

import npdsp_concepts;
import signals;
import derivators;
import integrators;
import <cmath>;
import <numbers>;
import <vector>;
import <utility>;
import utility_math;

export module inst_freq_computers;

namespace NP_DSP{
    namespace ONE_D{
        namespace InstFreqComputers{
            export
            enum class InstFreqDerivativeBasedKind {Momental, TimeAverage, DeriveAverage, DeriveDouble};

            export
            template<Signal DataT, Signal OutT, SignalWrapper OptFn, OptFn opt_fn,
                    Integrator IntegratorT, Derivator DerivatorT, InstFreqDerivativeBasedKind kind>
            struct DerivativeBased{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = OutT;

                constexpr static bool is_inst_freq_computer = true;
                constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

                using OptFunction = OptFn;
                OptFunction opt_function = opt_fn;

                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;

                static_assert(OutType::is_writable == true);

                DerivativeBased(IntegratorT integrator_o,
                                DerivatorT derivator_o)
                {
                    integrator = integrator_o;
                    derivator = derivator_o;
                }

                void compute(DataType data, OutType & out, AdditionalDataType & computer_buffer)
                {
                    if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental){
                        derivator.compute(data, computer_buffer, {});
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, {});
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(out.getValueByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage){
                        derivator.compute(data, computer_buffer, {});
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, {});
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(out.getValueByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(out, computer_buffer);
                        for (auto i = 0; i < data.getSize(); i++){
                            auto approx_answer = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer = approx_answer;
                            auto counter = 0;
                            while (approx_answer < 2.0 * std::numbers::pi){
                                counter++;

                                old_approx_answer = approx_answer;
                                approx_answer += std::abs(computer_buffer.interpolate(i+counter)-computer_buffer.interpolate(i+counter-1)) +
                                        std::abs(computer_buffer.interpolate(i-counter)-computer_buffer.interpolate(i-counter+1));
                            }
                            auto left_loss = std::numbers::pi * 2.0 - old_approx_answer;
                            auto right_loss = approx_answer - std::numbers::pi * 2.0;
                            auto sum_loss = left_loss + right_loss;
                            auto period = (static_cast<OutType::SampleType>(counter) - right_loss/sum_loss) * 2;
                            out.getRefByIdx(i) = static_cast<OutType::SampleType>(1.0/period);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage){
                        derivator.compute(data, computer_buffer, {});
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, {});
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(out.getValueByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(out, computer_buffer);
                        auto right_edge = computer_buffer.findIncr(computer_buffer.getByIdx(0) + std::numbers::pi, {0}, {});
                        auto left_edge = computer_buffer.findIncr(computer_buffer.getByIdx(0) - std::numbers::pi, {}, {0});
                        auto period = right_edge + left_edge;
                        out.getRefByIdx(0) = static_cast<OutType::SampleType>(1.0/period);
                        for (auto i = 1; i < data.getSize(); i++){
                            right_edge = computer_buffer.findIncr(computer_buffer.getByIdx(i) + std::numbers::pi, {right_edge}, {});
                            left_edge = computer_buffer.findIncr(computer_buffer.getByIdx(i) - std::numbers::pi, {left_edge}, {});
                            period = right_edge + left_edge;
                            out.getRefByIdx(i) = static_cast<OutType::SampleType>(1.0/period);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveDouble){
                        derivator.compute(data, computer_buffer, {});
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getValueByIdx(i));
                        }
                        auto out_ref_forward_getter_lamda = [&out](OutType::IdxType idx){
                            auto & ref = out.getRefByIdx(idx);
                            auto & out_ref = ref.forward;
                            return out_ref;
                        };
                        auto out_ref_backward_getter_lamda = [&out](OutType::IdxType idx){
                            auto & ref = out.getRefByIdx(idx);
                            auto & out_ref = ref.forward;
                            return out_ref;
                        };
                        auto out_val_getter_lamda = [&out](OutType::IdxType idx){
                            return out.getValueByIdx(idx).forward;
                        };
                        auto get_size_lamda = [out](){
                            return out.getSize();
                        };
                        auto signal_base = ExpressionWrapper<typename OutType::SampleType,
                            typename OutType::IdxType, decltype(out_val_getter_lamda), decltype(out_ref_forward_getter_lamda),
                            decltype(get_size_lamda), true>(out_val_getter_lamda, out_ref_forward_getter_lamda, get_size_lamda);
                        auto signal_for_compute = GenericSignal(signal_base);

                        derivator.compute(computer_buffer, signal_for_compute, {});
                        for (int i = 0; i < data.getSize(); i++){
                            out_ref_forward_getter_lamda(i) = std::abs(out_val_getter_lamda(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(signal_for_compute, computer_buffer);

                        auto right_edge = computer_buffer.findIncr(computer_buffer.getByIdx(0) + std::numbers::pi, {0}, {});
                        auto left_edge = computer_buffer.findIncr(computer_buffer.getByIdx(0) - std::numbers::pi, {}, {0});
                        out_ref_forward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/left_edge);
                        out_ref_backward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/right_edge);
                        for (auto i = 1; i < data.getSize(); i++){
                            right_edge = computer_buffer.findIncr(computer_buffer.getByIdx(i) + std::numbers::pi, {right_edge}, {});
                            left_edge = computer_buffer.findIncr(computer_buffer.getByIdx(i) - std::numbers::pi, {left_edge}, {});
                            out_ref_forward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/left_edge);
                            out_ref_backward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/right_edge);
                        }
                    }
                }
            };

            export
            template<Signal DataT, Signal OutT, SignalWrapper OptFn, OptFn opt_fn,
                    Integrator IntegratorT, Derivator DerivatorT, InstFreqDerivativeBasedKind kind>
            struct DerivativeBasedWithExternalOptParametr{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = OutT;

                constexpr static bool is_inst_freq_computer = true;
                constexpr static InstFreqDerivativeBasedKind counting_kind = kind;

                using OptFunction = OptFn;
                OptFunction opt_function = opt_fn;

                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;

                static_assert(OutType::is_writable == true);

                DerivativeBasedWithExternalOptParametr(IntegratorT integrator_o,
                                DerivatorT derivator_o)
                {
                    integrator = integrator_o;
                    derivator = derivator_o;
                }

                void compute(DataType data, OutType & out, AdditionalDataType & computer_buffer)
                {
                    if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage){
                        derivator.compute(data, out, {});
                        std::vector<typename OutType::SampleType> memory_buffer;
                        for (auto i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::atan(out.getValueByIdx(i)) * computer_buffer.getValueByIdx(i);
                            memory_buffer.push_back(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(out, computer_buffer, {});
                        for (auto i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(computer_buffer.getValueByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(out, computer_buffer);
                        for (auto i = 0; i < data.getSize(); i++){
                            auto approx_answer = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer = approx_answer;
                            auto counter = 0;
                            while (approx_answer < 2.0 * std::numbers::pi){
                                counter++;

                                old_approx_answer = approx_answer;
                                approx_answer += std::abs(computer_buffer.interpolate(i+counter)-computer_buffer.interpolate(i+counter-1)) +
                                                 std::abs(computer_buffer.interpolate(i-counter)-computer_buffer.interpolate(i-counter+1));
                            }
                            auto left_loss = std::numbers::pi * 2.0 - old_approx_answer;
                            auto right_loss = approx_answer - std::numbers::pi * 2.0;
                            auto sum_loss = left_loss + right_loss;
                            auto period = (static_cast<OutType::SampleType>(counter) - right_loss/sum_loss) * 2;
                            out.getRefByIdx(i) = static_cast<OutType::SampleType>(1.0/period);
                        }
                        for (auto i = 0; i < memory_buffer.size(); i++){
                            computer_buffer.getRefByIdx(i) = memory_buffer[i];
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage){
                        derivator.compute(data, out, {});
                        std::vector<typename OutType::SampleType> memory_buffer;
                        for (auto i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::atan(out.getValueByIdx(i)) * computer_buffer.getValueByIdx(i);
                            memory_buffer.push_back(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(out, computer_buffer, {});
                        for (auto i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(computer_buffer.getValueByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(out, computer_buffer);
                        auto right_edge = computer_buffer.findIncr(computer_buffer.getByIdx(0) + std::numbers::pi, {0}, {});
                        auto left_edge = computer_buffer.findIncr(computer_buffer.getByIdx(0) - std::numbers::pi, {}, {0});
                        auto period = right_edge + left_edge;
                        out.getRefByIdx(0) = static_cast<OutType::SampleType>(1.0/period);
                        for (auto i = 1; i < data.getSize(); i++){
                            right_edge = computer_buffer.findIncr(computer_buffer.getByIdx(i) + std::numbers::pi, {right_edge}, {});
                            left_edge = computer_buffer.findIncr(computer_buffer.getByIdx(i) - std::numbers::pi, {left_edge}, {});
                            period = right_edge + left_edge;
                            out.getRefByIdx(i) = static_cast<OutType::SampleType>(1.0/period);
                        }
                        for (auto i = 0; i < memory_buffer.size(); i++){
                            computer_buffer.getRefByIdx(i) = memory_buffer[i];
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveDouble){
                        std::vector<typename OutType::SampleType> memory_buffer;
                        auto out_ref_forward_getter_lamda = [&out](OutType::IdxType idx){
                            auto & ref = out.getRefByIdx(idx);
                            auto & out_ref = ref.forward;
                            return out_ref;
                        };
                        auto out_ref_backward_getter_lamda = [&out](OutType::IdxType idx){
                            auto & ref = out.getRefByIdx(idx);
                            auto & out_ref = ref.forward;
                            return out_ref;
                        };
                        auto out_val_getter_lamda = [&out](OutType::IdxType idx){
                            return out.getValueByIdx(idx).forward;
                        };
                        auto get_size_lamda = [out](){
                            return out.getSize();
                        };
                        auto signal_base = ExpressionWrapper<typename OutType::SampleType,
                                typename OutType::IdxType, decltype(out_val_getter_lamda), decltype(out_ref_forward_getter_lamda),
                                decltype(get_size_lamda), true>(out_val_getter_lamda, out_ref_forward_getter_lamda, get_size_lamda);
                        auto signal_for_compute = GenericSignal(signal_base);

                        derivator.compute(data, signal_for_compute, {});

                        for (int i = 0; i < data.getSize(); i++){
                            out_ref_forward_getter_lamda(i) = std::atan(out.out_ref_forward_getter_lamda(i)) * computer_buffer.getValueByIdx(i);
                            memory_buffer.push_back(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(signal_for_compute, computer_buffer, {});
                        for (int i = 0; i < data.getSize(); i++){
                            out_ref_forward_getter_lamda(i) = std::abs(computer_buffer.getValueByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(signal_for_compute, computer_buffer);

                        auto right_edge = computer_buffer.findIncr(computer_buffer.getByIdx(0) + std::numbers::pi, {0}, {});
                        auto left_edge = computer_buffer.findIncr(computer_buffer.getByIdx(0) - std::numbers::pi, {}, {0});
                        out_ref_forward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/left_edge);
                        out_ref_backward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/right_edge);
                        for (auto i = 1; i < data.getSize(); i++){
                            right_edge = computer_buffer.findIncr(computer_buffer.getByIdx(i) + std::numbers::pi, {right_edge}, {});
                            left_edge = computer_buffer.findIncr(computer_buffer.getByIdx(i) - std::numbers::pi, {left_edge}, {});
                            out_ref_forward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/left_edge);
                            out_ref_backward_getter_lamda(0) = static_cast<OutType::SampleType>(1.0/right_edge);
                        }
                        for (auto i = 0; i < memory_buffer.size(); i++){
                            computer_buffer.getRefByIdx(i) = memory_buffer[i];
                        }
                    }
                }
            };

            enum class ExtremumsBasedComputeInstFreqKind {Simple, Linear};

            template<Signal DataT, Signal OutT, ExtremumsBasedComputeInstFreqKind compute_kind>
            struct ExtremumsBased{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = GENERAL::Nil;

                constexpr static ExtremumsBasedComputeInstFreqKind kind = compute_kind;
                constexpr static bool is_inst_freq_computer = true;

                static_assert(OutT::is_writable == true);

                void compute(DataType data, OutType out, GENERAL::Nil nil){
                    std::vector<typename DataType::IdxType> extremums;
                    extremums.push_back(static_cast<typename DataType::IdxType>(0));
                    for (auto i = 1; i < data.getSize()-1; i++){
                        if ((data.getValueByIdx(i) >= data.getValueByIdx(i-1) &&
                            data.getValueByIdx(i) > data.getValueByIdx(i+1)) ||
                            (data.getValueByIdx(i) > data.getValueByIdx(i-1) &&
                            data.getValueByIdx(i) >= data.getValueByIdx(i+1)) ||
                            (data.getValueByIdx(i) <= data.getValueByIdx(i-1) &&
                             data.getValueByIdx(i) < data.getValueByIdx(i+1)) ||
                            (data.getValueByIdx(i) < data.getValueByIdx(i-1) &&
                             data.getValueByIdx(i) <= data.getValueByIdx(i+1)))
                        {
                            extremums.push_back(i);
                        }
                    }
                    extremums.push_back(static_cast<typename DataType::IdxType>(data.getSize()-1));
                    if constexpr (kind == ExtremumsBasedComputeInstFreqKind::Simple){
                        auto left = extremums[0];
                        auto right = extremums[1];
                        auto counter = 0;
                        for (auto i = 0; i < data.getSize(); i++){
                            if (i>right){
                                counter++;
                                left = extremums[counter];
                                right = extremums[counter+1];
                            }
                            out.getRefByIdx(i) = static_cast<OutType::SampleType>(0.5/(right - left));
                        }
                    }
                    else if constexpr (kind == ExtremumsBasedComputeInstFreqKind::Linear){
                        std::vector<std::pair<typename OutType::SampleType, typename OutType::SampleType>> points;
                        points.push_back({static_cast<OutType::SampleType>(0), static_cast<OutType::SampleType>(0.5/(extremums[1] - extremums[0]))});
                        for (auto i = 0; i < extremums.size()-1; i++){
                            points.push_back({static_cast<OutType::SampleType>((extremums[i+1] + extremums[i])/2.0),
                                              static_cast<OutType::SampleType>(0.5/(extremums[i+1] - extremums[i]))});
                        }
                        points.push_back({static_cast<OutType::SampleType>(data.getSize()-1), static_cast<OutType::SampleType>
                                    (0.5/(extremums[extremums.size()-1] - extremums[extremums.size()-2]))});

                        auto left = points[0];
                        auto right = points[1];
                        auto counter = 0;

                        for (auto i = 0; i < out.getSize(); i++){
                            if (i>right.first){
                                counter++;
                                left = points[counter];
                                right = points[counter+1];
                            }
                            out.getRefByIdx(i) =
                                    GENERAL::UTILITY_MATH::linearInterpolate(left, right,static_cast<OutType::SampleType>(i));
                        }
                    }
                    else{
                        std::unreachable();
                    }
                }
            };

        }
    }
}