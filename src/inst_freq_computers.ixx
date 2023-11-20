module;

import npdsp_concepts;
import signals;
import derivators;
import integrators;
import <cmath>;
import <numbers>;

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
                    auto nil = GENERAL::Nil{};
                    if constexpr (counting_kind == InstFreqDerivativeBasedKind::Momental){
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(out.getValueByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage){
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, nil);
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
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getValueByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, nil);
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
                            out_ref_forward_getter_lamda(i) = std::abs(out.getValueByIdx(i)) / (std::numbers::pi * 2.0);
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

                    if constexpr (!std::is_same_v<OptFunction, GENERAL::Nil>){
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = out.getValueByIdx(i) * opt_function.getValueByIdx(i);
                        }
                    }
                }
            };
        }
    }
}