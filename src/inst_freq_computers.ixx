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
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(out.getByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::TimeAverage){
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(out.getByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(out, computer_buffer);
                        //todo
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveAverage){
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(out.getByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(out, computer_buffer);
                        //todo
                    }
                    else if constexpr (counting_kind == InstFreqDerivativeBasedKind::DeriveDouble){
                        derivator.compute(data, computer_buffer, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            computer_buffer.getRefByIdx(i) = std::atan(computer_buffer.getByIdx(i));
                        }
                        derivator.compute(computer_buffer, out, nil);
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = std::abs(out.getByIdx(i)) / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(out, computer_buffer);
                        //todo
                    }

                    if constexpr (!std::is_same_v<OptFunction, GENERAL::Nil>){
                        for (int i = 0; i < data.getSize(); i++){
                            out.getRefByIdx(i) = out.getByIdx(i) * opt_function.getByIdx(i);
                        }
                    }
                }
            };
        }
    }
}