module;

#include <icecream.hpp>

export module filters;

import npdsp_concepts;
import signals;
import <utility>;
import integrators;
import <string>;

namespace NP_DSP{
    namespace ONE_D{
        namespace FILTERS{
            export enum class InstFreqKind{Average, Double}; 
            export enum class FilteringType{DerivativeBased, ValueBased};

            export
            template<Signal DataT, Signal OutT, Signal InstFreqT, FilteringType filtering_type_k,
                    Integrator IntegratorT, InstFreqKind inst_freq_k>
            struct NonOptPeriodBasedFilter{
                using DataType = DataT;
                using OutType = OutT;
                using InstFreqType = InstFreqT;
                using AdditionalDataType = GENERAL::Nil;

                InstFreqT * inst_freq;
                IntegratorT integrator;
                constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
                constexpr static bool is_filter = true;

                NonOptPeriodBasedFilter(InstFreqT & inst_freq, IntegratorT integrator){
                    this->inst_freq = &inst_freq;
                    this->integrator = integrator;
                }

                void compute(const DataType & data, OutType & out, GENERAL::Nil & additional_data){
                    if constexpr (filtering_type_k == FilteringType::DerivativeBased){
                        if constexpr (inst_freq_kind == InstFreqKind::Average){
                            auto size_expr = [&](){
                                return data.size();
                            };

                            GENERAL::Nil nil;
                            auto val_expression = [&](DataType::IdxType idx){
                                return (data.interpolate(idx + 0.5/inst_freq->interpolate(idx)) - data.interpolate(idx - 0.5/inst_freq->interpolate(idx))) *
                                inst_freq->interpolate(idx);
                            };

                            ExpressionWrapper<typename DataType::SampleType, typename DataType::IdxType, 
                                decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                    expr_wrapper (val_expression, size_expr);

                            GenericSignal<decltype(expr_wrapper), false> expr_signal(expr_wrapper);
                            
                            NP_DSP::ONE_D::INTEGRATORS::Riman<decltype(expr_signal), OutType, 
                                NP_DSP::ONE_D::INTEGRATORS::PolygonType::ByPoint> integrator_new;
                            
                            //todo fix its big crucher
                            integrator_new.compute(expr_signal, out, nil);
                        }
                        else if constexpr (inst_freq_kind == InstFreqKind::Double){
                            /*auto val_expression = [=](DataType::IdxType idx){
                                return data.interpolate(idx + 0.5/inst_freq.interpolate(idx).forward) * inst_freq.interpolate(idx).forward
                                - data.interpolate(idx - 0.5/inst_freq.interpolate(idx).backward) * inst_freq.interpolate(idx).backward;
                            };
                            integrator.compute(ExpressionWrapper<typename DataType::SampleType,typename DataType::IdxType,
                                    decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                                       (val_expression, GENERAL::Nil{}, size_expr), out, {});*/
                            std::unreachable();
                        }
                        else{
                            std::unreachable();
                        }
                    }
                    else if constexpr (filtering_type_k == FilteringType::ValueBased){
                        if constexpr (inst_freq_kind == InstFreqKind::Average){
                            auto val_expression = [&](DataType::IdxType idx){
                                return (data.interpolate(idx + 0.25/inst_freq->interpolate(idx)) 
                                    + data.interpolate(idx - 0.25/inst_freq->interpolate(idx))) / 2;
                            };
                            for (auto i = 0; i < data.size(); i++){
                                out[i] = val_expression(i);
                            }
                        }
                        else if constexpr (inst_freq_kind == InstFreqKind::Double){
                            /*auto val_expression = [=](DataType::IdxType idx){
                                return data.interpolate(idx + 0.5/inst_freq.interpolate(idx).forward) * inst_freq.interpolate(idx).forward
                                - data.interpolate(idx - 0.5/inst_freq.interpolate(idx).backward) * inst_freq.interpolate(idx).backward;
                            };
                            integrator.compute(ExpressionWrapper<typename DataType::SampleType,typename DataType::IdxType,
                                    decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                                       (val_expression, GENERAL::Nil{}, size_expr), out, {});*/
                            std::unreachable();
                        }
                        else{
                            std::unreachable();
                        }
                    }
                    else{
                        std::unreachable();
                    }
                }
            };


            export
            template<Signal DataT, Signal OutT, Signal InstFreqT, InstFreqT inst_freq,
                    Integrator IntegratorT, IntegratorT integrator, InstFreqKind inst_freq_k>
            struct OptPeriodBasedFilter{
                using DataType = DataT;
                using OutType = OutT;
                using InstFreqType = InstFreqT;
                using AdditionalDataType = GENERAL::Nil;

                constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
                constexpr static bool is_filter = true;

                void compute(DataType data, OutType & out, GENERAL::Nil & additional_data) {

                }
            };
        }
        
    }
}