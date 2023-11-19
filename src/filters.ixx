module;


import npdsp_concepts;
import signals;
import <utility>;

export module filters;

namespace NP_DSP{
    namespace ONE_D{
        export enum class InstFreqKind{Average, Double};

        export
        template<Signal DataT, Signal OutT, Signal InstFreqT, InstFreqT inst_freq,
                Integrator IntegratorT, IntegratorT integrator, InstFreqKind inst_freq_k>
        struct NonOptPeriodBasedFilter{
            using DataType = DataT;
            using OutType = OutT;
            using InstFreqType = InstFreqT;
            using AdditionalDataType = GENERAL::Nil;

            constexpr static InstFreqKind inst_freq_kind = inst_freq_k;
            constexpr static bool is_filter = true;

            void compute(DataType data, OutType & out, GENERAL::Nil & additional_data){
                auto nil = GENERAL::Nil{};

                auto size_expr = [=](){
                    return data.getSize();
                };

                if constexpr (inst_freq_kind == InstFreqKind::Average){
                    auto val_expression = [=](DataType::IdxType idx){
                        return (data.getValueByIdx(idx + 0.5/inst_freq.getValueByIdx(idx)) - data.getValueByIdx(idx - 0.5/inst_freq.getValueByIdx(idx))) *
                        inst_freq.getValueByIdx(idx);
                    };
                    integrator.compute(ExpressionWrapper<typename DataType::SampleType,typename DataType::IdxType,
                            decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                               (val_expression, GENERAL::Nil{}, size_expr), out, nil);
                }
                else if constexpr (inst_freq_kind == InstFreqKind::Double){
                    auto val_expression = [=](DataType::IdxType idx){
                        return data.getValueByIdx(idx + 0.5/inst_freq.getValueByIdx(idx).forward) * inst_freq.getValueByIdx(idx).forward
                        - data.getValueByIdx(idx - 0.5/inst_freq.getValueByIdx(idx).backward) * inst_freq.getValueByIdx(idx).backward;
                    };
                    integrator.compute(ExpressionWrapper<typename DataType::SampleType,typename DataType::IdxType,
                            decltype(val_expression), GENERAL::Nil, decltype(size_expr), false>
                                               (val_expression, GENERAL::Nil{}, size_expr), out, nil);
                }
                else{
                    std::unreachable();
                }
            }
        };


        export
        template<typename T>
        struct OptPeriodBasedFilter{
            //todo
        };
    }
}