module;

export module derivators;

import <cstddef>;
import <utility>;

import npdsp_concepts;

namespace NP_DSP{
    namespace ONE_D{
        namespace DERIVATORS{
            export enum class FinniteDifferenceType{Forward, Central, Backward};

            export
            template <Signal DataT, Signal DerivativeT, FinniteDifferenceType different_t>
            struct FinniteDifference{
                constexpr static bool is_derivator = true;
                constexpr static FinniteDifferenceType different_type = different_t;
                using DataType = DataT;
                using DerivativeType = DerivativeT;
                using AdditionalDataType = GENERAL::Nil;

                void compute(DataType data, DerivativeType & out, GENERAL::Nil & additional_data)
                {
                    //init first and last values
                    out[0] = data[1] - data[0];
                    out[data.size()-1] = data[data.size()-1] - data[data.size()-2];
                    if constexpr (different_type == FinniteDifferenceType::Backward){
                        for (std::size_t i = 1; i < data.size() - 1; i++){
                            out[i] = static_cast<typename DerivativeType::SampleType>(data[i] - data[i-1]);
                        }
                    }
                    else if constexpr (different_type == FinniteDifferenceType::Central){
                        for (std::size_t i = 1; i < data.size() - 1; i++){
                            out[i] = (data[i+1] - data[i-1]) / static_cast<typename DerivativeType::SampleType>(2.0);
                        }
                    }
                    else if constexpr (different_type == FinniteDifferenceType::Forward){
                        for (std::size_t i = 1; i < data.size() - 1; i++){
                            out[i] = static_cast<typename DerivativeType::SampleType>(data[i+1] - data[i]);
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
