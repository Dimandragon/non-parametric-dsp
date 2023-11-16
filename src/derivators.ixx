module;

import <cstddef>;
import <utility>;

import npdsp_concepts;

export module derivators;


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

                void compute(DataType data, DerivativeType & out)
                {
                    //init first and last values
                    out.getRefByIdx(0) = data.getByIdx(1) - data.getByIdx(0);
                    out.getRefByIdx(data.getSize()-1) = data.getByIdx(data.getSize()-1) - data.getByIdx(data.getSize()-2);
                    if constexpr (different_type == FinniteDifferenceType::Backward){
                        for (std::size_t i = 1; i < data.getSize() - 1; i++){
                            out.getRefByIdx(i) = static_cast<typename DerivativeType::DataType>(data.getByIdx(i) - data.getByIdx(i-1));
                        }
                    }
                    else if constexpr (different_type == FinniteDifferenceType::Central){
                        for (std::size_t i = 1; i < data.getSize() - 1; i++){
                            out.getRefByIdx(i) = (data.getByIdx(i+1) - data.getByIdx(i-1)) / static_cast<typename DerivativeType::DataType>(2.0);
                        }
                    }
                    else if constexpr (different_type == FinniteDifferenceType::Forward){
                        for (std::size_t i = 1; i < data.getSize() - 1; i++){
                            out.getRefByIdx(i) = static_cast<typename DerivativeType::DataType>(data.getByIdx(i+1) - data.getByIdx(i));
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
