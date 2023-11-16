module;

import <cstddef>;

import npdsp_concepts;

export module derivators;


namespace NP_DSP{
    namespace ONE_D{
        namespace DERIVATORS{
            export
            template <Signal DataT, Signal DerivativeT>
            struct FinniteDifference{
                constexpr static bool is_derivator = true;
                using DataType = DataT;
                using DerivativeType = DerivativeT;

                void compute(DataType data, DerivativeType & out)
                {
                    //todo
                }
            };
        }
    }
}
