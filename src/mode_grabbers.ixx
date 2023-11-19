module;

import npdsp_concepts;
import <cstddef>;

export module mode_grabbers;

namespace NP_DSP{
    namespace ONE_D{
        export
        template<Signal DataT, Signal ModeT>
        struct Linear{
            using DataType = DataT;
            using ModeType = ModeT;
            using AdditionalDataType = GENERAL::Nil;

            constexpr static bool is_mode_grabber = true;

            void compute(DataType & data, ModeType mode, GENERAL::Nil nil)
            {
                for (std::size_t i = 0; i < data.getSize(); i++){
                    data.getRefByIdx(i) = data.getValueByIdx(i) - mode.getValueByIdx(i);
                }
            }

        };
    }
}