module;

export module mode_grabbers;

import npdsp_concepts;
import <cstddef>;

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
                for (std::size_t i = 0; i < data.size(); i++){
                    data[i] = data[i] - mode[i];
                }
            }

        };
    }
}