module;

import <vector>;
import <tuple>;
import <cassert>;
import <cstddef>;

import npdsp_concepts;

export module signal_wrappers;


namespace NP_DSP
{
    namespace ONE_D{
        export
        template<typename T>
        struct SimpleVecWrapper{
            using DataType = T;
            using IdxType = size_t;
            constexpr static bool is_signal = true;

            std::vector<T> & vec;
            inline T& getRefByIdx(size_t idx){
                return vec[idx];
            }

            inline T getByIdx(size_t idx){
                return vec[idx];
            }

            size_t getSize(){
                return vec.size();
            }
        };

        static_assert(is_signal<SimpleVecWrapper<int>>);
    }
}