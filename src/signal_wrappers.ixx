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
            constexpr static size_t dims_count = 1;

            std::vector<T> & vec;
            inline T& getRefByIdx(size_t idx){
                return vec[idx];
            }

            inline T getByIdx(size_t idx){
                return vec[idx];
            }

            size_t getDimSize(size_t dim_number){

                return vec.size();
            }
        };

        static_assert(is_signal<SimpleVecWrapper<int>>);
    }

}