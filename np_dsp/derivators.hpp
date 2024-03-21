#pragma once

#include "icecream.hpp"


#include <cstddef>
#include <utility>

#include <npdsp_concepts.hpp>

namespace NP_DSP::ONE_D::DERIVATORS {
     enum class FinniteDifferenceType { Forward, Central, Backward };

    
    template<FinniteDifferenceType different_t>
    struct FinniteDifference {
        constexpr static bool is_derivator = true;
        constexpr static FinniteDifferenceType different_type = different_t;


        using AdditionalDataType = GENERAL::Nil;

        template<Signal DataType, Signal DerivativeType>
        static void compute(const DataType& data, DerivativeType& out, auto * nil) {
            using T = typename DerivativeType::SampleType;
            //additional_data can be nullptr
            //init first and last values

            out[0] = data[1] - data[0];
            out[data.size() - 1] = data[data.size() - 1] - data[data.size() - 2];
            if constexpr (different_type == FinniteDifferenceType::Backward) {
                for (auto i = 1; i < data.size() - 1; i++) {
                    out[i] = static_cast<T>(data[i] - data[i - 1]);
                }
            } else if constexpr (different_type == FinniteDifferenceType::Central) {
                for (auto i = 1; i < data.size() - 1; i++) {
                    out[i] = (data[i + 1] - data[i - 1]) / static_cast<T>(2.0);
                }
            } else if constexpr (different_type == FinniteDifferenceType::Forward) {
                for (auto i = 1; i < data.size() - 1; i++) {
                    out[i] = static_cast<T>(data[i + 1] - data[i]);
                }
            } else {
                std::unreachable();
            }
        }

        template<Signal DataType, Signal DerivativeType>
        static void compute(const DataType& data, DerivativeType& out, std::nullptr_t nil) {
            using T = typename DerivativeType::SampleType;
            //additional_data can be nullptr
            //init first and last values
            auto size_data = data.size();
            auto out_size = out.size();

            out[0] = data[1] - data[0];
            out[data.size() - 1] = data[data.size() - 1] - data[data.size() - 2];
            if constexpr (different_type == FinniteDifferenceType::Backward) {
                for (auto i = 1; i < data.size() - 1; i++) {
                    out[i] = static_cast<T>(data[i] - data[i - 1]);
                }
            } else if constexpr (different_type == FinniteDifferenceType::Central) {
                for (auto i = 1; i < data.size() - 1; i++) {
                    out[i] = (data[i + 1] - data[i - 1]) / static_cast<T>(2.0);
                }
            } else if constexpr (different_type == FinniteDifferenceType::Forward) {
                for (auto i = 1; i < data.size() - 1; i++) {
                    out[i] = static_cast<T>(data[i + 1] - data[i]);
                }
            } else {
                std::unreachable();
            }
        }
    };
}
