#pragma once

#include "utility_math.hpp"

#include <cstddef>
#include <utility>
#include <vector>
#include <complex>
#include <npdsp_concepts.hpp>

namespace NP_DSP::ONE_D::DERIVATORS {
    enum class FinniteDifferenceType { Forward, Central, Backward };
    
    template<FinniteDifferenceType different_t>
    struct FinniteDifference {
        constexpr static bool is_derivator = true;
        constexpr static FinniteDifferenceType different_type = different_t;


        using AdditionalDataType = GENERAL::Nil;

        template<typename DataType, typename DerivativeType>
        static void compute(const DataType& data, DerivativeType& out, auto * nil) {
            using T = double;
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
                /*std::unreachable();*/
            }
        }

        template<typename DataType, typename DerivativeType>
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
                /*std::unreachable();*/
            }
        }
    };

    struct FTBased{
        constexpr static bool is_derivator = true;
        using AdditionalDataType = GENERAL::Nil;
        double power = 1.0;

        template<typename DataT, typename OutT>
        void compute(const DataT & data, OutT & out, std::nullptr_t nil){
            std::vector<std::complex<double>> signal;
            std::vector<std::complex<double>> spectre;
            for (int i = 0; i < data.size(); i++){
                signal.push_back({data[i], 0.0});
                spectre.push_back({0.0, 0.0});
            }
            UTILITY_MATH::fftc2c<double>(signal, spectre);
            for (int i = 1; i < data.size(); i++){
                double imag = std::numbers::pi * i / (double)data.size() * 2.0;
                std::complex<double> muller = {0.0, std::numbers::pi * i / data.size()};
                muller = std::pow(muller, power);
                spectre[i] = spectre[i] * muller;
            }
            UTILITY_MATH::ifftc2c<double>(spectre, signal);
            for (int i = 0; i < data.size(); i++){
                out[i] = signal[i].real();
                //out[i] = std::sqrt(signal[i].real()*signal[i].real() + signal[i].imag()*signal[i].imag());
            }
        }
    };

    struct RiemanLoiouville{

    };
    
    struct Caputo{

    };
    // struct 
}
