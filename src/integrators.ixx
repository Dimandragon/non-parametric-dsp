module;
#include <icecream.hpp>

export module integrators;

import <cstddef>;
import <utility>;

import npdsp_concepts;

namespace NP_DSP{
    namespace ONE_D{
        namespace INTEGRATORS{
            export enum class PolygonType {ByPoint, ByAverage};

            export
            template <Signal DataT, Signal IntegralT, PolygonType polygon_t>
            struct Riman{
                constexpr static bool is_integrator = true;
                constexpr static PolygonType polygon_type = polygon_t;

                using DataType = DataT;
                using IntegralType = IntegralT;
                using AdditionalDataType = GENERAL::Nil;

                void compute(const DataType & data, IntegralType & out, GENERAL::Nil & additional_data)
                {
                    typename IntegralType::SampleType integral = static_cast<IntegralType::SampleType>(0.0);

                    if constexpr (polygon_type == PolygonType::ByPoint) {
                        for (size_t i = 0; i < data.size(); i++){
                            integral += static_cast<IntegralType::SampleType>(data[i]);
                            out[i] = integral;
                        }
                    }
                    else if constexpr (polygon_type == PolygonType::ByAverage){
                        integral += static_cast<IntegralType::SampleType>(data[1] + data[0])/static_cast<IntegralType::SampleType>(4.0);
                        out[0] = integral;
                        for (auto i = 1; i < data.size()-1; i++){
                            integral += static_cast<IntegralType::SampleType>(data[i-1] + data[i]*2 + data[i+1])
                                    / static_cast<IntegralType::SampleType>(4.0);
                            out[i] = integral;
                        }
                        integral += static_cast<IntegralType::SampleType>(data[data.size()-1] + data[data.size()-2])
                                    / static_cast<IntegralType::SampleType>(4.0);
                        out[data.size()-1] = integral;
                    }
                }
            };
        }
    }
}