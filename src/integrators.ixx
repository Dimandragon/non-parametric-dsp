module;

export module integrators;

import <cstddef>;
import <utility>;

import npdsp_concepts;

namespace NP_DSP::ONE_D::INTEGRATORS {
    export enum class PolygonType { ByPoint, ByAverage };

    export
    template<PolygonType polygon_t>
    struct Riman {
        constexpr static bool is_integrator = true;
        constexpr static PolygonType polygon_type = polygon_t;

        using AdditionalDataType = GENERAL::Nil;

        template<Signal DataType, Signal IntegralType>
        void compute(const DataType& data, IntegralType& out, auto * nil) {
            using T = typename IntegralType::SampleType;
            T integral = static_cast<T>(0.0);

            if constexpr (polygon_type == PolygonType::ByPoint) {
                for (size_t i = 0; i < data.size(); i++) {
                    integral += static_cast<T>(data[i]);
                    out[i] = integral;
                }
            } else if constexpr (polygon_type == PolygonType::ByAverage) {
                integral += static_cast<T>(data[1] + data[0]) / static_cast<T>(4.0);
                out[0] = integral;
                for (auto i = 1; i < data.size() - 1; i++) {
                    integral += static_cast<T>(data[i - 1] + data[i] * 2 + data[i + 1])
                            / static_cast<T>(4.0);
                    out[i] = integral;
                }
                integral += static_cast<T>(data[data.size() - 1] + data[data.size() - 2])
                        / static_cast<T>(4.0);
                out[data.size() - 1] = integral;
            }
        }

        template<Signal DataType, Signal IntegralType>
        void compute(const DataType& data, IntegralType& out, std::nullptr_t nil) {
            using T = typename IntegralType::SampleType;
            T integral = static_cast<T>(0.0);

            if constexpr (polygon_type == PolygonType::ByPoint) {
                for (size_t i = 0; i < data.size(); i++) {
                    integral += static_cast<T>(data[i]);
                    out[i] = integral;
                }
            } else if constexpr (polygon_type == PolygonType::ByAverage) {
                integral += static_cast<T>(data[1] + data[0]) / static_cast<T>(4.0);
                out[0] = integral;
                for (auto i = 1; i < data.size() - 1; i++) {
                    integral += static_cast<T>(data[i - 1] + data[i] * 2 + data[i + 1])
                            / static_cast<T>(4.0);
                    out[i] = integral;
                }
                integral += static_cast<T>(data[data.size() - 1] + data[data.size() - 2])
                        / static_cast<T>(4.0);
                out[data.size() - 1] = integral;
            }
        }
    };
}
