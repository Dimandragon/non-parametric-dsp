module;

import <cstddef>;
import <utility>;

import npdsp_concepts;

export module integrators;

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

                void compute(DataType data, IntegralType & out)
                {
                    typename IntegralType::DataType integral = static_cast<IntegralType>(0.0);

                    if constexpr (polygon_type == PolygonType::ByPoint) {
                        for (size_t i = 0; i < data.getSize(); i++){
                            integral += static_cast<IntegralType>(data.getByIdx(i));
                            out.getRefByIdx(i) = integral;
                        }
                    }
                    else if constexpr (polygon_type == PolygonType::ByAverage){
                        integral += static_cast<IntegralType>(data.getByIdx(1) + data.getByIdx(0))/static_cast<IntegralType>(4.0);
                        out.getRefByIdx(0) = integral;
                        for (size_t i = 1; i < data.getSize()-1; i++){
                            integral += static_cast<IntegralType>(data.getByIdx(i-1) + data.getByIdx(i)*2 + data.getByIdx(i+1))
                                    / static_cast<IntegralType>(4.0);
                            out.getRefByIdx(i) = integral;
                        }
                        integral += static_cast<IntegralType>(data.getByIdx(data.getSize()-1) + data.getByIdx(data.getSize()-2))
                                    / static_cast<IntegralType>(4.0);
                        out.getRefByIdx(data.getSize()-1) = integral;
                    }
                }
            };
        }
    }
}