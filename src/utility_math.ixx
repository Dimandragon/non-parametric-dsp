module;

import <utility>;

export module utility_math;

namespace NP_DSP{
    namespace GENERAL{
        namespace UTILITY_MATH{
            export
            template<typename xType, typename yType>
            yType linearInterpolate(std::pair<xType, yType> point1, std::pair<xType, yType> point2, xType x_in){
                auto dx = point2.first - point1.first;
                auto dy = point2.second - point1.second;
                return point1.second + dy*(x_in-point1.first)/dx;
            }

            export
            template<typename xType, typename yType>
            xType backLinearInterpolate(std::pair<xType, yType> point1, std::pair<xType, yType> point2, yType y_in){
                auto dx = point2.first - point1.first;
                auto dy = point2.second - point1.second;
                return point1.first + dx*(y_in-point1.second)/dy;
            }
        }
    }
}