#pragma once

#include "derivators.hpp"
#include <cstddef>

#include <npdsp_concepts.hpp>

namespace NP_DSP::ONE_D::INTEGRATORS {
enum class PolygonType { ByPoint, ByAverage };

// численная аппроксимация интеграла по Риману
template <PolygonType polygon_t> struct Riman {
  constexpr static bool is_integrator = true;
  constexpr static PolygonType polygon_type = polygon_t;

  using AdditionalDataType = GENERAL::Nil;

  template <typename DataType, typename IntegralType>
  void compute(const DataType &data, IntegralType &out, auto *nil) {
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
        integral += static_cast<T>(data[i - 1] + data[i] * 2 + data[i + 1]) /
                    static_cast<T>(4.0);
        out[i] = integral;
      }
      integral +=
          static_cast<T>(data[data.size() - 1] + data[data.size() - 2]) /
          static_cast<T>(4.0);
      out[data.size() - 1] = integral;
    }
  }

  template <typename DataType, typename IntegralType>
  void compute(const DataType &data, IntegralType &out, std::nullptr_t nil) {
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
        integral += static_cast<T>(data[i - 1] + data[i] * 2 + data[i + 1]) /
                    static_cast<T>(4.0);
        out[i] = integral;
      }
      integral +=
          static_cast<T>(data[data.size() - 1] + data[data.size() - 2]) /
          static_cast<T>(4.0);
      out[data.size() - 1] = integral;
    }
  }
};

// Численное интегро-диффиринцирование произвольных (в том числе и дробных)
// степеней через разложение в ряд фурье
using FTDerivativeKind = DERIVATORS::FTDerivativeKind;

template <FTDerivativeKind kind_e> using FTBased = DERIVATORS::FTBased<kind_e>;
} // namespace NP_DSP::ONE_D::INTEGRATORS
