#pragma once

#include <cstdint>

#include "interpolation.h"
#include <complex>
#include <makima.hpp>
#include <npdsp_concepts.hpp>
#include <npdsp_config.hpp>
#include <numbers>
#include <optional>
#include <pchip.hpp>
#include <random>
#include <string>
#include <utility>
#include <utility_math.hpp>
#include <vector>

namespace NP_DSP::ONE_D::APPROX {
enum class FSApproxKind { Simple, Positive };

// Аппроксимация рядом фурье с возможностью оптимизации функции ошибки
template <typename LossFunc, typename StopPointFunc, FSApproxKind kind_v,
          typename BySampleLoss>
struct FourierSeriesBased {
  using Loss = LossFunc;
  using StopPoint = StopPointFunc;
  static constexpr bool is_signal_approximator = true;

  static constexpr FSApproxKind kind = kind_v;

  Loss *loss;
  BySampleLoss *bySampleLoss = nullptr;

  StopPoint *stopPont;

  size_t signal_size;
  bool is_actual = false;
  int polynoms_count;

  int polynoms_count_on_tile;
  int tile_size;

  std::vector<std::complex<double>> fourier_series;
  std::vector<std::complex<double>> approximated_data;

  double max_value = 10000000000000000000000000.;

  template <typename SignalType>
  FourierSeriesBased(Loss &lossFn, SignalType &signal_in,
                     StopPointFunc &stop_point) {
    loss = &lossFn;
    signal_size = signal_in.size();
    fourier_series = std::vector<std::complex<double>>(signal_in.size());
    approximated_data = std::vector<std::complex<double>>(signal_in.size());
    stopPont = &stop_point;
    polynoms_count = signal_in.size() / 2;
    polynoms_count_on_tile = signal_in.size() / 2;
    tile_size = signal_in.size();
  }

  template <typename SignalType> void computeFSFromData(SignalType &signal_in) {
    for (size_t i = 0; i < signal_in.size(); i++) {
      approximated_data[i] = {signal_in[i], 0.0};
    }
    is_actual = true;
    for (auto i = 0; i < approximated_data.size() / tile_size; i++) {
      auto const pad = i * tile_size;
      UTILITY_MATH::fftc2c(approximated_data, fourier_series, tile_size, pad);
    }
    auto const pad = approximated_data.size() / tile_size * tile_size;
    UTILITY_MATH::fftc2c(approximated_data, fourier_series,
                         approximated_data.size() - pad, pad);
  }

  size_t mirrorIdx(size_t const i) const {
    return i / tile_size * tile_size + i / tile_size * tile_size + tile_size -
           i;
  }

  void applyMirror(size_t const i) {
    if (i % tile_size != 0) {
      auto const mirror_idx = mirrorIdx(i);
      fourier_series[mirror_idx] = {fourier_series[i].real() / 2.0,
                                    -fourier_series[i].imag() / 2.0};
      fourier_series[i] = fourier_series[i] / 2.0;
    }
  }

  void computeData() {
    is_actual = true;
    for (auto i = 0; i < approximated_data.size() / tile_size; i++) {
      auto const pad = i * tile_size;
      UTILITY_MATH::ifftc2c(fourier_series, approximated_data, tile_size, pad);
    }
    auto const pad = approximated_data.size() / tile_size * tile_size;
    UTILITY_MATH::ifftc2c(fourier_series, approximated_data,
                          approximated_data.size() - pad, pad);
  }

  void computeTile(size_t const idx) {
    auto const i = idx / tile_size;
    auto const pad = i * tile_size;
    if (approximated_data.size() > pad + tile_size) {
      UTILITY_MATH::ifftc2c(fourier_series, approximated_data, tile_size, pad);
    } else {
      UTILITY_MATH::ifftc2c(fourier_series, approximated_data,
                            approximated_data.size() - pad, pad);
    }
  }

  void setApproxOrderRatio(double const ratio) {
    polynoms_count_on_tile =
        static_cast<int>(static_cast<double>(tile_size) * 0.5 * ratio);
  }

  template <typename IdxT> std::complex<double> computeSample(IdxT idx) {
    using SampleType = double;
    std::complex<double> accum = {0.0, 0.0};

    if (idx <= approximated_data.size() && idx >= 0) {
      auto tile_first_idx = static_cast<int>(idx) / tile_size * tile_size;
      if (tile_first_idx + tile_size < approximated_data.size()) {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(tile_size);
        for (auto i = 0; i < tile_size; i++) {
          accum += fourier_series[i + tile_first_idx] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      } else {
        SampleType w =
            std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
            static_cast<SampleType>(approximated_data.size() - tile_first_idx);
        for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
          accum += fourier_series[i + tile_first_idx] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      }
    } else if (idx > approximated_data.size()) {
      auto tile_first_idx = approximated_data.size() / tile_size * tile_size;
      SampleType w =
          std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
          static_cast<SampleType>(approximated_data.size() - tile_first_idx);
      for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
        accum += fourier_series[i + tile_first_idx] *
                 std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
      }
    } else {
      if (approximated_data.size() > tile_size) {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(tile_size);
        for (auto i = 0; i < tile_size; i++) {
          accum += fourier_series[i] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      } else {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(approximated_data.size());
        for (auto i = 0; i < approximated_data.size(); i++) {
          accum += fourier_series[i] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      }
    }

    return accum;
  }

  template <typename SampleType, typename IdxType>
  SampleType compute(IdxType idx) {
    if (std::abs(static_cast<double>(idx) - static_cast<int64_t>(idx)) == 0) {
      if (is_actual) {
        return approximated_data[idx].real();
      } else {
        computeData();
        return approximated_data[idx].real();
      }
    } else {
      return computeSample(idx).real();
    }
  }

  template <typename SampleType, typename IdxType>
  std::complex<SampleType> computeComplex(IdxType idx) {
    if (std::abs(static_cast<double>(idx) - static_cast<int64_t>(idx)) == 0) {
      if (is_actual) {
        return approximated_data[idx];
      } else {
        computeData();
        return approximated_data[idx];
      }
    } else {
      return computeSample(idx);
    }
  }

  void fineTrainIter() {
    using SampleType = double;
    computeData();
    // todo for tiles and for positive
    for (auto i = 0; i < fourier_series.size(); i++) {
      if (i % tile_size > tile_size / 2) {
        // auto idx_mirror = mirrorIdx(i);
        // fourier_series[i] = {fourier_series[idx_mirror].real() / 2.0,
        // -fourier_series[idx_mirror].imag() / 2.0}; fourier_series[idx_mirror]
        // = fourier_series[idx_mirror] / 2.0;
        continue;
      }
      if (i % tile_size >= polynoms_count_on_tile) {
        continue;
      }

      auto check_loss = [&](std::pair<SampleType, SampleType> data) {
        auto const complex_sample =
            UTILITY_MATH::convertFSampleT2C<SampleType>(data);
        fourier_series[i] = complex_sample;
        applyMirror(i);
        computeTile(i);
        SampleType loss1 = 0.;
        if (bySampleLoss) {
          size_t pad = i / tile_size * tile_size;
          for (auto idx = 0; idx < tile_size; idx++) {
            if (idx + pad >= approximated_data.size()) {
              break;
            }
            loss1 += (*bySampleLoss)(*this, idx + pad);
          }
        } else {
          loss1 = (*loss)(*this);
        }
        fourier_series[i] = {complex_sample.real(), -complex_sample.imag()};
        applyMirror(i);
        computeTile(i);
        SampleType loss2 = 0.;
        if (bySampleLoss) {
          size_t const pad = i / tile_size * tile_size;
          for (auto idx = 0; idx < tile_size; idx++) {
            if (idx + pad >= approximated_data.size()) {
              break;
            }
            loss2 += (*bySampleLoss)(*this, idx + pad);
          }
        } else {
          loss2 = (*loss)(*this);
        }
        fourier_series[i] = {-complex_sample.real(), complex_sample.imag()};
        applyMirror(i);
        computeTile(i);
        SampleType loss3 = 0.;
        if (bySampleLoss) {
          size_t pad = i / tile_size * tile_size;
          for (auto idx = 0; idx < tile_size; idx++) {
            if (idx + pad >= approximated_data.size()) {
              break;
            }
            loss3 += (*bySampleLoss)(*this, idx + pad);
          }
        } else {
          loss3 = (*loss)(*this);
        }
        fourier_series[i] = {-complex_sample.real(), -complex_sample.imag()};
        applyMirror(i);
        computeTile(i);
        SampleType loss4 = 0.;
        if (bySampleLoss) {
          size_t pad = i / tile_size * tile_size;
          for (auto idx = 0; idx < tile_size; idx++) {
            if (idx + pad >= approximated_data.size()) {
              break;
            }
            loss4 += (*bySampleLoss)(*this, idx + pad);
          }
        } else {
          loss4 = (*loss)(*this);
        }

        if (loss1 <= loss2 && loss1 <= loss3 && loss1 <= loss4) {
          fourier_series[i] = complex_sample;
          applyMirror(i);
          computeTile(i);
          return loss1;
        }
        if (loss2 <= loss1 && loss2 <= loss3 && loss2 <= loss4) {
          fourier_series[i] = {complex_sample.real(), -complex_sample.imag()};
          applyMirror(i);
          computeTile(i);
          return loss2;
        }
        if (loss3 <= loss2 && loss3 <= loss1 && loss3 <= loss4) {
          fourier_series[i] = {-complex_sample.real(), complex_sample.imag()};
          applyMirror(i);
          computeTile(i);
          return loss3;
        }

        fourier_series[i] = {-complex_sample.real(), -complex_sample.imag()};
        applyMirror(i);
        computeTile(i);
        return loss4;
      };

      auto trigonometric_sample =
          UTILITY_MATH::convertFSampleC2T(fourier_series[i]);

      auto theta_max = trigonometric_sample.second + std::numbers::pi / 10.0;
      auto theta_central = trigonometric_sample.second;
      auto theta_min = trigonometric_sample.second - std::numbers::pi / 10.0;
      if (theta_max > std::numbers::pi / 2.0) {
        theta_max = std::numbers::pi / 2.0;
      }
      if (theta_min < -std::numbers::pi / 2.0) {
        theta_min = -std::numbers::pi / 2.0;
      }

      auto ampl = trigonometric_sample.first;

      auto left_loss = check_loss({static_cast<SampleType>(1.), theta_min});
      auto central_loss =
          check_loss({static_cast<SampleType>(1.), theta_central});
      auto right_loss = check_loss({static_cast<SampleType>(1.), theta_max});

      auto period_opt_iter = 0;

      while (!(*stopPont)(std::abs(right_loss - central_loss) +
                              std::abs(left_loss - central_loss),
                          *this)) {
        auto left_diff = left_loss - central_loss;
        auto right_diff = right_loss - central_loss;
        if (left_diff > 0 && right_diff <= 0) {
          left_loss = central_loss;
          theta_min = theta_central;
          theta_central = (theta_min + theta_max) / 2;
          trigonometric_sample.second = theta_central;
          central_loss = check_loss(trigonometric_sample);
        } else if (left_diff <= 0 && right_diff > 0) {
          right_loss = central_loss;
          theta_max = theta_central;
          theta_central = (theta_min + theta_max) / 2;
          trigonometric_sample.second = theta_central;
          central_loss = check_loss(trigonometric_sample);
        } else if (left_diff > right_diff) {
          // left branch is higher
          auto theta_left_avg = (theta_min + theta_central) / 2;
          trigonometric_sample.second = theta_left_avg;
          auto new_left_loss = check_loss(trigonometric_sample);
          if (new_left_loss < left_loss && new_left_loss > central_loss) {
            left_loss = new_left_loss;
            theta_min = theta_left_avg;
            theta_central = (theta_min + theta_max) / 2;
            trigonometric_sample.second = theta_central;
            central_loss = check_loss(trigonometric_sample);
          } else {
            right_loss = central_loss;
            theta_max = theta_central;
            theta_central = theta_left_avg;
            central_loss = new_left_loss;
          }
        } else if (left_diff <= right_diff) {
          // right branch is higher
          auto theta_right_avg = (theta_max + theta_central) / 2;
          trigonometric_sample.second = theta_right_avg;
          auto new_right_loss = check_loss(trigonometric_sample);
          if (new_right_loss < right_loss && new_right_loss > central_loss) {
            right_loss = new_right_loss;
            theta_max = theta_right_avg;
            theta_central = (theta_max + theta_min) / 2;
            trigonometric_sample.second = theta_central;
            central_loss = check_loss(trigonometric_sample);
          } else {
            left_loss = central_loss;
            theta_min = theta_central;
            central_loss = new_right_loss;
            theta_central = theta_right_avg;
          }
        }
        period_opt_iter++;
      }
      trigonometric_sample.second = theta_central;

      SampleType ampl_twenty_procent = ampl * 0.2;
      SampleType ampl_left = ampl;
      SampleType ampl_right = ampl;
      SampleType ampl_central = ampl;
      if (ampl_twenty_procent < 2.) {
        ampl_left -= 2.;
        ampl_right += 2.;
      } else {
        ampl_left -= ampl_twenty_procent;
        ampl_right += ampl_twenty_procent;
      }
      if (ampl_left < 0.) {
        ampl_left = 0.;
      }
      trigonometric_sample.first = ampl_left;
      auto loss_left = check_loss(trigonometric_sample);
      trigonometric_sample.first = ampl_right;
      auto loss_right = check_loss(trigonometric_sample);
      trigonometric_sample.first = ampl_central;
      auto loss_central = check_loss(trigonometric_sample);
      auto ampl_opt_iter = 0;
      auto errors_counter = 0;
      while (!(*stopPont)(std::abs(loss_right - loss_central) +
                              std::abs(loss_left - loss_central),
                          *this)) {

        auto left_diff = loss_left - loss_central;
        auto right_diff = loss_right - loss_central;
        if (left_diff > 0 && right_diff <= 0) {
          loss_left = loss_central;
          ampl_left = ampl_central;
          ampl_central = (ampl_left + ampl_right) / 2;
          trigonometric_sample.first = ampl_central;
          loss_central = check_loss(trigonometric_sample);
        } else if (left_diff <= 0 && right_diff > 0) {
          loss_right = loss_central;
          ampl_right = ampl_central;
          ampl_central = (ampl_left + ampl_right) / 2;
          trigonometric_sample.first = ampl_central;
          loss_central = check_loss(trigonometric_sample);
        } else if (left_diff > right_diff) {
          // left branch is higher
          auto ampl_left_avg = (ampl_left + ampl_central) / 2;
          trigonometric_sample.first = ampl_left_avg;
          auto new_loss_left = check_loss(trigonometric_sample);
          if (new_loss_left < loss_left && new_loss_left > central_loss) {
            loss_left = new_loss_left;
            ampl_left = ampl_left_avg;
            ampl_central = (ampl_left + ampl_right) / 2;
            trigonometric_sample.first = ampl_central;
            loss_central = check_loss(trigonometric_sample);
          } else {
            ampl_right = (ampl_central + ampl_right) / 2;
            trigonometric_sample.first = ampl_right;
            loss_right = check_loss(trigonometric_sample);
            ampl_central = (ampl_right + ampl_left) / 2;
            trigonometric_sample.first = ampl_central;
            loss_central = check_loss(trigonometric_sample);
          }
        } else if (left_diff <= right_diff) {
          auto ampl_right_avg = (ampl_right + ampl_central) / 2;
          trigonometric_sample.first = ampl_right_avg;
          auto new_loss_right = check_loss(trigonometric_sample);
          if (new_loss_right < loss_right && new_loss_right > loss_central) {
            loss_right = new_loss_right;
            ampl_right = ampl_right_avg;
            ampl_central = (ampl_left + ampl_right) / 2;
            trigonometric_sample.first = ampl_central;
            loss_central = check_loss(trigonometric_sample);
          } else {
            ampl_left = (ampl_central + ampl_left) / 2;
            trigonometric_sample.first = ampl_left;
            loss_left = check_loss(trigonometric_sample);
            ampl_central = (ampl_right + ampl_left) / 2;
            trigonometric_sample.first = ampl_central;
            loss_central = check_loss(trigonometric_sample);
          }
        }
        if (left_diff == (loss_left - loss_central)) {
          if (right_diff == (loss_right - loss_central)) {

            errors_counter++;

            if (errors_counter > 10) {
              if constexpr (CONFIG::debug) {
                for (;;) {
                }
              }
            }
          }
        }
        ampl_opt_iter++;
        if (ampl_central - ampl_left < 0.000001 ||
            ampl_right - ampl_central < 0.000001) {
          break;
        }
      }
      if (loss_central <= loss_right && loss_central <= loss_left) {
        trigonometric_sample.first = ampl_central;
      } else if (loss_right < loss_left) {
        trigonometric_sample.first = ampl_right;
      } else {
        trigonometric_sample.first = ampl_left;
      }

      check_loss(trigonometric_sample);
    }
  }

  void train() {
    using SampleType = double;
    computeData();
    for (auto i = 0; i < fourier_series.size(); i++) {
      if (i % tile_size > tile_size / 2) {
        continue;
      }

      if (i % tile_size >= polynoms_count_on_tile) {
        continue;
      }
      std::pair<SampleType, SampleType> trigonometric_sample;
      trigonometric_sample.first = 1.;

      auto check_loss = [&](std::pair<SampleType, SampleType> data) {
        auto const complex_sample =
            UTILITY_MATH::convertFSampleT2C<SampleType>(data);
        fourier_series[i] = complex_sample;
        computeTile(i);
        SampleType loss1 = 0.;
        if (bySampleLoss) {
          const size_t pad = i / tile_size * tile_size;
          for (auto idx = 0; idx < tile_size; idx++) {
            if (idx + pad >= approximated_data.size()) {
              continue;
            }
            loss1 += (*bySampleLoss)(*this, idx + pad);
          }
        } else {
          loss1 = (*loss)(*this);
        }
        fourier_series[i] = {complex_sample.real(), -complex_sample.imag()};
        computeTile(i);
        SampleType loss2 = 0.;
        if (bySampleLoss) {
          const size_t pad = i / tile_size * tile_size;
          for (auto idx = 0; idx < tile_size; idx++) {
            if (idx + pad >= approximated_data.size()) {
              continue;
            }
            loss2 += (*bySampleLoss)(*this, idx + pad);
          }
        } else {
          loss2 = (*loss)(*this);
        }
        fourier_series[i] = {-complex_sample.real(), complex_sample.imag()};
        computeTile(i);
        SampleType loss3 = 0.;
        if (bySampleLoss) {
          const size_t pad = i / tile_size * tile_size;
          for (auto idx = 0; idx < tile_size; idx++) {
            if (idx + pad >= approximated_data.size()) {
              continue;
            }
            loss3 += (*bySampleLoss)(*this, idx + pad);
          }
        } else {
          loss3 = (*loss)(*this);
        }
        fourier_series[i] = {-complex_sample.real(), -complex_sample.imag()};
        computeTile(i);
        SampleType loss4 = 0.;
        if (bySampleLoss) {
          const size_t pad = i / tile_size * tile_size;
          for (auto idx = 0; idx < tile_size; idx++) {
            if (idx + pad >= approximated_data.size()) {
              continue;
            }
            loss4 += (*bySampleLoss)(*this, idx + pad);
          }
        } else {
          loss4 = (*loss)(*this);
        }

        if (loss1 <= loss2 && loss1 <= loss3 && loss1 <= loss4) {
          fourier_series[i] = complex_sample;
          computeTile(i);
          return loss1;
        }
        if (loss2 <= loss1 && loss2 <= loss3 && loss2 <= loss4) {
          fourier_series[i] = {complex_sample.real(), -complex_sample.imag()};
          computeTile(i);
          return loss2;
        }
        if (loss3 <= loss2 && loss3 <= loss1 && loss3 <= loss4) {
          fourier_series[i] = {-complex_sample.real(), complex_sample.imag()};
          computeTile(i);
          return loss3;
        }

        fourier_series[i] = {-complex_sample.real(), -complex_sample.imag()};
        computeTile(i);
        return loss4;
      };

      auto theta_min = static_cast<SampleType>(-std::numbers::pi / 2.0 + 0.01);
      auto theta_max = static_cast<SampleType>(std::numbers::pi / 2.0 + 0.01);
      auto theta_central = static_cast<SampleType>(0.0);
      auto max_ampl = 1.;
      if constexpr (kind == FSApproxKind::Positive) {
        if (i % tile_size != 0) {
          auto thr_sample = UTILITY_MATH::convertFSampleC2T<SampleType>(
              fourier_series[i / tile_size * tile_size]);

          max_ampl = thr_sample.first * std::cos(thr_sample.second);
        }
      }
      auto left_loss =
          check_loss({static_cast<SampleType>(max_ampl), theta_min});
      auto right_loss =
          check_loss({static_cast<SampleType>(max_ampl), theta_max});
      auto central_loss =
          check_loss({static_cast<SampleType>(max_ampl), theta_central});

      auto period_opt_iter = 0;


      while (!(*stopPont)(std::abs(right_loss - central_loss) +
                              std::abs(left_loss - central_loss),
                          *this)) {

        auto left_diff = left_loss - central_loss;
        auto right_diff = right_loss - central_loss;
        if (left_diff > 0 && right_diff <= 0) {
          left_loss = central_loss;
          theta_min = theta_central;
          theta_central = (theta_min + theta_max) / 2;
          trigonometric_sample.second = theta_central;
          central_loss = check_loss(trigonometric_sample);
        } else if (left_diff <= 0 && right_diff > 0) {
          right_loss = central_loss;
          theta_max = theta_central;
          theta_central = (theta_min + theta_max) / 2;
          trigonometric_sample.second = theta_central;
          central_loss = check_loss(trigonometric_sample);
        } else if (left_diff > right_diff) {
          // left branch is higher
          auto theta_left_avg = (theta_min + theta_central) / 2;
          trigonometric_sample.second = theta_left_avg;
          auto new_left_loss = check_loss(trigonometric_sample);
          if (new_left_loss < left_loss && new_left_loss > central_loss) {
            left_loss = new_left_loss;
            theta_min = theta_left_avg;
            theta_central = (theta_min + theta_max) / 2;
            trigonometric_sample.second = theta_central;
            central_loss = check_loss(trigonometric_sample);
          } else {
            right_loss = central_loss;
            theta_max = theta_central;
            theta_central = theta_left_avg;
            central_loss = new_left_loss;
          }
        } else if (left_diff <= right_diff) {
          // right branch is higher
          auto theta_right_avg = (theta_max + theta_central) / 2;
          trigonometric_sample.second = theta_right_avg;
          auto new_right_loss = check_loss(trigonometric_sample);
          if (new_right_loss < right_loss && new_right_loss > central_loss) {
            right_loss = new_right_loss;
            theta_max = theta_right_avg;
            theta_central = (theta_max + theta_min) / 2;
            trigonometric_sample.second = theta_central;
            central_loss = check_loss(trigonometric_sample);
          } else {
            left_loss = central_loss;
            theta_min = theta_central;
            central_loss = new_right_loss;
            theta_central = theta_right_avg;
          }
        }
        period_opt_iter++;
      }
      trigonometric_sample.second = theta_central;

      SampleType ampl_left = 0.0;
      SampleType ampl_right = max_value;
      if constexpr (kind == FSApproxKind::Positive) {
        if (i % tile_size != 0) {
          ampl_right = max_ampl;
        }
      }
      SampleType ampl_central = (ampl_left + ampl_right) / 2;
      trigonometric_sample.first = ampl_left;
      auto loss_left = check_loss(trigonometric_sample);
      trigonometric_sample.first = ampl_right;
      auto loss_right = check_loss(trigonometric_sample);
      trigonometric_sample.first = ampl_central;
      auto loss_central = check_loss(trigonometric_sample);
      auto ampl_opt_iter = 0;
      auto errors_counter = 0;
      while (!((*stopPont)(std::abs(loss_right - loss_central) +
                               std::abs(loss_left - loss_central),
                           *this))) {

        auto left_diff = loss_left - loss_central;
        auto right_diff = loss_right - loss_central;
        if (left_diff > 0 && right_diff <= 0) {
          loss_left = loss_central;
          ampl_left = ampl_central;
          ampl_central = (ampl_left + ampl_right) / 2;
          trigonometric_sample.first = ampl_central;
          loss_central = check_loss(trigonometric_sample);
        } else if (left_diff <= 0 && right_diff > 0) {
          loss_right = loss_central;
          ampl_right = ampl_central;
          ampl_central = (ampl_left + ampl_right) / 2;
          trigonometric_sample.first = ampl_central;
          loss_central = check_loss(trigonometric_sample);
        } else if (left_diff > right_diff) {
          // left branch is higher
          auto ampl_left_avg = (ampl_left + ampl_central) / 2;
          trigonometric_sample.first = ampl_left_avg;
          auto new_loss_left = check_loss(trigonometric_sample);
          if (new_loss_left < loss_left && new_loss_left > central_loss) {
            loss_left = new_loss_left;
            ampl_left = ampl_left_avg;
            ampl_central = (ampl_left + ampl_right) / 2;
            trigonometric_sample.first = ampl_central;
            loss_central = check_loss(trigonometric_sample);
          } else {
            ampl_right = (ampl_central + ampl_right) / 2;
            trigonometric_sample.first = ampl_right;
            loss_right = check_loss(trigonometric_sample);
            ampl_central = (ampl_right + ampl_left) / 2;
            trigonometric_sample.first = ampl_central;
            loss_central = check_loss(trigonometric_sample);
          }
        } else if (left_diff <= right_diff) {
          auto ampl_right_avg = (ampl_right + ampl_central) / 2;
          trigonometric_sample.first = ampl_right_avg;
          auto new_loss_right = check_loss(trigonometric_sample);
          if (new_loss_right < loss_right && new_loss_right > loss_central) {
            loss_right = new_loss_right;
            ampl_right = ampl_right_avg;
            ampl_central = (ampl_left + ampl_right) / 2;
            trigonometric_sample.first = ampl_central;
            loss_central = check_loss(trigonometric_sample);
          } else {
            ampl_left = (ampl_central + ampl_left) / 2;
            trigonometric_sample.first = ampl_left;
            loss_left = check_loss(trigonometric_sample);
            ampl_central = (ampl_right + ampl_left) / 2;
            trigonometric_sample.first = ampl_central;
            loss_central = check_loss(trigonometric_sample);
          }
        }
        if (left_diff == loss_left - loss_central) {
          if (right_diff == (loss_right - loss_central)) {

            errors_counter++;

            if (errors_counter > 10) {
              if constexpr (CONFIG::debug) {
                for (;;) {
                }
              }
            }
          }
        }
        ampl_opt_iter++;
        if (ampl_central - ampl_left < 0.000001 ||
            ampl_right - ampl_central < 0.000001) {
          break;
        }
      }
      if (loss_central <= loss_right && loss_central <= loss_left) {
        trigonometric_sample.first = ampl_central;
      } else if (loss_right < loss_left) {
        trigonometric_sample.first = ampl_right;
      } else {
        trigonometric_sample.first = ampl_left;
      }

      check_loss(trigonometric_sample);
    }
  }
};

// Аппроксимация рядом фурье
struct FourierSeriesBasedWithNoTrain {
  // static constexpr bool is_signal_approximator = true;
  size_t signal_size;
  bool is_actual = false;
  int polynoms_count;

  int polynoms_count_on_tile;
  int tile_size;

  std::vector<std::complex<double>> fourier_series;
  std::vector<std::complex<double>> approximated_data;

  template <typename SignalType>
  FourierSeriesBasedWithNoTrain(const SignalType &signal_in) {
    signal_size = signal_in.size();
    fourier_series = std::vector<std::complex<double>>(signal_in.size());
    approximated_data = std::vector<std::complex<double>>(signal_in.size());
    polynoms_count = signal_in.size();
    polynoms_count_on_tile = signal_in.size();
    tile_size = signal_in.size();
  }

  template <typename SignalType> void computeFSFromData(SignalType &signal_in) {
    for (size_t i = 0; i < signal_in.size(); i++) {
      approximated_data[i] = {signal_in[i], 0.0};
    }
    is_actual = true;
    for (auto i = 0; i < approximated_data.size() / tile_size; i++) {
      auto const pad = i * tile_size;
      UTILITY_MATH::fftc2c(approximated_data, fourier_series, tile_size, pad);
    }
    auto const pad = approximated_data.size() / tile_size * tile_size;
    UTILITY_MATH::fftc2c(approximated_data, fourier_series,
                         approximated_data.size() - pad, pad);
  }

  size_t mirrorIdx(size_t const i) const {
    return i / tile_size * tile_size + i / tile_size * tile_size + tile_size -
           i;
  }

  void applyMirror(size_t const i) {
    if (i % tile_size != 0) {
      auto const mirror_idx = mirrorIdx(i);
      fourier_series[mirror_idx] = {fourier_series[i].real() / 2.0,
                                    -fourier_series[i].imag() / 2.0};
      fourier_series[i] = fourier_series[i] / 2.0;
    }
  }

  void computeData() {
    is_actual = true;
    for (auto i = 0; i < approximated_data.size() / tile_size; i++) {
      auto const pad = i * tile_size;
      UTILITY_MATH::ifftc2c(fourier_series, approximated_data, tile_size, pad);
    }
    auto const pad = approximated_data.size() / tile_size * tile_size;
    UTILITY_MATH::ifftc2c(fourier_series, approximated_data,
                          approximated_data.size() - pad, pad);
  }

  void computeTile(size_t const idx) {
    auto const i = idx / tile_size;
    auto const pad = i * tile_size;
    if (approximated_data.size() > pad + tile_size) {
      UTILITY_MATH::ifftc2c(fourier_series, approximated_data, tile_size, pad);
    } else {
      UTILITY_MATH::ifftc2c(fourier_series, approximated_data,
                            approximated_data.size() - pad, pad);
    }
  }

  void setApproxOrderRatio(double const ratio) {
    polynoms_count_on_tile =
        static_cast<int>(static_cast<double>(tile_size) * 0.5 * ratio);
  }

  template <typename IdxT> std::complex<double> computeSample(IdxT idx) {
    using SampleType = double;
    std::complex<double> accum = {0.0, 0.0};

    if (idx <= approximated_data.size() && idx >= 0) {
      auto tile_first_idx = idx / tile_size * tile_size;
      if (tile_first_idx + tile_size < approximated_data.size()) {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(tile_size);
        for (auto i = 0; i < tile_size; i++) {
          accum += fourier_series[i + tile_first_idx] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      } else {
        SampleType w =
            std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
            static_cast<SampleType>(approximated_data.size() - tile_first_idx);
        for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
          accum += fourier_series[i + tile_first_idx] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      }
    } else if (idx > approximated_data.size()) {
      auto tile_first_idx = approximated_data.size() / tile_size * tile_size;
      SampleType w =
          std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
          static_cast<SampleType>(approximated_data.size() - tile_first_idx);
      for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
        accum += fourier_series[i + tile_first_idx] *
                 std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
      }
    } else {
      if (approximated_data.size() > tile_size) {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(tile_size);
        for (auto i = 0; i < tile_size; i++) {
          accum += fourier_series[i] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      } else {
        SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                       static_cast<SampleType>(approximated_data.size());
        for (auto i = 0; i < approximated_data.size(); i++) {
          accum += fourier_series[i] *
                   std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
        }
      }
    }

    return accum;
  }

  template <typename SampleType, typename IdxType>
  SampleType compute(IdxType idx) {
    if (std::abs(static_cast<double>(idx) - static_cast<int64_t>(idx)) == 0) {
      if (is_actual) {
        return approximated_data[idx].real();
      } else {
        computeData();
        return approximated_data[idx].real();
      }
    } else {
      return computeSample(idx).real();
    }
  }

  template <typename SampleType, typename IdxType>
  std::complex<SampleType> computeComplex(IdxType idx) {
    if (std::abs(static_cast<double>(idx) - static_cast<int64_t>(idx)) == 0) {
      if (is_actual) {
        return approximated_data[idx];
      } else {
        computeData();
        return approximated_data[idx];
      }
    } else {
      return computeSample(idx);
    }
  }
};

// Интерполяция модифицированным сплайном Акимы
// see https://www.mathworks.com/help/matlab/ref/makima.html
template <typename T> struct ModifiedAkimaBasedWithNoTrain {
  // using boost::math::interpolators::makima;
  std::optional<boost::math::interpolators::makima<std::vector<double>>>
      spline = {};
  std::optional<UTILITY_MATH::SquarePolynome> square_polynom = {};
  std::optional<UTILITY_MATH::Linear> linear = {};

  double min_bound = 0.0;
  double max_bound = 0.0;

  void loadData(const T &x, const T &y) {
    std::vector<double> x_(x.size());
    std::vector<double> y_(x.size());
    for (auto i = 0; i < x.size(); i++) {
      auto const &el = x[i];
      x_[i] = el;
    }
    for (auto i = 0; i < x.size(); i++) {
      auto const &el = y[i];
      y_[i] = el;
    }
    double min = 9999999999999999999999.0;  // todo
    double max = -9999999999999999999999.0; // todo
    for (auto const &el : x_) {
      if (el < min) {
        min = el;
      }
      if (el > max) {
        max = el;
      }
    }
    min_bound = min;
    max_bound = max;

    if (x_.size() > 3) {
      spline = boost::math::interpolators::makima<std::vector<double>>(
          std::move(x_), std::move(y_));
    } else if (x_.size() == 3) {
      square_polynom = UTILITY_MATH::SquarePolynome{};
      square_polynom->solve(x_[0], x_[1], x_[2], y_[0], y_[1], y_[2]);
    } else if (x_.size() == 2) {
      linear = UTILITY_MATH::Linear{};
      linear->solve(x_[0], x_[1], y_[0], y_[1]);
    } else {
      // todo error
    }
  }

  void loadData(const T &y) {
    std::vector<double> x_(y.size());
    std::vector<double> y_(y.size());
    for (int i = 0; i < y.size(); i++) {
      x_[i] = (double)i;
    }
    for (int i = 0; i < y.size(); i++) {
      auto const &el = y[i];
      y_[i] = el;
    }
    min_bound = 0.0;
    max_bound = y.size() - 1;

    if (x_.size() > 3) {
      spline = boost::math::interpolators::makima<std::vector<double>>(
          std::move(x_), std::move(y_));
    } else if (x_.size() == 3) {
      square_polynom = UTILITY_MATH::SquarePolynome{};
      square_polynom->solve(x_[0], x_[1], x_[2], y_[0], y_[1], y_[2]);
    } else if (x_.size() == 2) {
      linear = UTILITY_MATH::Linear{};
      linear->solve(x_[0], x_[1], y_[0], y_[1]);
    } else {
      // todo error
    }
  }

  template <typename IdxT> double compute(IdxT idx) {
    if (idx >= min_bound && idx <= max_bound) {
      if (spline) {
        return (*spline)(idx);
      } else if (square_polynom) {
        return square_polynom->compute(idx);
      } else if (linear) {
        return linear->compute(idx);
      } else {
        return 0.0;
        // todo error
      }
    } else if (idx < min_bound) {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      int64_t _size = static_cast<int64_t>(size); // todo size_
      double new_idx =
          size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
      new_idx -= idx_;
      if (spline) {
        return (*spline)(new_idx);
      } else if (square_polynom) {
        return square_polynom->compute(new_idx);
      } else if (linear) {
        return linear->compute(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    } else {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      auto _size = static_cast<int64_t>(size); // todo size_
      double new_idx = _idx % _size + idx_;
      if (spline) {
        return (*spline)(new_idx);
      } else if (square_polynom) {
        return square_polynom->compute(new_idx);
      } else if (linear) {
        return linear->compute(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    }
  }

  template <typename IdxT> double computeDerive(IdxT idx) {
    if (idx >= min_bound && idx <= max_bound) {
      if (spline) {
        return spline->prime(idx);
      } else if (square_polynom) {
        return square_polynom->derive(idx);
      } else if (linear) {
        return linear->derive(idx);
      } else {
        return 0.0;
        // todo error
      }
    } else if (idx < min_bound) {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      int64_t _size = static_cast<int64_t>(size); // todo size_
      double new_idx =
          size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
      new_idx -= idx_;
      if (spline) {
        return spline->prime(new_idx);
      } else if (square_polynom) {
        return square_polynom->derive(new_idx);
      } else if (linear) {
        return linear->derive(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    } else {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      auto _size = static_cast<int64_t>(size); // todo size_
      double new_idx = _idx % _size + idx_;
      if (spline) {
        return spline->prime(new_idx);
      } else if (square_polynom) {
        return square_polynom->derive(new_idx);
      } else if (linear) {
        return linear->derive(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    }
  }

  template <typename IdxT> double computeDerive(IdxT idx, double a) {
    if (idx >= min_bound && idx <= max_bound) {
      if (spline) {
        return spline->derive(idx, a);
      } else if (square_polynom) {
        return square_polynom->derive(idx); // todo
      } else if (linear) {
        return linear->derive(idx); // todo
      } else {
        return 0.0;
        // todo error
      }
    } else if (idx < min_bound) {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      int64_t _size = static_cast<int64_t>(size); // todo size_
      double new_idx =
          size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
      new_idx -= idx_;
      if (spline) {
        return spline->derive(idx, a);
      } else if (square_polynom) {
        return square_polynom->derive(new_idx); // todo
      } else if (linear) {
        return linear->derive(new_idx); // todo
      } else {
        return 0.0;
        // todo error
      }
    } else {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      auto _size = static_cast<int64_t>(size); // todo size_
      double new_idx = _idx % _size + idx_;
      if (spline) {
        return spline->prime(new_idx);
      } else if (square_polynom) {
        return square_polynom->derive(new_idx);
      } else if (linear) {
        return linear->derive(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    }
  }
};

// Кусочно-полиномиальная монотонная сплайн интерполяция (pchip spline)
template <typename T> struct PiecewiseCubicHermitePolynomialBasedWithNoTrain {
  // using boost::math::interpolators::makima;
  std::optional<boost::math::interpolators::pchip<std::vector<double>>> spline =
      {};
  double min_bound = 0.0;
  double max_bound = 0.0;

  std::optional<UTILITY_MATH::SquarePolynome> square_polynom = {};
  std::optional<UTILITY_MATH::Linear> linear = {};

  void loadData(const T &x, const T &y) {
    square_polynom = {};
    linear = {};
    std::vector<double> x_(x.size());
    std::vector<double> y_(x.size());
    for (auto i = 0; i < x.size(); i++) {
      auto const &el = x[i];
      x_[i] = el;
    }
    for (auto i = 0; i < x.size(); i++) {
      auto const &el = y[i];
      y_[i] = el;
    }
    double min = 9999999999999999999999.0;  // todo
    double max = -9999999999999999999999.0; // todo
    for (auto const &el : x_) {
      if (el < min) {
        min = el;
      }
      if (el > max) {
        max = el;
      }
    }
    min_bound = min;
    max_bound = max;
    if (x_.size() > 3) {
      spline = boost::math::interpolators::pchip<std::vector<double>>(
          std::move(x_), std::move(y_));
    } else if (x_.size() == 3) {
      square_polynom = UTILITY_MATH::SquarePolynome{};
      square_polynom->solve(x_[0], x_[1], x_[2], y_[0], y_[1], y_[2]);
    } else if (x_.size() == 2) {
      linear = UTILITY_MATH::Linear{};
      linear->solve(x_[0], x_[1], y_[0], y_[1]);
    } else {
      // todo error
    }
  }

  void loadData(const T &y) {
    std::vector<double> x_(y.size());
    std::vector<double> y_(y.size());
    for (int i = 0; i < y.size(); i++) {
      x_[i] = (double)i;
    }
    for (int i = 0; i < y.size(); i++) {
      auto const &el = y[i];
      y_[i] = el;
    }
    min_bound = 0.0;
    max_bound = y.size() - 1;

    if (x_.size() > 3) {
      spline = boost::math::interpolators::pchip<std::vector<double>>(
          std::move(x_), std::move(y_));
    } else if (x_.size() == 3) {
      square_polynom = UTILITY_MATH::SquarePolynome{};
      square_polynom->solve(x_[0], x_[1], x_[2], y_[0], y_[1], y_[2]);
    } else if (x_.size() == 2) {
      linear = UTILITY_MATH::Linear{};
      linear->solve(x_[0], x_[1], y_[0], y_[1]);
    } else {
      // todo error
    }
  }

  template <typename IdxT> double compute(IdxT idx) {
    if (idx >= min_bound && idx <= max_bound) {
      if (spline) {
        return (*spline)(idx);
      } else if (square_polynom) {
        return square_polynom->compute(idx);
      } else if (linear) {
        return linear->compute(idx);
      } else {
        return 0.0;
        // todo error
      }
    } else if (idx < min_bound) {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      int64_t _size = static_cast<int64_t>(size); // todo size_
      double new_idx =
          size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
      new_idx -= idx_;
      if (spline) {
        return (*spline)(new_idx);
      } else if (square_polynom) {
        return square_polynom->compute(new_idx);
      } else if (linear) {
        return linear->compute(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    } else {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      auto _size = static_cast<int64_t>(size); // todo size_
      double new_idx = _idx % _size + idx_;
      if (spline) {
        return (*spline)(new_idx);
      } else if (square_polynom) {
        return square_polynom->compute(new_idx);
      } else if (linear) {
        return linear->compute(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    }
  }

  template <typename IdxT> double computeDerive(IdxT idx) {
    if (idx >= min_bound && idx <= max_bound) {
      if (spline) {
        return spline->prime(idx);
      } else if (square_polynom) {
        return square_polynom->derive(idx);
      } else if (linear) {
        return linear->derive(idx);
      } else {
        return 0.0;
        // todo error
      }
    } else if (idx < min_bound) {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      int64_t _size = static_cast<int64_t>(size); // todo size_
      double new_idx =
          size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
      new_idx -= idx_;
      if (spline) {
        return spline->prime(new_idx);
      } else if (square_polynom) {
        return square_polynom->derive(new_idx);
      } else if (linear) {
        return linear->derive(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    } else {
      int64_t _idx = static_cast<int64_t>(idx);
      double idx_ = idx - _idx;
      auto size = max_bound - min_bound;
      auto _size = static_cast<int64_t>(size); // todo size_
      double new_idx = _idx % _size + idx_;
      if (spline) {
        return spline->prime(new_idx);
      } else if (square_polynom) {
        return square_polynom->derive(new_idx);
      } else if (linear) {
        return linear->derive(new_idx);
      } else {
        return 0.0;
        // todo error
      }
    }
  }
};

// Аппроксимация методом взвешенных обратных расстояний
// see https://www.alglib.net/inverse-distance-weighting/
struct InverseDistanceWeightingBasedWithNoTrain {
  alglib::idwbuilder builder = alglib::idwbuilder();
  alglib::idwmodel model = alglib::idwmodel();
  alglib::idwreport report = alglib::idwreport();

  int layers = 15;
  double search_radius = 100;

  int n_dims_x = 1;
  int n_dims_y = 1;

  template <typename T> void loadData(const T &x, const T &y) {
    model = alglib::idwmodel();
    report = alglib::idwreport();
    builder = alglib::idwbuilder();

    alglib::idwbuildercreate(1, 1, builder);
    alglib::idwbuildersetnlayers(builder, layers);
    alglib::idwbuildersetalgomstab(builder, search_radius);
    int N = x.size();

    alglib::real_2d_array data;
    data.setlength(N, 2);

    for (int i = 0; i < N; i++) {
      data(i, 0) = x[i];
      data(i, 1) = y[i];
    }
    try{
      alglib::idwbuildersetpoints(builder, data, N);
    }
    catch (const alglib::ap_error ap_error){
      std::cout << ap_error.msg << std::endl;
      std::unreachable();
    }
    
    alglib::idwfit(builder, model, report);

    n_dims_x = 1;
    n_dims_y = 1;
  }

  template <typename T> void loadData(const T &y) {
    model = alglib::idwmodel();
    report = alglib::idwreport();
    builder = alglib::idwbuilder();

    alglib::idwbuildercreate(1, 1, builder);
    alglib::idwbuildersetnlayers(builder, layers);
    alglib::idwbuildersetalgomstab(builder, search_radius);
    int N = y.size();

    alglib::real_2d_array data;
    data.setlength(N, 2);

    for (int i = 0; i < N; i++) {
      data(i, 0) = i;
      data(i, 1) = y[i];
    }

    alglib::idwbuildersetpoints(builder, data, N);
    alglib::idwfit(builder, model, report);

    n_dims_x = 1;
    n_dims_y = 1;
  }

  // U is 2d array
  // x and y are n_elems*(n_dims_x + n_dims_y) 2d arrays
  template <typename U>
  void loadNDData(const U &x, const U &y, int n_dims_x, int n_dims_y,
                  int n_elems) {
    model = alglib::idwmodel();
    report = alglib::idwreport();
    builder = alglib::idwbuilder();

    alglib::real_2d_array data;
    int n_dims = n_dims_x + n_dims_y;
    data.setlength(n_elems, n_dims);

    for (int i = 0; i < n_elems; i++) {
      for (int j = 0; j < n_dims_x; j++) {
        data(i, j) = x[i][j];
      }
      for (int j = 0; j < n_dims_y; j++) {
        data(i, j + n_dims_x) = y[i][j];
      }
    }

    alglib::idwbuildercreate(n_dims_x, n_dims_y, builder);
    alglib::idwbuildersetnlayers(builder, layers);
    alglib::idwbuildersetalgomstab(builder, search_radius);

    alglib::idwbuildersetpoints(builder, data, n_elems);
    alglib::idwfit(builder, model, report);

    this->n_dims_x = n_dims_x;
    this->n_dims_y = n_dims_y;
  }

  template <typename IdxT> double compute(IdxT idx) {
    alglib::real_1d_array idx_;
    idx_.setlength(1);
    idx_[0] = idx;
    alglib::real_1d_array value;
    value.setlength(1);
    alglib::idwcalc(model, idx_, value);
    return value[0];
  }

  // idx and val are arrays of n_dims_x and n_dims_y sizes
  template <typename IdxT, typename ValueT>
  void compute(const IdxT &idx, ValueT &val) {
    alglib::real_1d_array idx_;
    idx_.setlength(n_dims_x);
    alglib::real_1d_array val_;
    val_.setlength(n_dims_y);
    for (int i = 0; i < n_dims_x; i++) {
      idx_[i] = idx[i];
    }
    alglib::idwcalc(model, idx_, val_);
    for (int i = 0; i < n_dims_y; i++) {
      val[i] = val_[i];
    }
  }
};

// Аппроксимация радиальными базисными функциями
enum class RBFKind { TPS, Gaussian, Bell, Multiquadric, MultiquadricAuto };
enum class LinTermKind { None, Const, Linear };
struct RBFBasedWithNoTrain {
  alglib::rbfmodel model = alglib::rbfmodel(); 
  alglib::rbfreport report = alglib::rbfreport();

  LinTermKind linterm_kind = LinTermKind::Linear;

  int n_dims_x = 1;
  int n_dims_y = 1;

  double lambda_v = 0.0;
  /*
  Only for TPS, Multiquadric or MultiquadricAuto

  lambda_v -   smoothing parameter, LambdaV>=0, defaults to 0.0:
          * LambdaV=0 means that no smoothing is applied,  i.e.  the
            spline tries to pass through all dataset points exactly
          * LambdaV>0 means that a smoothing thin  plate  spline  is
            built, with larger LambdaV corresponding to models  with
            less nonlinearities. Smoothing spline reproduces  target
            values at nodes with small error; from the  other  side,
            it is much more stable.
            Recommended values:
            * 1.0E-6 for minimal stability improving smoothing
            * 1.0E-3 a good value to start experiments; first results
              are visible
            * 1.0 for strong smoothing
  */

  double r_base = 20;
  double n_layers = 20;
  double lambda_n_s = 0.0;
  double search_r = 0.4;

  /*************************************************************************
  This function sets support radius parameter  of  hierarchical  (version 2)
  RBF constructor.

  Hierarchical RBF model achieves great speed-up  by removing from the model
  excessive (too dense) nodes. Say, if you have RBF radius equal to 1 meter,
  and two nodes are just 1 millimeter apart, you  may  remove  one  of  them
  without reducing model quality.

  Support radius parameter is used to justify which points need removal, and = alglib::rbfreport()
  which do not. If two points are less than  SUPPORT_R*CUR_RADIUS  units  of
  distance apart, one of them is removed from the model. The larger  support
  radius  is, the faster model  construction  AND  evaluation are.  However,
  too large values result in "bumpy" models.

  search_r       -   support radius coefficient, >=0.
              Recommended values are [0.1,0.4] range, with 0.1 being
              default value.

  *************************************************************************/

  /*
  Only for Gaussian or Bell
  S       -   RBF model, initialized by rbfcreate() call
  RBase   -   RBase parameter, RBase>0
  NLayers -   NLayers parameter, NLayers>0, recommended value  to  start
              with - about 5.
  LambdaNS-   >=0, nonlinearity penalty coefficient, negative values are
              not allowed. This parameter adds controllable smoothing to
              the problem, which may reduce noise. Specification of non-
              zero lambda means that in addition to fitting error solver
              will  also  minimize   LambdaNS*|S''(x)|^2  (appropriately
              generalized to multiple dimensions.

              Specification of exactly zero value means that no  penalty
              is added  (we  do  not  even  evaluate  matrix  of  second
              derivatives which is necessary for smoothing).

              Calculation of nonlinearity penalty is costly - it results
              in  several-fold  increase  of  model  construction  time.
              Evaluation time remains the same.

              Optimal  lambda  is  problem-dependent and requires  trial
              and  error.  Good  value to  start  from  is  1e-5...1e-6,
              which corresponds to slightly noticeable smoothing  of the
              function.  Value  1e-2  usually  means  that  quite  heavy
              smoothing is applied.

  TUNING ALGORITHM

  In order to use this algorithm you have to choose three parameters:
  * initial radius RBase
  * number of layers in the model NLayers
  * penalty coefficient LambdaNS

  Initial radius is easy to choose - you can pick any number  several  times
  larger  than  the  average  distance between points. Algorithm won't break
  down if you choose radius which is too large (model construction time will
  increase, but model will be built correctly).

  Choose such number of layers that RLast=RBase/2^(NLayers-1)  (radius  used
  by  the  last  layer)  will  be  smaller than the typical distance between
  points.  In  case  model  error  is  too large, you can increase number of
  layers.  Having  more  layers  will make model construction and evaluation
  proportionally slower, but it will allow you to have model which precisely
  fits your data. From the other side, if you want to  suppress  noise,  you
  can DECREASE number of layers to make your model less flexible (or specify
  non-zero LambdaNS).

  TYPICAL ERRORS

  1. Using too small number of layers - RBF models with large radius are not
     flexible enough to reproduce small variations in the  target  function.
     You  need  many  layers  with  different radii, from large to small, in
     order to have good model.

  2. Using  initial  radius  which  is  too  small.  You will get model with
     "holes" in the areas which are too far away from interpolation centers.
     However, algorithm will work correctly (and quickly) in this case.

  */

  double alpha = 10.0;
  /*
  for Multiquadric only
  f(r)=sqrt(r^2+Alpha^2) - rbf
  */

  bool v3tol = true;
  double tol = 0.0001;
  /*
  As of ALGLIB 3.20.0, version 3 models include biharmonic RBFs, thin  plate
  splines, multiquadrics.

  Version 3 models are fit  with  specialized  domain  decomposition  method
  which splits problem into smaller  chunks.  Models  with  size  less  than
  the DDM chunk size are computed nearly exactly in one step. Larger  models
  are built with an iterative linear solver. This function controls accuracy
  of the solver.

  desired precision:
          * must be non-negative
          * should be somewhere between 0.001 and 0.000001
          * values higher than 0.001 make little sense   -  you  may
            lose a lot of precision with no performance gains.
          * values below 1E-6 usually require too much time to converge,
            so they are silenly replaced by a 1E-6 cutoff value. Thus,
            zero can be used to denote 'maximum precision'.
  */

  RBFKind kind = RBFKind::TPS;

  ~RBFBasedWithNoTrain(){}

  template <typename T> void loadData(const T &x, const T &y) {
    model = alglib::rbfmodel();
    report = alglib::rbfreport();
    alglib::rbfcreate(1, 1, model);

    if (linterm_kind == LinTermKind::Linear) {
      alglib::rbfsetlinterm(model);
    } else if (linterm_kind == LinTermKind::Const) {
      alglib::rbfsetconstterm(model);
    } else if (linterm_kind == LinTermKind::None) {
      alglib::rbfsetzeroterm(model);
    }

    if (kind == RBFKind::TPS) {
      alglib::rbfsetalgothinplatespline(model, lambda_v);
    } else if (kind == RBFKind::Gaussian) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
    } else if (kind == RBFKind::Bell) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
      alglib::rbfsetv2bf(model, 1);
    } else if (kind == RBFKind::Multiquadric) {
      alglib::rbfsetalgomultiquadricmanual(model, alpha, lambda_v);
    } else if (kind == RBFKind::MultiquadricAuto) {
      alglib::rbfsetalgomultiquadricauto(model, lambda_v);
    }
    if (v3tol) {
      alglib::rbfsetv3tol(model, tol);
    }

    int N = x.size();

    alglib::real_2d_array data;
    data.setlength(N, 2);

    for (int i = 0; i < N; i++) {
      data(i, 0) = x[i];
      data(i, 1) = y[i];
    }

    alglib::rbfsetpoints(model, data, N);
    alglib::rbfbuildmodel(model, report);

    n_dims_x = 1;
    n_dims_y = 1;
  }

  template <typename T> void loadData(const T &y) {
    model = alglib::rbfmodel();
    report = alglib::rbfreport();
    alglib::rbfcreate(1, 1, model);

    if (linterm_kind == LinTermKind::Linear) {
      alglib::rbfsetlinterm(model);
    } else if (linterm_kind == LinTermKind::Const) {
      alglib::rbfsetconstterm(model);
    } else if (linterm_kind == LinTermKind::None) {
      alglib::rbfsetzeroterm(model);
    }

    if (kind == RBFKind::TPS) {
      alglib::rbfsetalgothinplatespline(model, lambda_v);
    } else if (kind == RBFKind::Gaussian) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
    } else if (kind == RBFKind::Bell) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
      alglib::rbfsetv2bf(model, 1);
    } else if (kind == RBFKind::Multiquadric) {
      alglib::rbfsetalgomultiquadricmanual(model, alpha, lambda_v);
    } else if (kind == RBFKind::MultiquadricAuto) {
      alglib::rbfsetalgomultiquadricauto(model, lambda_v);
    }
    if (v3tol) {
      alglib::rbfsetv3tol(model, tol);
    }

    int N = y.size();

    alglib::real_2d_array data;
    data.setlength(N, 2);

    for (int i = 0; i < N; i++) {
      data(i, 0) = i;
      data(i, 1) = y[i];
    }

    alglib::rbfsetpoints(model, data, N);
    alglib::rbfbuildmodel(model, report);

    n_dims_x = 1;
    n_dims_y = 1;
  }

  // U is 2d array
  // x and y are n_elems*n_dims_x and n_elems*n_dims_y 2d arrays
  template <typename U>
  void loadNDData(const U &x, const U &y, int n_dims_x, int n_dims_y,
                  int n_elems) {
    model = alglib::rbfmodel();
    report = alglib::rbfreport();
    alglib::real_2d_array data;
    int n_dims = n_dims_x + n_dims_y;
    data.setlength(n_elems, n_dims);

    for (int i = 0; i < n_elems; i++) {
      for (int j = 0; j < n_dims_x; j++) {
        data(i, j) = x[i][j];
      }
      for (int j = 0; j < n_dims_y; j++) {
        data(i, j + n_dims_x) = y[i][j];
      }
    }

    alglib::rbfcreate(n_dims_x, n_dims_y, model);

    if (linterm_kind == LinTermKind::Linear) {
      alglib::rbfsetlinterm(model);
    } else if (linterm_kind == LinTermKind::Const) {
      alglib::rbfsetconstterm(model);
    } else if (linterm_kind == LinTermKind::None) {
      alglib::rbfsetzeroterm(model);
    }

    if (kind == RBFKind::TPS) {
      alglib::rbfsetalgothinplatespline(model, lambda_v);
    } else if (kind == RBFKind::Gaussian) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
    } else if (kind == RBFKind::Bell) {
      alglib::rbfsetalgohierarchical(model, r_base, n_layers, lambda_n_s);
      alglib::rbfsetv2supportr(model, search_r);
      alglib::rbfsetv2bf(model, 1);
    } else if (kind == RBFKind::Multiquadric) {
      alglib::rbfsetalgomultiquadricmanual(model, alpha, lambda_v);
    } else if (kind == RBFKind::MultiquadricAuto) {
      alglib::rbfsetalgomultiquadricauto(model, lambda_v);
    }
    if (v3tol) {
      alglib::rbfsetv3tol(model, tol);
    }

    alglib::rbfsetpoints(model, data, n_elems);
    alglib::rbfbuildmodel(model, report);

    this->n_dims_x = n_dims_x;
    this->n_dims_y = n_dims_y;
  }

  template <typename IdxT> double compute(IdxT idx) {
    alglib::real_1d_array idx_;
    idx_.setlength(1);
    idx_[0] = idx;
    alglib::real_1d_array value;
    value.setlength(1);
    alglib::rbfcalc(model, idx_, value);
    return value[0];
  }

  // idx and val are arrays of n_dims_x and n_dims_y sizes
  template <typename IdxT, typename ValueT>
  void compute(const IdxT &idx, ValueT &val) {
    alglib::real_1d_array idx_;
    idx_.setlength(n_dims_x);
    alglib::real_1d_array val_;
    val_.setlength(n_dims_y);
    for (int i = 0; i < n_dims_x; i++) {
      idx_[i] = idx[i];
    }
    alglib::rbfcalc(model, idx_, val_);
    for (int i = 0; i < n_dims_y; i++) {
      val[i] = val_[i];
    }
  }

  template <typename IdxT> double computeDerive(const IdxT &idx) {
    // todo
    alglib::real_1d_array idx_;
    idx_.setlength(1);
    idx_[0] = idx;
    alglib::real_1d_array value;
    value.setlength(1);
    alglib::real_1d_array derivative;
    derivative.setlength(1);
    alglib::rbfdiff(model, idx_, value, derivative);
    return derivative[0];
  }

  // idx and val are arrays of n_dims_x and n_dims_y sizes
  template <typename IdxT, typename DerivativeT>
  void computeDerive(const IdxT &idx, DerivativeT &der) {
    alglib::real_1d_array idx_;
    idx_.setlength(n_dims_x);
    alglib::real_1d_array val_, der_;
    val_.setlength(n_dims_y);
    der_.setlength(n_dims_y * n_dims_x);
    for (int i = 0; i < n_dims_x; i++) {
      idx_[i] = idx[i];
    }
    alglib::rbfdiff(model, idx_, val_, der_);
    for (int i = 0; i < n_dims_y * n_dims_x; i++) {
      der[i] = der_[i];
    }
  }
};
} // namespace NP_DSP::ONE_D::APPROX
