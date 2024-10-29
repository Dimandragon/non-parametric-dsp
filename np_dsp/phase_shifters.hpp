#include "signals.hpp"
#include <approximators.hpp>
#include <cstdio>
#include <npdsp_concepts.hpp>
#include <utility>
#include <utility_math.hpp>
#include <vector>
#include <derivators.hpp>

namespace NP_DSP::ONE_D::PHASE_SHIFTERS {
struct HTBased {
  double phase_shift = 0.5 * std::numbers::pi;

  template <typename DataT, typename OutT>
  void compute(const DataT &data, OutT &out) {
    std::vector<std::complex<double>> spectre;
    std::vector<std::complex<double>> signal;

    for (int i = 0; i < data.size(); i++) {
      signal.push_back({data[i], 0.0});
      spectre.push_back({0.0, 0.0});
    }

    UTILITY_MATH::fftc2c<double>(signal, spectre);
    for (int i = 1; i < data.size(); i++) {
      std::pair<double, double> t_sample =
          UTILITY_MATH::convertFSampleC2T(spectre[i]);
      double theta = t_sample.second;
      double ampl = t_sample.first;
      // IC(theta, phase_shift);
      theta = theta + phase_shift;
      if (theta > std::numbers::pi) {
        theta = theta - std::numbers::pi;
      }
      if (theta < -std::numbers::pi) {
        theta = theta + std::numbers::pi;
      }
      // IC(theta);
      // IC(spectre[i]);
      // IC(UTILITY_MATH::convertFSampleT2C<double>({ampl, theta}));
      spectre[i] = UTILITY_MATH::convertFSampleT2C<double>({ampl, theta});
      // IC(spectre[i]);
    }
    UTILITY_MATH::ifftc2c<double>(spectre, signal);
    for (int i = 0; i < data.size(); i++) {
      // out[i] = std::sqrt(signal[i].real()*signal[i].real() +
      // signal[i].imag()*signal[i].imag());
      out[i] = signal[i].real();
      // IC(signal[i]);
    }
  }
};

/*struct FTBased{
    double phase_shift = 0.5 * std::numbers::pi;

    template<typename DataT, typename OutT>
    void compute(const DataT & data, OutT & out){
        std::vector<double> spectre;
        std::vector<double> signal;

        for (int i = 0; i < data.size(); i++){
            signal.push_back({data[i], 0.0});
            spectre.push_back({0.0, 0.0});
        }

        UTILITY_MATH::fftc2c<double>(signal, spectre);
        for (int i = 0; i < data.size(); i++){
            double imag = std::numbers::pi * i / (double)data.size() * 2.0;
            std::complex<double> muller = (0.0, std::numbers::pi * i /
data.size()); muller = std::pow(muller, power); spectre[i] = spectre[i] *
muller;
        }
        UTILITY_MATH::ifftc2c<double>(spectre, signal);
    }
};*/

struct NaiveExtremumsPhaseShifter {
  double phase_shift = 0.5 * std::numbers::pi;

  template <typename DataT, typename OutT, typename PhaseT>
  void compute(const DataT &data, OutT &out, PhaseT &phase) {}
};

template <typename DerivatorT> struct FracDiffsBasedSimple {
  double phase_shift = 0.5 * std::numbers::pi;
  DerivatorT *derivator;

  template <typename DataT, typename OutT>
  void compute(const DataT &data, OutT &out, std::nullptr_t nil) {
    IC(phase_shift);
    derivator->power = phase_shift / std::numbers::pi * 2.0;
    derivator->compute(data, out, nullptr);
    UTILITY_MATH::normalizeSTD(data, out);
  }
};

template <typename PhaseShifterT> struct WithOversampling {
  double phase_shift = 0.5 * std::numbers::pi;
  PhaseShifterT *phase_shifter;
  double oversampling_ratio = 1.0;

  template <typename DataT, typename OutT>
  void compute(const DataT &data, OutT &out, std::nullptr_t nil) {
    phase_shifter->phase_shift = phase_shift;
    APPROX::ModifiedAkimaBasedWithNoTrain<DataT> data_approx;
    data_approx.loadData(data);
    GenericSignal<SimpleVecWrapper<double>, true> data_resampled;
    data_resampled.has_ovnership = true;
    GenericSignal<SimpleVecWrapper<double>, true> out_resampled;
    out_resampled.has_ovnership = true;

    for (int i = 0; i < data.size() * oversampling_ratio; i++) {
      data_resampled.base->vec->push_back(
          data_approx.compute(double(i) / oversampling_ratio));
      out_resampled.base->vec->push_back(0.0);
    }

    phase_shifter->phase_shift = phase_shift;
    phase_shifter->compute(data_resampled, out_resampled, nullptr);

    APPROX::ModifiedAkimaBasedWithNoTrain<decltype(out_resampled)>
        out_res_approx;
    out_res_approx.loadData(out_resampled);

    for (int i = 0; i < data.size(); i++) {
      out[i] = out_res_approx.compute(double(i) * oversampling_ratio);
    }
  }
};

enum class RotateKind
{
    Naive,
    HTBased, 
    NaiveFTFracDir,
};


struct ExtremumsRotator {
    double oversampling_ratio_for_ft_der = 1.0;
    std::vector<double> extremums;
    std::vector<double> rotated_extremums;
    double phase_shift;
    RotateKind kind_e;

    void rotateExtremums(double phase_shift) {
        APPROX::PiecewiseCubicHermitePolynomialBasedWithNoTrain<std::vector<double>> approx;
        rotated_extremums.clear();
        std::vector<double> phase_x;
        std::vector<double> phase_y;
        for (int i = 0; i < extremums.size(); i++) {
          phase_x.push_back(extremums[i]);
          phase_y.push_back(i * std::numbers::pi);
        }
        approx.loadData(phase_y, phase_x);

        rotated_extremums.push_back(extremums[0]);

        for (int i = 0; i < phase_y.size() - 1; i++) {
          double temp = approx.compute(phase_y[i] + phase_shift);
          // rotated_extremums.push_back(approx.compute(phase_y[i] + phase_shift));
          if (rotated_extremums[rotated_extremums.size() - 1] == temp) {
            continue;
          }
          rotated_extremums.push_back(temp);
        }
        if (rotated_extremums[rotated_extremums.size() - 1] !=
            extremums[extremums.size() - 1]) {
          rotated_extremums.push_back(extremums[extremums.size() - 1]);
        }
    }

    template<typename DataT>
    void rotateSignalsExtremums(const DataT & data, double phase_shift){
        if (kind_e == RotateKind::Naive){
            extremums.clear();
            rotated_extremums.clear();
            UTILITY_MATH::computeExtremums<decltype(data), double>(
                data, extremums, UTILITY_MATH::ExtremumsKind::Simple);
            IC(extremums.size());
            rotateExtremums(phase_shift);
            IC(extremums.size(), rotated_extremums.size());
        }
        else if (kind_e == RotateKind::NaiveFTFracDir){
            GenericSignal<SimpleVecWrapper<double>, true> out;
            DERIVATORS::FTBased<DERIVATORS::FTDerivativeKind::Naive> ft_based1;
            PHASE_SHIFTERS::FracDiffsBasedSimple<decltype(ft_based1)> phase_shifter1;
            phase_shifter1.derivator = &ft_based1;
            PHASE_SHIFTERS::WithOversampling<decltype(phase_shifter1)> phase_shifter2;
            phase_shifter2.phase_shifter = &phase_shifter1;
            phase_shifter2.oversampling_ratio = oversampling_ratio_for_ft_der;
            phase_shifter2.phase_shift = phase_shift;

            for (int i = 0; i < data.size(); i++){
                out.base->vec->push_back(0.0);
            }

            phase_shifter2.compute(data, out, nullptr);

            UTILITY_MATH::computeExtremums<decltype(out), double>(
                out, rotated_extremums, UTILITY_MATH::ExtremumsKind::Simple);
        }
    }
};

} // namespace NP_DSP::ONE_D::PHASE_SHIFTERS
