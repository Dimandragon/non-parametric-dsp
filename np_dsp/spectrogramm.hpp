#pragma once

#include <approximators.hpp>
#include <cmath>
#include <functional>
#include <math.h>
#include <matplot/matplot.h>
#include <tuple>
#include <utility>
#include <vector>
#include <gaussian.hpp>


template <typename ExtractorT>
void spectrogramm(ExtractorT &extractor, size_t x_size, size_t y_size) {
  std::vector<std::vector<double>> data;
  for (int i = 0; i < y_size; i++) {
    data.push_back(std::vector<double>{});
    for (int j = 0; j < x_size; j++) {
      data[i].push_back(0.0);
    }
  }
  auto freqToYIdx = [&](double freq) -> int {
    auto my_log = [&](double arg) { return log2(arg) * log2(arg); };

    int idx = y_size -
              my_log(freq * extractor.getDataSize()) /
                  my_log(0.5 * extractor.getDataSize()) * (y_size - 1) -
              1;
    if (idx < 0) {
      idx = 0;
    }
    if (idx > y_size - 1) {
      idx = y_size - 1;
    }
    return idx;
  };

  auto idxToXIdx = [&](double idx) -> int {
    return idx / extractor.getDataSize() * (x_size - 1);
  };
  auto amplToZ = [&](double ampl) -> double { return std::log10(ampl); };

  for (int i = 0; i < extractor.getModesCount(); i++) {
    std::vector<double> ampl_vec = extractor.getInstAmpl(i);
    std::vector<double> freq_vec = extractor.getInstFreq(i);
    for (int j = 0; j < extractor.getDataSize(); j++) {
      data[freqToYIdx(freq_vec[j])][idxToXIdx(j)] = amplToZ(ampl_vec[j]);
    }
  }
}

enum class AmplModifier { sqrt, log2, loge, log10, linear };

// класс спектрограммы, получаеморй из экстрактора
struct Spectrogramm {
private:
  std::vector<std::tuple<double, double, double>> points = {};
  std::vector<std::vector<double>> matrix = {};
public:

  int x_size = 500;
  int y_size = 500;

  double time_size_ratio;
  double freq_size_ratio;

  double max_time;
  double max_freq;

  AmplModifier modifier = AmplModifier::linear;

  const std::vector<std::vector<double>> & getMatrix(){
    return matrix;
  }

  void loadPoints(const std::vector<double> &time,
                  const std::vector<double> &freqency,
                  const std::vector<double> &amplitude) {
    for (int i = 0; i < time.size(); i++) {
      points.push_back(std::make_tuple(time[i], freqency[i], amplitude[i]));
    }
  }

  void setBounds(double max_time, double max_freq) {
    this->max_freq = max_freq;
    this->max_time = max_time;

    time_size_ratio = max_time / x_size;
    freq_size_ratio = max_freq / y_size;
  }

  void setAxis(double x_size, double y_size) {
    this->x_size = x_size;
    this->y_size = y_size;

    time_size_ratio = max_time / x_size;
    freq_size_ratio = max_freq / y_size;
  }

  void computeMatrix() {
    if (matrix.size() != 0){
      matrix.clear();
    }
    
    for (int i = 0; i < x_size; i++) {
      matrix.push_back({});
      for (int j = 0; j < y_size; j++) {
        matrix[i].push_back(0.0);
      }
    }
    time_size_ratio = max_time / x_size;
    freq_size_ratio = max_freq / y_size;

    for (auto const &point : points) {
      double time = std::get<0>(point);
      double freq = std::get<1>(point);
      double ampl = std::abs(std::get<2>(point));

      int time_idx = time / time_size_ratio;
      int freq_idx = freq / freq_size_ratio;

      if (time_idx >= x_size) {
        time_idx = x_size - 1;
      }
      if (freq_idx >= y_size) {
        freq_idx = y_size - 1;
      }
      if (time_idx < 0) {
        time_idx = 0;
      }
      if (freq_idx < 0) {
        freq_idx = 0;
      }

      switch (modifier) {
      case AmplModifier::linear:
        break;
      case AmplModifier::log10:
        ampl = std::log10(ampl);
        break;
      case AmplModifier::log2:
        ampl = std::log2(ampl);
        break;
      case AmplModifier::loge:
        ampl = std::log(ampl);
        break;
      case AmplModifier::sqrt:
        ampl = std::sqrt(ampl);
        break;
      }

      matrix[time_idx][freq_idx] += ampl;
    }
  }

  void computeRBFBluredMatrix(double radius = 10.0, double n_layers = 20,
                                double lambda_n_s = 0.0,
                                double search_r = 1.0) {
    std::vector<std::vector<double>> x_for_load;
    std::vector<std::vector<double>> y_for_load;

    /*for (int i = 0; i < points.size(); i++){
        x_for_load.push_back({std::get<0>(point), std::get<1>(point)});
        y_for_load.push_back({std::get<2>(point)});
    }*/

    for (int i = 0; i < matrix.size(); i++) {
      for (int j = 0; j < matrix[i].size(); j++) {
        if (matrix[i][j] != 0.0) {
          x_for_load.push_back({static_cast<double>(i), static_cast<double>(j)});
          y_for_load.push_back({matrix[i][j]});
        }
      }
    }

    NP_DSP::ONE_D::APPROX::RBFBasedWithNoTrain approximator;
    approximator.kind = NP_DSP::ONE_D::APPROX::RBFKind::Gaussian;
    approximator.linterm_kind = NP_DSP::ONE_D::APPROX::LinTermKind::None;
    approximator.r_base = radius;
    approximator.n_layers = n_layers;
    approximator.lambda_n_s = lambda_n_s;
    approximator.search_r = search_r;

    approximator.loadNDData<decltype(x_for_load)>(x_for_load, y_for_load, 2, 1,
                                                  points.size());

    for (int i = 0; i < matrix.size(); i++) {
      for (int j = 0; j < matrix[i].size(); j++) {
        std::vector<double> idx = {static_cast<double>(i), static_cast<double>(j)};
        std::vector<double> value = {0.0};
        approximator.compute<decltype(idx), decltype(value)>(idx, value);
        matrix[i][j] = value[0];
      }
    }
  }

  void loadFromExtractor(auto const &extractor) {
    points.clear();
    auto max_time = extractor.getDataSize();
    auto max_freq = 0.6 * max_time;
    setAxis(extractor.getDataSize(), extractor.getDataSize());
    setBounds(max_time, max_freq);
    for (int i = 0; i < extractor.getModesCount(); i++) {
      auto freqs = extractor.getInstFreq(i);
      auto ampls = extractor.getInstAmpl(i);
      for (int j = 0; j < extractor.getDataSize(); j++) {
        auto freq = freqs[j] * max_time;
        auto ampl = ampls[j];
        auto time = j;

        points.push_back(std::make_tuple(time, freq, ampl));
      }
    }
  }

  void computeUsingExtractor(auto &extractor, const auto &data) {
    extractor.data = data;
    extractor.load(data);
    extractor.compute();
    loadFromExtractor(extractor);
  }

  double computeL2Distance(auto &other_matrix) {
    double distance = 0.0;
    for (int i = 0; i < matrix.size(); i++) {
      double distance_temp = 0.0;
      for (int j = 0; j < matrix[i].size(); j++) {
        distance_temp += (matrix[i][j] - other_matrix[i][j]) *
                         (matrix[i][j] - other_matrix[i][j]);
      }
      distance += distance_temp;
    }
    distance = std::sqrt(distance);
    return distance;
  }

  double computeRMSDistance(auto &other_matrix) {
    double distance = 0.0;
    for (int i = 0; i < matrix.size(); i++) {
      double distance_temp = 0.0;
      for (int j = 0; j < matrix[i].size(); j++) {
        distance_temp += std::sqrt((matrix[i][j] - other_matrix[i][j]) *
                                   (matrix[i][j] - other_matrix[i][j])) /
                         matrix[i].size();
      }
      distance += distance_temp / matrix.size();
    }
    return distance;
  }

  auto generateGaussian(double radius, double pow) {
    double c = radius * 2 / 2.35482;
    double a = 1.0 / (2.506628275 * c);

    int width_i = static_cast<int>(radius * 2) + 2;
    double b = 0;

    return [=](double r) {
      return a * std::pow(std::numbers::e, -((r - b) * (r - b) / (2 * c * c)));
    };
  }

  void smoothMatrix(double smoothing_radius, double gaussian_pow) {
    auto gaussian = generateGaussian(smoothing_radius, gaussian_pow);
    // todo
    std::vector<std::vector<double>> matrix_new;

    for (int i = 0; i < matrix.size(); i++) {
      matrix_new.push_back({});
      for (int j = 0; j < matrix[i].size(); j++) {
        matrix_new[i].push_back(0.0);
      }
    }

    for (int i = 0; i < matrix.size(); i++) {
      for (int j = 0; j < matrix[i].size(); j++) {
        for (int idx = -smoothing_radius; idx < smoothing_radius; idx++) {
          for (int jdx = -smoothing_radius; jdx < smoothing_radius; jdx++) {
            if (i + idx >= 0 && i + idx < matrix.size()) {
              if (j + jdx >= 0 && j + jdx < matrix[i].size()) {
                auto r = std::sqrt(idx * idx + jdx * jdx);
                matrix_new[i + idx][j + jdx] += matrix[i][j] * gaussian(r);
              }
            }
          }
        }
      }
    }
    for(int i = 0; i < matrix.size(); i++){
      for(int j = 0; j < matrix[i].size(); j++){
        matrix[i][j] = matrix_new[i][j];
      }
    }
    //std::swap(matrix, matrix_new);
  }

  void smoothMatrixFast(double smoothing_sigma, int passes){
    std::vector<std::vector<double>> matrix_new;

    for (int i = 0; i < matrix.size(); i++) {
      matrix_new.push_back({});
      for (int j = 0; j < matrix[i].size(); j++) {
        matrix_new[i].push_back(0.0);
      }
    }

    gaussianBlur<decltype(matrix)>(matrix, matrix_new, smoothing_sigma, passes);

    for(int i = 0; i < matrix.size(); i++){
      for(int j = 0; j < matrix[i].size(); j++){
        matrix[i][j] = matrix_new[i][j];
      }
    }
  }
};