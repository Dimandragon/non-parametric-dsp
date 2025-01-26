#include "filters.hpp"
#include "matplot/freestanding/plot.h"
#include "npdsp_concepts.hpp"
#include <approximators.hpp>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <magic_enum.hpp>
#include <map>
#include <math.h>
#include <matplot/matplot.h>
#include <memory>
#include <modes_extractors.hpp>
#include <numbers>
#include <ostream>
#include <spectrogramm.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <icecream.hpp>
#include <nlohmann/json.hpp>
#include <fstream>

template <typename T>
void addRandomTrigMode(T &data, double min_ampl, double max_ampl) {
  auto size = data.size();
  //IC(size);
  double period = static_cast<double>(std::rand() % (size * 100 / 2)) / 100;
  if (period < 2.0) {
    period = 2.0;
  }
  double d_phase = 6.28 / period;

  double phase = 0.0;
  double ampl =
      static_cast<double>(std::rand() % static_cast<int>((max_ampl * 100))) /
      100;
  if (ampl < min_ampl) {
    ampl = min_ampl;
  }
  phase += static_cast<double>(std::rand()) / 100;
  for (int i = 0; i < data.size(); i++) {
    data[i] += std::sin(phase) * ampl;
    phase += d_phase;
  }
}

template <typename T>
void addRandomTrigModeWithNonStationaryAmpl(T &data, double min_ampl,
                                            double max_ampl,
                                            int ampl_points_count) {
  auto size = data.size();
  double period = static_cast<double>(std::rand() % (size * 100 / 2)) / 100;
  if (period < 2.0) {
    period = 2.0;
  }
  double d_phase = 6.28 / period;

  double phase = 0.0;

  std::vector<double> ampl_points_x = {};
  std::vector<double> ampl_points_y = {};
  ampl_points_x.push_back(0.0);
  ampl_points_y.push_back((std::rand() % int(max_ampl * 100)) / 100.0);

  for (int i = 1; i < ampl_points_count; i++) {
    ampl_points_x.push_back(ampl_points_x[i - 1] + std::rand());
    ampl_points_y.push_back((std::rand() % int(max_ampl * 100)) / 100.0);
  }

  for (int i = 0; i < ampl_points_x.size(); i++) {
    if (ampl_points_y[i] < min_ampl) {
      ampl_points_y[i] = min_ampl;
    }
    ampl_points_x[i] = ampl_points_x[i] /
                       ampl_points_x[ampl_points_x.size() - 1] * data.size();
  }
  NP_DSP::ONE_D::APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>>
      ampl_approximator;
  ampl_approximator.loadData(ampl_points_x, ampl_points_y);

  phase += static_cast<double>(std::rand()) / 100;
  for (int i = 0; i < data.size(); i++) {
    data[i] += std::sin(phase) * ampl_approximator.compute(i);
    phase += d_phase;
  }
}

template <typename T>
void addRandomTrigModeWithNonStationaryFreq(T &data, double min_ampl,
                                            double max_ampl,
                                            int freq_points_count) {
  auto size = data.size();
  double d_phase_max = 6.28 / (2);
  double d_phase_min = 6.28 / (size / 2);
  double ampl =
      static_cast<double>(std::rand() % static_cast<int>((max_ampl * 100))) /
      100;
  if (ampl < min_ampl) {
    ampl = min_ampl;
  }

  auto generateDPHASE = [&]() {
    double period = static_cast<double>(std::rand() % (size * 100 / 2)) / 100;
    if (period < 2) {
      period = 2.0;
    }
    double d_phase = 6.28 / period;
    return d_phase;
  };

  std::vector<double> d_phase_points_x = {0.0};
  std::vector<double> d_phase_points_y = {generateDPHASE()};

  for (int i = 1; i < freq_points_count; i++) {
    d_phase_points_x.push_back(d_phase_points_x[i - 1] + std::rand());
    d_phase_points_y.push_back(generateDPHASE());
  }

  for (int i = 0; i < freq_points_count; i++) {
    d_phase_points_x[i] = d_phase_points_x[i] /
                          d_phase_points_x[d_phase_points_x.size() - 1] *
                          data.size();
  }

  NP_DSP::ONE_D::APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>>
      d_phase_approx;
  d_phase_approx.loadData(d_phase_points_x, d_phase_points_y);

  double phase = 0.0;
  if (ampl < min_ampl) {
    ampl = min_ampl;
  }

  phase += static_cast<double>(std::rand()) / 100;
  for (int i = 0; i < data.size(); i++) {
    data[i] += std::sin(phase) * ampl;
    phase += d_phase_approx.compute(i);
  }
}

template <typename T>
void addRandomTrigModeWithNonStationaryAmplAndFreq(T &data, double min_ampl,
                                                   double max_ampl,
                                                   int ampl_points_count,
                                                   int freq_points_count) {
  auto size = data.size();
  double d_phase_max = 6.28 / (2);
  double d_phase_min = 6.28 / (size / 2);

  auto generateDPHASE = [&]() {
    double period = static_cast<double>(std::rand() % (size * 100 / 2)) / 100;
    if (period < 2) {
      period = 2.0;
    }
    double d_phase = 6.28 / period;
    return d_phase;
  };

  std::vector<double> d_phase_points_x = {0.0};
  std::vector<double> d_phase_points_y = {generateDPHASE()};

  for (int i = 1; i < freq_points_count; i++) {
    d_phase_points_x.push_back(d_phase_points_x[i - 1] + std::rand());
    d_phase_points_y.push_back(generateDPHASE());
  }

  for (int i = 0; i < freq_points_count; i++) {
    d_phase_points_x[i] = d_phase_points_x[i] /
                          d_phase_points_x[d_phase_points_x.size() - 1] *
                          data.size();
  }

  NP_DSP::ONE_D::APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>>
      d_phase_approx;
  d_phase_approx.loadData(d_phase_points_x, d_phase_points_y);

  double phase = 0.0;
  double ampl =
      static_cast<double>(std::rand() % static_cast<int>((max_ampl * 100))) /
      100;
  if (ampl < min_ampl) {
    ampl = min_ampl;
  }

  std::vector<double> ampl_points_x = {};
  std::vector<double> ampl_points_y = {};
  ampl_points_x.push_back(0.0);
  ampl_points_y.push_back((std::rand() % int(max_ampl * 100)) / 100.0);

  for (int i = 1; i < ampl_points_count; i++) {
    ampl_points_x.push_back(ampl_points_x[i - 1] + std::rand());
    ampl_points_y.push_back((std::rand() % int(max_ampl * 100)) / 100.0);
  }

  for (int i = 0; i < ampl_points_x.size(); i++) {
    if (ampl_points_y[i] < min_ampl) {
      ampl_points_y[i] = min_ampl;
    }
    ampl_points_x[i] = ampl_points_x[i] /
                       ampl_points_x[ampl_points_x.size() - 1] * data.size();
  }
  NP_DSP::ONE_D::APPROX::ModifiedAkimaBasedWithNoTrain<std::vector<double>>
      ampl_approximator;
  ampl_approximator.loadData(ampl_points_x, ampl_points_y);

  phase += static_cast<double>(std::rand()) / 100;
  for (int i = 0; i < data.size(); i++) {
    data[i] += std::sin(phase) * ampl_approximator.compute(i);
    phase += d_phase_approx.compute(i);
  }
}

void copyMatrix(const std::vector<std::vector<double>> & source, std::vector<std::vector<double>> & dist){
  dist.clear();
  for (int i = 0; i < source.size(); i++){
    dist.push_back({});
    for (int j = 0; j < source[i].size(); j++){
      dist[i].push_back(source[i][j]);
    }
  }
}

void copyMatrix(const std::vector<std::shared_ptr<std::vector<double>>> & source, std::vector<std::vector<double>> & dist){
  dist.clear();
  for (int i = 0; i < source.size(); i++){
    dist.push_back({});
    for (int j = 0; j < source[i]->size(); j++){
      dist[i].push_back((*source[i])[j]);
    }
  }
}


int main() {
  int SIZE = 1000;
  int smoothing_sigma = 3;

  /*NP_DSP::ONE_D::FILTERS::InterpolationKind some_enum =
      NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFMultiquadricAuto;

  std::cout << magic_enum::enum_name(some_enum) << " "
            << static_cast<int>(some_enum) << std::endl;*/
  std::map<std::string,
           std::shared_ptr<NP_DSP::ONE_D::MODES_EXTRACTORS::SOTAEMD>>
      extractors;

  for (int phase_shiting_kind = 0; phase_shiting_kind < 2;
       phase_shiting_kind++) {
    for (int interpolation_kind = 0; interpolation_kind < 7;
         interpolation_kind++) {
      for (int full_pi_flag = 1; full_pi_flag <= 1; full_pi_flag++) {
        for (int filtering_iter_numbers = 0; filtering_iter_numbers < 10;
             filtering_iter_numbers++) {
          for (int phase_shifts_count = 1; phase_shifts_count < 10;
               phase_shifts_count++) {
            if (interpolation_kind == (int)NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFBell){
              continue;
            }
            if (interpolation_kind == (int)NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFGaussian){
              continue;
            }
            if (interpolation_kind == (int)NP_DSP::ONE_D::FILTERS::InterpolationKind::IDW){
              continue;
            }
            if (interpolation_kind == (int)NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFMultiquadricManual){
              continue;
            }
            //if (interpolation_kind != (int)NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFTPS){
            //  continue;
            //}
            //if (phase_shifts_count != 4){
            //  continue;
            //}
            //if (filtering_iter_numbers != 2){
            //  continue;
            //}
            //if (phase_shiting_kind != 0 || phase_shifts_count != 5 || interpolation_kind != (int)NP_DSP::ONE_D::FILTERS::InterpolationKind::RBFTPS || full_pi_flag != 0 || filtering_iter_numbers != 1){
            //  continue;
            //}
            auto extractor =
                std::make_shared<NP_DSP::ONE_D::MODES_EXTRACTORS::SOTAEMD>();
            extractor->setPhaseShiftsByCount(
                phase_shifts_count,
                (std::numbers::pi + std::numbers::pi * full_pi_flag) / 2.0,
                !full_pi_flag);
            extractor->max_iter_number_for_filter = filtering_iter_numbers;
            extractor->interpolation_kind =
                static_cast<NP_DSP::ONE_D::FILTERS::InterpolationKind>(
                    interpolation_kind);
            extractor->oversampling_ratio_for_ft_der = 10.0;
            extractor->extremums_rotation_kind_e =
                static_cast<NP_DSP::ONE_D::PHASE_SHIFTERS::RotateKind>(
                    phase_shiting_kind);

            std::stringstream name_stream;
            name_stream << magic_enum::enum_name(extractor->interpolation_kind)
                        << "_iters" << extractor->max_iter_number_for_filter
                        << "_phase_shifts" << extractor->phase_shifts.size()
                        << "_"
                        << magic_enum::enum_name(
                               extractor->extremums_rotation_kind_e)
                        << "_full_pi_is_" << full_pi_flag;
            std::string name = name_stream.str();
            extractors.insert({name, extractor});
          }
        }
      }
    }
  }

  std::vector<std::vector<double>> stationary_modes{};
  for (int i = 0; i < 10; i++) {
    stationary_modes.push_back({});
    for (int j = 0; j < SIZE; j++) {
      stationary_modes[i].push_back(0.0);
    }
    //IC(stationary_modes[i].size(), SIZE);
    addRandomTrigMode(stationary_modes[i], 10, 1000);
    //matplot::plot(stationary_modes[i]);
    //matplot::show();
  }

  std::vector<std::vector<double>> nonstationary_ampl_modes{};
  for (int i = 0; i < 10; i++) {
    nonstationary_ampl_modes.push_back({});
    for (int j = 0; j < SIZE; j++) {
      nonstationary_ampl_modes[i].push_back(0.0);
    }
    addRandomTrigModeWithNonStationaryAmpl(nonstationary_ampl_modes[i], 10,
                                           1000, 10);
    //matplot::plot(nonstationary_ampl_modes[i]);
    //matplot::show();
  }

  std::vector<std::vector<double>> nonstationary_freq_modes{};
  for (int i = 0; i < 10; i++) {
    nonstationary_freq_modes.push_back({});
    for (int j = 0; j < SIZE; j++) {
      nonstationary_freq_modes[i].push_back(0.0);
    }
    addRandomTrigModeWithNonStationaryFreq(nonstationary_freq_modes[i], 10,
                                           1000, 10);
    //matplot::plot(nonstationary_freq_modes);
    //matplot::show();
  }

  std::vector<std::vector<double>> nonstationary_modes{};
  for (int i = 0; i < 10; i++) {
    nonstationary_modes.push_back({});
    for (int j = 0; j < SIZE; j++) {
      nonstationary_modes[i].push_back(0.0);
    }
    addRandomTrigModeWithNonStationaryAmplAndFreq(nonstationary_modes[i], 10, 1000,
                                           10, 10);
    //matplot::plot(nonstationary_modes[i]);
    //matplot::show();
  }
  std::vector<std::vector<double>> random_modes{};
  for (int i = 0; i < 10; i++) {
    random_modes.push_back({});
    for (int j = 0; j < SIZE; j++) {
      random_modes[i].push_back(random() % 1000);
    }
    //matplot::plot(nonstationary_modes[i]);
    //matplot::show();
  }


  std::vector<std::pair<std::string, std::vector<double>>> results;

  int counter = 0;
  for (auto &[key, extractor] : extractors) {
    std::cout << key;
    std::cout << counter << "/" << extractors.size() << std::endl;
    counter++;
    std::vector<double> errors;
    std::pair<std::string, std::vector<double>> res = {key, {}};

    std::vector<std::shared_ptr<std::vector<std::vector<double>>>> matrixes_stationary;
    std::vector<std::shared_ptr<std::vector<std::vector<double>>>> matrixes_nonstationary;
    std::vector<std::shared_ptr<std::vector<std::vector<double>>>> matrixes_nonstationary_freq;
    std::vector<std::shared_ptr<std::vector<std::vector<double>>>> matrixes_nonstationary_ampl;
    std::vector<std::shared_ptr<std::vector<std::vector<double>>>> matrixes_random_modes;
    Spectrogramm spectrogramm;

    auto fillMatrixesFromModes = [&](auto const & modes, auto & matrixes){
      for (auto &mode : modes) {
        //IC(mode.size());
        //std::cout << "compute stationary_matrix" << std::endl;
        extractor->compute(mode);
        IC(extractor->getModesCount());

        spectrogramm.setAxis(SIZE, SIZE);
        spectrogramm.loadFromExtractor(*extractor);
        spectrogramm.computeMatrix();
        //spectrogramm.smoothMatrix(smoothing_radius, 1);
        spectrogramm.smoothMatrixFast(smoothing_sigma, 10);
        std::shared_ptr<std::vector<std::vector<double>>> matrix = std::make_shared<std::vector<std::vector<double>>>();
        copyMatrix(spectrogramm.getMatrix(), *matrix);
        matrixes.push_back(matrix);
        //matplot::plot(mode);
        //matplot::show();
        //spectrogramm.plot();
      }
    };

    std::cout << "compute stationary_matrix" << std::endl;
    fillMatrixesFromModes(stationary_modes, matrixes_stationary);
    std::cout << "compute nonstationary_matrix" << std::endl;
    fillMatrixesFromModes(nonstationary_modes, matrixes_nonstationary);
    std::cout << "compute nonstationary_ampl_matrix" << std::endl;
    fillMatrixesFromModes(nonstationary_ampl_modes, matrixes_nonstationary_ampl);
    std::cout << "compute nonstationary_freq_matrix" << std::endl;
    fillMatrixesFromModes(nonstationary_freq_modes, matrixes_nonstationary_freq);
    std::cout << "compute random_matrix" << std::endl;
    fillMatrixesFromModes(random_modes, matrixes_random_modes);

    std::vector<double> data;
    std::shared_ptr<std::vector<std::vector<double>>> matrix_temp_linear = std::make_shared<std::vector<std::vector<double>>>();

    auto computeErrors = [&](auto const & modes, auto & matrixes){
      data = {};
      for (int i = 0; i < SIZE; i++){
        data.push_back(modes[0][i]);
      }
      copyMatrix(*(matrixes[0]), *matrix_temp_linear);
      for (int i = 1; i < modes.size(); i++){
        std::cout << "mode_linearty iter " << i << std::endl;
        for (int idx = 0; idx < matrix_temp_linear->size(); idx++){
          for (int jdx = 0; jdx < (*matrix_temp_linear)[idx].size(); jdx++){
            (*matrix_temp_linear)[idx][jdx] = (*matrix_temp_linear)[idx][jdx] + (*matrixes[i])[idx][jdx];
          }
        }
        for (int idx = 0; idx < data.size(); idx++){
          data[idx] = data[idx] + modes[i][idx];
        }
        //matplot::plot(data);
        //matplot::show();
        //extractor->debug = true;
        extractor->compute(data);
        //extractor->data.show(NP_DSP::ONE_D::PlottingKind::Simple);
        IC(extractor->getModesCount());
        spectrogramm.setAxis(SIZE, SIZE);
        spectrogramm.loadFromExtractor(*extractor);
        spectrogramm.computeMatrix();
        //spectrogramm.smoothMatrix(smoothing_radius, 1);
        spectrogramm.smoothMatrixFast(smoothing_sigma, 10);
        double error = spectrogramm.computeRMSDistance(*matrix_temp_linear);
        errors.push_back(error);

        std::vector<std::vector<double>> matrix_transposed = {};

        int x_size = matrix_temp_linear->size();
        int y_size = (*matrix_temp_linear)[0].size();
        for (int i = 0; i < y_size; i++) {
          matrix_transposed.push_back({});
          for (int j = 0; j < x_size; j++) {
            matrix_transposed[i].push_back((*matrix_temp_linear)[j][i]);
          }
        }
        IC(error);
        //matplot::plot(data);
        //matplot::show();
        //matplot::image(matrix_transposed, true);
        //matplot::colorbar();
        //matplot::show();
        //spectrogramm.plot();
      }
    };

    std::cout << "process stationary_modes" << std::endl;
    computeErrors(stationary_modes, matrixes_stationary);
    std::cout << "process nonstationary_modes" << std::endl;
    computeErrors(nonstationary_modes, matrixes_nonstationary);
    std::cout << "process nonstationary_ampl_modes" << std::endl;
    computeErrors(nonstationary_ampl_modes, matrixes_nonstationary_ampl);
    std::cout << "process nonstationary_freq_modes" << std::endl;
    computeErrors(nonstationary_freq_modes, matrixes_nonstationary_freq);
    std::cout << "process random_modes" << std::endl;
    computeErrors(random_modes, matrixes_random_modes);

    data = {};
    for (int i = 0; i < SIZE; i++){
      data.push_back(random_modes[0][i]);
    }
    copyMatrix(*(matrixes_random_modes[0]), *matrix_temp_linear);
    for (int i = 1; i < random_modes.size(); i++){
      std::cout << "compute nonstationary_freq_mode_linearty iter " << i << std::endl;
      for (int idx = 0; idx < matrix_temp_linear->size(); idx++){
        for (int jdx = 0; jdx < (*matrix_temp_linear)[idx].size(); jdx++){
          (*matrix_temp_linear)[idx][jdx] = (*matrix_temp_linear)[idx][jdx]
            + (*matrixes_random_modes[i])[idx][jdx];
        }
      }
      for (int idx = 0; idx < data.size(); idx++){
        data[idx] = data[idx] + random_modes[i][idx];
      }
      //matplot::plot(data);
      //matplot::show();
      extractor->compute(data);
      IC(extractor->getModesCount());
      spectrogramm.setAxis(SIZE, SIZE);
      spectrogramm.loadFromExtractor(*extractor);
      spectrogramm.computeMatrix();
      //spectrogramm.smoothMatrix(smoothing_radius, 1);
      spectrogramm.smoothMatrixFast(smoothing_sigma, 10);
      double error = spectrogramm.computeL2Distance(*matrix_temp_linear);
      errors.push_back(error);
    }
    res.second = errors;
    results.push_back(res);
    //break;
  }


  std::vector<double> average_errors;
  for (int i = 0; i < results[0].second.size(); i++){
    average_errors.push_back(0.0);
  }
  for (auto & [key, errors]: results){
    for (int i = 0; i < errors.size(); i++){
      average_errors[i] += errors[i];
    }
  }

  for (int i = 0; i < average_errors.size(); i++){
    average_errors[i] = average_errors[i] / results.size();
  }

  for (auto & [key, errors]: results){
    for (int i = 0; i < errors.size(); i++){
      errors[i] = errors[i]/average_errors[i];
    }
  }

  auto comparator = [&] (std::pair<std::string, std::vector<double>> a,
    std::pair<std::string, std::vector<double>> b){
    double L2_1 = 0.0;
    double L2_2 = 0.0;

    for(int i = 0; i < a.second.size(); i++){
      L2_1 += a.second[i]*a.second[i];
      L2_2 += b.second[i]*b.second[i];
    }
    //L2_1 = sqrt(L2_1);
    //L2_2 = sqrt(L2_2);

    return L2_1 < L2_2;
  };

  std::sort(results.begin(), results.end(), comparator);

  nlohmann::json report;

  for (int i = 0; i < results.size(); i++){
    nlohmann::json report_extractor;
    report_extractor["name"] = results[i].first;
    nlohmann::json errors;

    //std::cout << results[i].first << " ";
    
    for (int j = 0; j < results[i].second.size(); j++){
      //std::cout << results[i].second[j] << " ";
      errors.push_back(results[i].second[j]);
    }
    report_extractor["errors"] = errors;
    report.push_back(report_extractor);
    //std::cout << std::endl;
  }
  std::ofstream o("/home/dmitry/projects/non-parametric-dsp/experimants_results3.json");
  o << std::setw(4) << report << std::endl;

  return 0;
}
