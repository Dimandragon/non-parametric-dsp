#pragma once

//#include <icecream.hpp>

#include <derivators.hpp>
#include <filters.hpp>
#include <inst_ampl_computers.hpp>
#include <inst_freq_computers.hpp>
#include <integrators.hpp>
#include <npdsp_concepts.hpp>
#include <phase_computers.hpp>
#include <signals.hpp>
#include <vector>

namespace NP_DSP::ONE_D::MODES_EXTRACTORS {

enum class InstAmplKind { HTBased, ExtremumsBased };

struct SOTAEMD {
  NP_DSP::ONE_D::INST_FREQ_COMPUTERS::InstFreqComputingKind inst_freq_computing_kind = 
    INST_FREQ_COMPUTERS::InstFreqComputingKind::ExtremumsBasedMultiquadric;

  double max_iter_number_for_filter = 10;
  NP_DSP::ONE_D::PHASE_SHIFTERS::RotateKind extremums_rotation_kind_e =
      NP_DSP::ONE_D::PHASE_SHIFTERS::RotateKind::Naive;
  double oversampling_ratio_for_ft_der = 1.0;

  std::vector<double> phase_shifts{
      0.0 * std::numbers::pi, 0.1 * std::numbers::pi, 0.2 * std::numbers::pi,
      0.3 * std::numbers::pi, 0.4 * std::numbers::pi, 0.5 * std::numbers::pi,
      0.6 * std::numbers::pi, 0.7 * std::numbers::pi, 0.8 * std::numbers::pi,
      0.9 * std::numbers::pi};

  int idw_layers = 15;
  double idw_search_radius = 100;

  double rbf_r_base = 20;
  double rbf_n_layers = 20;
  double rbf_lambda_n_s = 0.0;
  double rbf_search_r = 0.4;
  bool rbf_v3tol = true;

  double rbf_lambda_v = 0.0;

  double rbf_alpha = 10.0;
  NP_DSP::ONE_D::FILTERS::InterpolationKind interpolation_kind_e = NP_DSP::ONE_D::FILTERS::InterpolationKind::Makima;

  using DataType = GenericSignal<SimpleVecWrapper<double>, true>;
  DataType data;
  DataType data_buffer;
  DataType compute_buffer;
  DataType compute_buffer2;
  std::vector<DataType *> modes;
  std::vector<DataType *> inst_freqs;
  std::vector<DataType *> inst_ampls;
  std::vector<DataType *> phases;
  std::vector<double> freq_conv;
  std::vector<double> freq_conv_image;

  bool debug = false;

  INTEGRATORS::Riman<INTEGRATORS::PolygonType::ByPoint> integrator;
  DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward>
      derivator;
  PHASE_COMPUTERS::ExtremumsBasedNonOpt<
      double, PHASE_COMPUTERS::ExtremumsKind::DerArctg, decltype(derivator)>
      phase_computer_der_atan;
  PHASE_COMPUTERS::ExtremumsBasedNonOpt<
      double, PHASE_COMPUTERS::ExtremumsKind::Simple, decltype(derivator)>
      phase_computer_simple;

  INST_FREQ_COMPUTERS::SOTAInstFreqComputer inst_freq_computer;

  INST_AMPL_COMPUTERS::HilbertTransformBased<UTILITY_MATH::HTKind::Mull>
      inst_ampl_computer;

  FILTERS::RecursiveFilter<double,
                           FILTERS::LocalFilteringType::InterpolationExtremums>
      filter;

  template <typename DataInT> void load(const DataInT &data_in) {
    data.base->vec->clear();
    for (auto i = 0; i < data_in.size(); i++) {
      data.base->vec->push_back(data_in[i]);
    }
    for (int i = 0; i < modes.size(); i++) {
      modes[i]->base->vec->clear();
    }
    for (int i = 0; i < inst_freqs.size(); i++) {
      inst_freqs[i]->base->vec->clear();
    }
    for (int i = 0; i < phases.size(); i++) {
      phases[i]->base->vec->clear();
    }
    for (int i = 0; i < inst_ampls.size(); i++) {
      inst_ampls[i]->base->vec->clear();
    }
  }

  void setPhaseShiftsByCount(int count, double max_shift, bool including_max_shift){
    if (count == 0){
      phase_shifts = {0.0};
    }
    else{
      double d_phase_shift = max_shift / (count);
      if (including_max_shift){
        if (count != 1){
          d_phase_shift = max_shift / (count + 1);
        }
      }
      phase_shifts.clear();
      for (int i = 0; i < count; i++){
        phase_shifts.push_back(i * d_phase_shift);
      }
    }
  }

  template <typename DataT> void compute(const DataT &data_in) {
    inst_freq_computer.kind = inst_freq_computing_kind;
    filter.filter.idw_layers = idw_layers;
    filter.filter.idw_search_radius = idw_search_radius;
    filter.filter.rbf_r_base = rbf_r_base;
    filter.filter.rbf_n_layers = rbf_n_layers;
    filter.filter.rbf_lambda_n_s = rbf_lambda_n_s;
    filter.filter.rbf_search_r = rbf_search_r;
    filter.filter.rbf_v3tol = rbf_v3tol;
    filter.filter.rbf_lambda_v = rbf_lambda_v;
    filter.filter.rbf_alpha = rbf_alpha;
    filter.filter.interpolation_kind = interpolation_kind_e;

    filter.filter.phase_shifts = phase_shifts;

    auto phase_shifts_local = phase_shifts;
    filter.max_iters = max_iter_number_for_filter;
    filter.filter.extremums_rotation_kind_e = extremums_rotation_kind_e;
    filter.filter.oversampling_ratio_for_ft_der = oversampling_ratio_for_ft_der;
    filter.debug = debug;
    filter.filter.debug = debug;

    size_t iter_number = 0;

    DataType non_resampled_data;

    

    auto prepare_memory_ext = [&]() {
      for (int i = 0; i < data_in.size(); i++) {
        non_resampled_data.base->vec->push_back(data_in[i]);
      }
      if (data.size() != data_in.size()) {
        data.base->vec->clear();
        for (int i = 0; i < data_in.size(); i++) {
          data.base->vec->push_back(data_in[i]);
        }
      }
      else{
        for (int i = 0; i < data_in.size(); i++) {
          data[i] = data_in[i];
        }
      }
      if (data_buffer.size() != data_in.size()) {
        data_buffer.base->vec->clear();
        for (int i = 0; i < data_in.size(); i++) {
          data_buffer.base->vec->push_back(data_in[i]);
        }
      }
      if (compute_buffer.size() != data_in.size()) {
        compute_buffer.base->vec->clear();
        for (int i = 0; i < data_in.size(); i++) {
          compute_buffer.base->vec->push_back(0.0);
        }
      }
      if (compute_buffer2.size() != data_in.size()) {
        compute_buffer2.base->vec->clear();
        for (int i = 0; i < data_in.size(); i++) {
          compute_buffer2.base->vec->push_back(0.0);
        }
      }
      if (freq_conv.size() != data.size()) {
        freq_conv.clear();
        for (int i = 0; i < data.size(); i++) {
          freq_conv.push_back(1.0);
        }
      }
      if (freq_conv_image.size() != data.size()) {
        freq_conv_image.clear();
        for (int i = 0; i < data.size(); i++) {
          freq_conv_image.push_back(1.0);
        }
      }
    };

    prepare_memory_ext();

    auto prepare_memory_int = [&]() {
      if (iter_number + 1 > modes.size()) {
        auto *modes_new = new GenericSignal<SimpleVecWrapper<double>, true>;
        modes.push_back(modes_new);
        for (int i = 0; i < data.size(); i++) {
          modes[iter_number]->base->vec->push_back(0.0);
        }
      }
      if (iter_number + 1 > inst_freqs.size()) {
        auto *inst_freqs_new =
            new GenericSignal<SimpleVecWrapper<double>, true>;
        inst_freqs.push_back(inst_freqs_new);
        for (int i = 0; i < data.size(); i++) {
          inst_freqs[iter_number]->base->vec->push_back(0.0);
        }
      }
      if (iter_number + 1 > inst_ampls.size()) {
        auto *inst_ampls_new =
            new GenericSignal<SimpleVecWrapper<double>, true>;
        inst_ampls.push_back(inst_ampls_new);
        for (int i = 0; i < data.size(); i++) {
          inst_ampls[iter_number]->base->vec->push_back(0.0);
        }
      }
      if (iter_number + 1 > phases.size()) {
        auto *phase_new = new GenericSignal<SimpleVecWrapper<double>, true>;
        phases.push_back(phase_new);
        for (int i = 0; i < data.size(); i++) {
          phases[iter_number]->base->vec->push_back(0.0);
        }
      }
      if (modes[iter_number]->size() != data.size()) {
        modes[iter_number]->base->vec->clear();
        for (int i = 0; i < data.size(); i++) {
          modes[iter_number]->base->vec->push_back(0.0);
        }
      }
      if (inst_freqs[iter_number]->size() != data.size()) {
        inst_freqs[iter_number]->base->vec->clear();
        for (int i = 0; i < data.size(); i++) {
          inst_freqs[iter_number]->base->vec->push_back(0.0);
        }
      }
      if (inst_ampls[iter_number]->size() != data.size()) {
        inst_ampls[iter_number]->base->vec->clear();
        for (int i = 0; i < data.size(); i++) {
          inst_ampls[iter_number]->base->vec->push_back(0.0);
        }
      }
      if (phases[iter_number]->size() != data.size()) {
        phases[iter_number]->base->vec->clear();
        for (int i = 0; i < data.size(); i++) {
          phases[iter_number]->base->vec->push_back(0.0);
        }
      }
    };

    while (true) {
      prepare_memory_int();
      phase_computer_simple.compute(data, *phases[iter_number], nullptr);

      if ((*phases[iter_number])[data.size() - 1] > 6.28) {
        filter.compute(data, data_buffer, &compute_buffer);
        for (int i = 0; i < data.size(); i++) {
          (*modes[iter_number])[i] = data[i] - data_buffer[i];
          data[i] = data_buffer[i]; // data is filtered signal
                                    // data_buffer is mode
        }
        phase_computer_simple.compute(*modes[iter_number], *phases[iter_number],
                                      nullptr);
        inst_freq_computer.compute(*modes[iter_number],
                                   *inst_freqs[iter_number], nullptr);
        inst_ampl_computer.compute(*modes[iter_number],
                                   *inst_ampls[iter_number], nullptr);
        iter_number++;
      } else {
        for (auto i = 0; i < data.size(); i++) {
          (*modes[iter_number])[i] = data[i];
        }
        phase_computer_simple.compute(*modes[iter_number], *phases[iter_number],
                                      nullptr);
        inst_freq_computer.compute(*modes[iter_number],
                                   *inst_freqs[iter_number], nullptr);
        inst_ampl_computer.compute(*modes[iter_number],
                                   *inst_ampls[iter_number], nullptr);

        for(int i = iter_number + 1; i < modes.size(); i++){
          delete modes[i];
        }
        for(int i = iter_number + 1; i < phases.size(); i++){
          delete phases[i];
        }
        for(int i = iter_number + 1; i < inst_freqs.size(); i++){
          delete inst_freqs[i];
        }
        for(int i = iter_number + 1; i < inst_ampls.size(); i++){
          delete inst_ampls[i];
        }
        while (modes.size() - 1 > iter_number){
          modes.pop_back();
        }
        while (phases.size() - 1 > iter_number){
          phases.pop_back();
        }
        while (inst_ampls.size() - 1 > iter_number){
          inst_ampls.pop_back();
        }
        while (inst_freqs.size() - 1 > iter_number){
          inst_freqs.pop_back();
        }
        iter_number++;
        return;
      }
    }
  }

  int getModesCount() const { return static_cast<int>(modes.size()); }

  int getDataSize() const { return static_cast<int>(modes[0]->size()); }
  std::vector<double> getMode(int idx) const {
    return *(modes[idx]->base->vec);
  }
  std::vector<double> getInstFreq(int idx) const {
    return *(inst_freqs[idx]->base->vec);
  }
  std::vector<double> getInstAmpl(int idx) const {
    return *(inst_ampls[idx]->base->vec);
  }
  std::vector<double> getPhase(int idx) const {
    return *(phases[idx]->base->vec);
  }

  ~SOTAEMD(){
    for (auto pointer: modes){
      delete pointer;
    }
    modes.clear();
    for (auto pointer: inst_freqs){
      delete pointer;
    }
    inst_freqs.clear();
    for (auto pointer: inst_ampls){
      delete pointer;
    }
    inst_ampls.clear();
    for (auto pointer: phases){
      delete pointer;
    }
    phases.clear();
  }
};
} // namespace NP_DSP::ONE_D::MODES_EXTRACTORS
