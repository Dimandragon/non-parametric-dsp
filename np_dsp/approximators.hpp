#pragma once

#include <cstdint>
#include <icecream.hpp>

#include <npdsp_concepts.hpp>
#include <npdsp_config.hpp>
#include <vector>
#include <complex>
#include <utility_math.hpp>
#include <numbers>
#include <string>
#include <matplot/matplot.h>
#include <makima.hpp>
#include <pchip.hpp>
#include <optional>

namespace NP_DSP::ONE_D::APPROX {
     enum class FSApproxKind { Simple, Positive };

    
    template<typename LossFunc, typename StopPointFunc, FSApproxKind kind_v, typename BySampleLoss>
    struct FourierSeriesBased {
        using Loss = LossFunc;
        using StopPoint = StopPointFunc;
        static constexpr bool is_signal_approximator = true;

        static constexpr FSApproxKind kind = kind_v;

        Loss* loss;
        BySampleLoss* bySampleLoss = nullptr;

        StopPoint* stopPont;

        size_t signal_size;
        bool is_actual = false;
        int polynoms_count;

        int polynoms_count_on_tile;
        int tile_size;

        std::vector<std::complex<double>> fourier_series;
        std::vector<std::complex<double>> approximated_data;

        double max_value = 10000000000000000000000000.;

        template<Signal SignalType>
        FourierSeriesBased(Loss& lossFn, SignalType& signal_in, StopPointFunc& stop_point) {
            loss = &lossFn;
            signal_size = signal_in.size();
            fourier_series = std::vector<std::complex<double>>(signal_in.size());
            approximated_data = std::vector<std::complex<double>>(signal_in.size());
            stopPont = &stop_point;
            polynoms_count = signal_in.size() / 2;
            polynoms_count_on_tile = signal_in.size() / 2;
            tile_size = signal_in.size();
        }

        template<Signal SignalType>
        void computeFSFromData(SignalType& signal_in){
            for (size_t i = 0; i < signal_in.size(); i++){
                approximated_data[i] = {signal_in[i], 0.0};
            }
            is_actual = true;
            for (auto i = 0; i < approximated_data.size() / tile_size; i++) {
                auto const pad = i * tile_size;
                UTILITY_MATH::fftc2c(approximated_data, fourier_series, tile_size, pad);
            }
            auto const pad = approximated_data.size() / tile_size * tile_size;
            UTILITY_MATH::fftc2c(approximated_data, fourier_series, approximated_data.size() - pad, pad);
        }

        size_t mirrorIdx(size_t const i) const {
            return i / tile_size * tile_size + i / tile_size * tile_size + tile_size - i;
        }

        void applyMirror(size_t const i) {
            if (i % tile_size != 0) {
                auto const mirror_idx = mirrorIdx(i);
                fourier_series[mirror_idx] = {fourier_series[i].real() / 2.0, -fourier_series[i].imag() / 2.0};
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
            UTILITY_MATH::ifftc2c(fourier_series, approximated_data, approximated_data.size() - pad, pad);
        }

        void computeTile(size_t const idx) {
            auto const i = idx / tile_size;
            auto const pad = i * tile_size;
            if (approximated_data.size() > pad + tile_size) {
                UTILITY_MATH::ifftc2c(fourier_series, approximated_data, tile_size, pad);
            } else {
                UTILITY_MATH::ifftc2c(fourier_series, approximated_data, approximated_data.size() - pad, pad);
            }
        }

        void setApproxOrderRatio(double const ratio) {
            polynoms_count_on_tile = static_cast<int>(static_cast<double>(tile_size) * 0.5 * ratio);
        }

        template<typename IdxT>
        std::complex<double> computeSample(IdxT idx) {
            using SampleType = double;
            std::complex<double> accum = {0.0, 0.0};

            if (idx <= approximated_data.size() && idx >= 0) {
                auto tile_first_idx = static_cast<int>(idx) / tile_size * tile_size;
                if (tile_first_idx + tile_size < approximated_data.size()) {
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) / static_cast<SampleType>(
                                       tile_size);
                    for (auto i = 0; i < tile_size; i++) {
                        accum += fourier_series[i + tile_first_idx] * std::complex<SampleType>{
                            std::cos(w * i), std::sin(w * i)
                        };
                    }
                } else {
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                                   static_cast<SampleType>(approximated_data.size() - tile_first_idx);
                    for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
                        accum += fourier_series[i + tile_first_idx] * std::complex<SampleType>{
                            std::cos(w * i), std::sin(w * i)
                        };
                    }
                }
            } else if (idx > approximated_data.size()) {
                auto tile_first_idx = approximated_data.size() / tile_size * tile_size;
                SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                               static_cast<SampleType>(approximated_data.size() - tile_first_idx);
                for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
                    accum += fourier_series[i + tile_first_idx] * std::complex<SampleType>{
                        std::cos(w * i), std::sin(w * i)
                    };
                }
            } else {
                if (approximated_data.size() > tile_size) {
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                                   static_cast<SampleType>(tile_size);
                    for (auto i = 0; i < tile_size; i++) {
                        accum += fourier_series[i] * std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
                    }
                } else {
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                                   static_cast<SampleType>(approximated_data.size());
                    for (auto i = 0; i < approximated_data.size(); i++) {
                        accum += fourier_series[i] * std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
                    }
                }
            }

            return accum;
        }

        template<typename SampleType, typename IdxType>
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

        template<typename SampleType, typename IdxType>
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
            //todo for tiles and for positive
            for (auto i = 0; i < fourier_series.size(); i++) {
                if (i % tile_size > tile_size / 2) {
                    //auto idx_mirror = mirrorIdx(i);
                    //fourier_series[i] = {fourier_series[idx_mirror].real() / 2.0, -fourier_series[idx_mirror].imag() / 2.0};
                    //fourier_series[idx_mirror] = fourier_series[idx_mirror] / 2.0;
                    continue;
                }
                if (i % tile_size >= polynoms_count_on_tile) {
                    continue;
                }


                auto check_loss = [&](std::pair<SampleType, SampleType> data) {
                    auto const complex_sample = UTILITY_MATH::convertFSampleT2C<SampleType>(data);
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

                auto trigonometric_sample = UTILITY_MATH::convertFSampleC2T(fourier_series[i]);

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
                auto central_loss = check_loss({static_cast<SampleType>(1.), theta_central});
                auto right_loss = check_loss({static_cast<SampleType>(1.), theta_max});

                if constexpr (CONFIG::debug) {
                    std::string mark = "end compute first losses";
                    IC(mark);

                    IC(theta_min, left_loss, theta_central, central_loss, theta_max, right_loss);
                }

                auto period_opt_iter = 0;

                if constexpr (CONFIG::debug) {
                    std::string mark = "start optimize period";
                    IC(mark);
                }


                while (!(*stopPont)(std::abs(right_loss - central_loss) + std::abs(left_loss - central_loss), *this)) {
                    if constexpr (CONFIG::debug) {
                        IC(period_opt_iter, left_loss, central_loss, right_loss, theta_min, theta_central, theta_max);
                        IC(right_loss-central_loss, left_loss-central_loss, right_loss-left_loss);
                        IC(theta_central - theta_min, theta_max-theta_central);
                    }

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
                        //left branch is higher
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
                        //right branch is higher
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


                if constexpr (CONFIG::debug) {
                    std::string mark = "start optimize ampl";
                    IC(mark);
                }
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
                //while(!(*stopPont)(std::abs(ampl_right - ampl_left), *this)) {
                auto errors_counter = 0;
                while (!(*stopPont)(std::abs(loss_right - loss_central) + std::abs(loss_left - loss_central), *this)) {
                    if constexpr (CONFIG::debug) {
                        IC(ampl_opt_iter, loss_left, loss_central, loss_right, ampl_left, ampl_central, ampl_right);
                        IC(loss_right-loss_central, loss_left-loss_central, loss_right-loss_left);
                        IC(ampl_central - ampl_left, ampl_right - ampl_central);
                    }

                    auto left_diff = loss_left - loss_central;
                    auto right_diff = loss_right - loss_central;
                    if (left_diff > 0 && right_diff <= 0) {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "ampl 1 branch";
                            IC(point_mark);
                        }
                        loss_left = loss_central;
                        ampl_left = ampl_central;
                        ampl_central = (ampl_left + ampl_right) / 2;
                        trigonometric_sample.first = ampl_central;
                        loss_central = check_loss(trigonometric_sample);
                    } else if (left_diff <= 0 && right_diff > 0) {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "ampl 2 branch";
                            IC(point_mark);
                        }
                        loss_right = loss_central;
                        ampl_right = ampl_central;
                        ampl_central = (ampl_left + ampl_right) / 2;
                        trigonometric_sample.first = ampl_central;
                        loss_central = check_loss(trigonometric_sample);
                    } else if (left_diff > right_diff) {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "ampl 3 branch";
                            IC(point_mark);
                        }
                        //left branch is higher
                        auto ampl_left_avg = (ampl_left + ampl_central) / 2;
                        trigonometric_sample.first = ampl_left_avg;
                        auto new_loss_left = check_loss(trigonometric_sample);
                        if (new_loss_left < loss_left && new_loss_left > central_loss) {
                            if constexpr (CONFIG::debug) {
                                std::string point_mark = "ampl 3.1 branch";
                                IC(point_mark);
                            }
                            loss_left = new_loss_left;
                            ampl_left = ampl_left_avg;
                            ampl_central = (ampl_left + ampl_right) / 2;
                            trigonometric_sample.first = ampl_central;
                            loss_central = check_loss(trigonometric_sample);
                        } else {
                            if constexpr (CONFIG::debug) {
                                std::string point_mark = "ampl 3.2 branch";
                                IC(point_mark);
                            }
                            ampl_right = (ampl_central + ampl_right) / 2;
                            trigonometric_sample.first = ampl_right;
                            loss_right = check_loss(trigonometric_sample);
                            ampl_central = (ampl_right + ampl_left) / 2;
                            trigonometric_sample.first = ampl_central;
                            loss_central = check_loss(trigonometric_sample);
                        }
                    } else if (left_diff <= right_diff) {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "ampl 4 branch";
                            IC(point_mark);
                        }
                        auto ampl_right_avg = (ampl_right + ampl_central) / 2;
                        trigonometric_sample.first = ampl_right_avg;
                        auto new_loss_right = check_loss(trigonometric_sample);
                        if (new_loss_right < loss_right && new_loss_right > loss_central) {
                            if constexpr (CONFIG::debug) {
                                std::string point_mark = "ampl 4.1 branch";
                                IC(point_mark);
                            }
                            loss_right = new_loss_right;
                            ampl_right = ampl_right_avg;
                            ampl_central = (ampl_left + ampl_right) / 2;
                            trigonometric_sample.first = ampl_central;
                            loss_central = check_loss(trigonometric_sample);
                        } else {
                            if constexpr (CONFIG::debug) {
                                std::string point_mark = "ampl 4.2 branch";
                                IC(point_mark);
                            }
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
                            if constexpr (CONFIG::debug) {
                                IC(left_diff, loss_left-loss_central, right_diff, loss_right-loss_central);
                            }

                            errors_counter++;
                            if constexpr (CONFIG::debug) {
                                IC(errors_counter);
                            }

                            if (errors_counter > 10) {
                                if constexpr (CONFIG::debug) {
                                    for (;;) {
                                    }
                                }
                            }
                        } else {
                            if constexpr (CONFIG::debug) {
                                IC(right_diff == loss_right-loss_central);
                                IC(left_diff, loss_left-loss_central, right_diff, loss_right-loss_central);
                            }
                        }
                    }
                    ampl_opt_iter++;
                    if (ampl_central - ampl_left < 0.000001 || ampl_right - ampl_central < 0.000001) {
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
                if constexpr (CONFIG::debug) {
                    std::string mark = "comopute fs cf";
                    IC(mark);
                    IC(i);
                }

                if constexpr (CONFIG::debug) {
                    std::string mark = "compute first losses";
                    IC(mark);
                }

                auto check_loss = [&](std::pair<SampleType, SampleType> data) {
                    auto const complex_sample = UTILITY_MATH::convertFSampleT2C<SampleType>(data);
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
                auto left_loss = check_loss({static_cast<SampleType>(max_ampl), theta_min});
                auto right_loss = check_loss({static_cast<SampleType>(max_ampl), theta_max});
                auto central_loss = check_loss({static_cast<SampleType>(max_ampl), theta_central});

                if constexpr (CONFIG::debug) {
                    std::string mark = "end compute first losses";
                    IC(mark);

                    IC(theta_min, left_loss, theta_central, central_loss, theta_max, right_loss);
                }

                auto period_opt_iter = 0;

                if constexpr (CONFIG::debug) {
                    std::string mark = "start optimize period";
                    IC(mark);
                }


                while (!(*stopPont)(std::abs(right_loss - central_loss) + std::abs(left_loss - central_loss),
                                     *this)) {
                    if constexpr (CONFIG::debug) {
                        IC(period_opt_iter, left_loss, central_loss, right_loss, theta_min, theta_central, theta_max);
                        IC(right_loss-central_loss, left_loss-central_loss, right_loss-left_loss);
                        IC(theta_central - theta_min, theta_max-theta_central);
                    }

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
                        //left branch is higher
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
                        //right branch is higher
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


                if constexpr (CONFIG::debug) {
                    std::string mark = "start optimize ampl";
                    IC(mark);
                }

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
                while (!((*stopPont)(std::abs(loss_right - loss_central) + std::abs(loss_left - loss_central),
                                     *this))) {
                    if constexpr (CONFIG::debug) {
                        IC(ampl_opt_iter, loss_left, loss_central, loss_right, ampl_left, ampl_central, ampl_right);
                        IC(loss_right-loss_central, loss_left-loss_central, loss_right-loss_left);
                        IC(ampl_central - ampl_left, ampl_right - ampl_central);
                    }

                    auto left_diff = loss_left - loss_central;
                    auto right_diff = loss_right - loss_central;
                    if (left_diff > 0 && right_diff <= 0) {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "ampl 1 branch";
                            IC(point_mark);
                        }
                        loss_left = loss_central;
                        ampl_left = ampl_central;
                        ampl_central = (ampl_left + ampl_right) / 2;
                        trigonometric_sample.first = ampl_central;
                        loss_central = check_loss(trigonometric_sample);
                    } else if (left_diff <= 0 && right_diff > 0) {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "ampl 2 branch";
                            IC(point_mark);
                        }
                        loss_right = loss_central;
                        ampl_right = ampl_central;
                        ampl_central = (ampl_left + ampl_right) / 2;
                        trigonometric_sample.first = ampl_central;
                        loss_central = check_loss(trigonometric_sample);
                    } else if (left_diff > right_diff) {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "ampl 3 branch";
                            IC(point_mark);
                        }
                        //left branch is higher
                        auto ampl_left_avg = (ampl_left + ampl_central) / 2;
                        trigonometric_sample.first = ampl_left_avg;
                        auto new_loss_left = check_loss(trigonometric_sample);
                        if (new_loss_left < loss_left && new_loss_left > central_loss) {
                            if constexpr (CONFIG::debug) {
                                std::string point_mark = "ampl 3.1 branch";
                                IC(point_mark);
                            }
                            loss_left = new_loss_left;
                            ampl_left = ampl_left_avg;
                            ampl_central = (ampl_left + ampl_right) / 2;
                            trigonometric_sample.first = ampl_central;
                            loss_central = check_loss(trigonometric_sample);
                        } else {
                            if constexpr (CONFIG::debug) {
                                std::string point_mark = "ampl 3.2 branch";
                                IC(point_mark);
                            }
                            ampl_right = (ampl_central + ampl_right) / 2;
                            trigonometric_sample.first = ampl_right;
                            loss_right = check_loss(trigonometric_sample);
                            ampl_central = (ampl_right + ampl_left) / 2;
                            trigonometric_sample.first = ampl_central;
                            loss_central = check_loss(trigonometric_sample);
                        }
                    } else if (left_diff <= right_diff) {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "ampl 4 branch";
                            IC(point_mark);
                        }
                        auto ampl_right_avg = (ampl_right + ampl_central) / 2;
                        trigonometric_sample.first = ampl_right_avg;
                        auto new_loss_right = check_loss(trigonometric_sample);
                        if (new_loss_right < loss_right && new_loss_right > loss_central) {
                            if constexpr (CONFIG::debug) {
                                std::string point_mark = "ampl 4.1 branch";
                                IC(point_mark);
                            }
                            loss_right = new_loss_right;
                            ampl_right = ampl_right_avg;
                            ampl_central = (ampl_left + ampl_right) / 2;
                            trigonometric_sample.first = ampl_central;
                            loss_central = check_loss(trigonometric_sample);
                        } else {
                            if constexpr (CONFIG::debug) {
                                std::string point_mark = "ampl 4.2 branch";
                                IC(point_mark);
                            }
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
                            if constexpr (CONFIG::debug) {
                                IC(left_diff, loss_left-loss_central, right_diff, loss_right-loss_central);
                            }

                            errors_counter++;
                            if constexpr (CONFIG::debug) {
                                IC(errors_counter);
                            }

                            if (errors_counter > 10) {
                                if constexpr (CONFIG::debug) {
                                    for (;;) {
                                    }
                                }
                            }
                        } else {
                            if constexpr (CONFIG::debug) {
                                IC(right_diff == (loss_right-loss_central));
                                IC(left_diff, loss_left-loss_central, right_diff, loss_right-loss_central);
                            }
                        }
                    } else {
                        if constexpr (CONFIG::debug) {
                            std::string point_mark = "aaaa fuck";
                            IC(point_mark);
                            IC(left_diff, loss_left-loss_central, right_diff, loss_right-loss_central);
                        }
                    }
                    ampl_opt_iter++;
                    if (ampl_central - ampl_left < 0.000001 || ampl_right - ampl_central < 0.000001) {
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

        void show(const PlottingKind kind) {
            using SampleType = double;
            is_actual = false;
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < signal_size; i++) {
                    plotting_data.push_back(compute<double, size_t>(i));
                }
                matplot::plot(plotting_data);
                matplot::show();
            }
        }

        void show(const PlottingKind kind, const std::string& filename, const std::string& format) {
            is_actual = false;
            using SampleType = double;
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < signal_size; i++) {
                    plotting_data.push_back(compute<double, size_t>(i));
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename, format);
            }
        }

        void show(const PlottingKind kind, const std::string& filename) {
            is_actual = false;
            using SampleType = double;
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < signal_size; i++) {
                    plotting_data.push_back(compute<double, size_t>(i));
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename);
            }
        }
    };


    struct FourierSeriesBasedWithNoTrain {
        //static constexpr bool is_signal_approximator = true;
        size_t signal_size;
        bool is_actual = false;
        int polynoms_count;

        int polynoms_count_on_tile;
        int tile_size;

        std::vector<std::complex<double>> fourier_series;
        std::vector<std::complex<double>> approximated_data;


        template<Signal SignalType>
        FourierSeriesBasedWithNoTrain(const SignalType & signal_in) {
            signal_size = signal_in.size();
            fourier_series = std::vector<std::complex<double>>(signal_in.size());
            approximated_data = std::vector<std::complex<double>>(signal_in.size());
            polynoms_count = signal_in.size();
            polynoms_count_on_tile = signal_in.size();
            tile_size = signal_in.size();
        }

        template<Signal SignalType>
        void computeFSFromData(SignalType& signal_in){
            for (size_t i = 0; i < signal_in.size(); i++){
                approximated_data[i] = {signal_in[i], 0.0};
            }
            is_actual = true;
            for (auto i = 0; i < approximated_data.size() / tile_size; i++) {
                auto const pad = i * tile_size;
                UTILITY_MATH::fftc2c(approximated_data, fourier_series, tile_size, pad);
            }
            auto const pad = approximated_data.size() / tile_size * tile_size;
            UTILITY_MATH::fftc2c(approximated_data, fourier_series, approximated_data.size() - pad, pad);
        }

        size_t mirrorIdx(size_t const i) const {
            return i / tile_size * tile_size + i / tile_size * tile_size + tile_size - i;
        }

        void applyMirror(size_t const i) {
            if (i % tile_size != 0) {
                auto const mirror_idx = mirrorIdx(i);
                fourier_series[mirror_idx] = {fourier_series[i].real() / 2.0, -fourier_series[i].imag() / 2.0};
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
            UTILITY_MATH::ifftc2c(fourier_series, approximated_data, approximated_data.size() - pad, pad);
        }

        void computeTile(size_t const idx) {
            auto const i = idx / tile_size;
            auto const pad = i * tile_size;
            if (approximated_data.size() > pad + tile_size) {
                UTILITY_MATH::ifftc2c(fourier_series, approximated_data, tile_size, pad);
            } else {
                UTILITY_MATH::ifftc2c(fourier_series, approximated_data, approximated_data.size() - pad, pad);
            }
        }

        void setApproxOrderRatio(double const ratio) {
            polynoms_count_on_tile = static_cast<int>(static_cast<double>(tile_size) * 0.5 * ratio);
        }

        template<typename IdxT>
        std::complex<double> computeSample(IdxT idx) {
            using SampleType = double;
            std::complex<double> accum = {0.0, 0.0};

            if (idx <= approximated_data.size() && idx >= 0) {
                auto tile_first_idx = idx / tile_size * tile_size;
                if (tile_first_idx + tile_size < approximated_data.size()) {
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) / static_cast<SampleType>(
                                       tile_size);
                    for (auto i = 0; i < tile_size; i++) {
                        accum += fourier_series[i + tile_first_idx] * std::complex<SampleType>{
                            std::cos(w * i), std::sin(w * i)
                        };
                    }
                } else {
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                                   static_cast<SampleType>(approximated_data.size() - tile_first_idx);
                    for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
                        accum += fourier_series[i + tile_first_idx] * std::complex<SampleType>{
                            std::cos(w * i), std::sin(w * i)
                        };
                    }
                }
            } else if (idx > approximated_data.size()) {
                auto tile_first_idx = approximated_data.size() / tile_size * tile_size;
                SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                               static_cast<SampleType>(approximated_data.size() - tile_first_idx);
                for (auto i = 0; i < approximated_data.size() - tile_first_idx; i++) {
                    accum += fourier_series[i + tile_first_idx] * std::complex<SampleType>{
                        std::cos(w * i), std::sin(w * i)
                    };
                }
            } else {
                if (approximated_data.size() > tile_size) {
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                                   static_cast<SampleType>(tile_size);
                    for (auto i = 0; i < tile_size; i++) {
                        accum += fourier_series[i] * std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
                    }
                } else {
                    SampleType w = std::numbers::pi * 2.0 * static_cast<SampleType>(idx) /
                                   static_cast<SampleType>(approximated_data.size());
                    for (auto i = 0; i < approximated_data.size(); i++) {
                        accum += fourier_series[i] * std::complex<SampleType>{std::cos(w * i), std::sin(w * i)};
                    }
                }
            }

            return accum;
        }

        template<typename SampleType, typename IdxType>
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

        template<typename SampleType, typename IdxType>
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

        void show(const PlottingKind kind) {
            using SampleType = double;
            is_actual = false;
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < signal_size; i++) {
                    plotting_data.push_back(compute<double, size_t>(i));
                }
                matplot::plot(plotting_data);
                matplot::show();
            }
        }

        void show(const PlottingKind kind, const std::string& filename, const std::string& format) {
            is_actual = false;
            using SampleType = double;
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < signal_size; i++) {
                    plotting_data.push_back(compute<double, size_t>(i));
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename, format);
            }
        }

        void show(const PlottingKind kind, const std::string& filename) {
            is_actual = false;
            using SampleType = double;
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < signal_size; i++) {
                    plotting_data.push_back(compute<double, size_t>(i));
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename);
            }
        }
    };

    template<typename T>
    struct ModifiedAkimaBasedWithNoTrain {
        //using boost::math::interpolators::makima;
        std::optional<boost::math::interpolators::makima<std::vector<double>>> spline = {};
        double min_bound = 0.0;
        double max_bound = 0.0;

        void loadData(const T & x, const T & y){
            std::vector<double> x_(x.size());
            std::vector<double> y_(x.size());
            for (auto i = 0; i < x.size(); i++){
                auto const & el = x[i];
                x_[i] = el;
            }
            for (auto i = 0; i < x.size(); i++){
                auto const & el = y[i];
                y_[i] = el;
            }
            double min = 9999999999999999999999.0; //todo
            double max = -9999999999999999999999.0; //todo
            for (auto const & el : x_){
                if (el < min){
                    min = el;
                }
                if (el > max){
                    max = el;
                }
            }
            min_bound = min;
            max_bound = max;
            spline = boost::math::interpolators::makima
                <std::vector<double>>(std::move(x_), std::move(y_));
        }

        void loadData(const T & y){
            //IC(y.size());
            std::vector<double> x_(y.size());
            std::vector<double> y_(y.size());
            for (int i = 0; i < y.size(); i++){
                x_[i] = (double)i;
            }
            for (int i = 0; i < y.size(); i++){
                auto const & el = y[i];
                y_[i] = el;
            }
            //IC(y_.size());
            min_bound = 0.0;
            max_bound = y.size() - 1;

            spline = boost::math::interpolators::makima
                <std::vector<double>>(std::move(x_), std::move(y_));
        }

        template<typename IdxT>
        double compute(IdxT idx){
            if (idx >= min_bound && idx <= max_bound){
                return (*spline)(idx);
            }
            else if (idx < min_bound){
                int64_t _idx = static_cast<int64_t>(idx);
                double idx_ = idx - _idx;
                auto size = max_bound - min_bound;
                int64_t _size = static_cast<int64_t>(size); //todo size_
                double new_idx = size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
                new_idx -= idx_;
                return (*spline)(new_idx);
            }
            else{
                int64_t _idx = static_cast<int64_t>(idx);
                double idx_ = idx - _idx;
                auto size = max_bound - min_bound;
                auto _size = static_cast<int64_t>(size); //todo size_
                double new_idx = _idx % _size + idx_;
                return (*spline)(new_idx);
            }
        }

        template<typename IdxT>
        double computeDerive(IdxT idx){
            if (idx >= min_bound && idx <= max_bound){
                return spline->prime(idx);
            }
            else if (idx < min_bound){
                int64_t _idx = static_cast<int64_t>(idx);
                double idx_ = idx - _idx;
                auto size = max_bound - min_bound;
                int64_t _size = static_cast<int64_t>(size); //todo size_
                double new_idx = size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
                new_idx -= idx_;
                return spline->prime(new_idx);
            }
            else{
                int64_t _idx = static_cast<int64_t>(idx);
                double idx_ = idx - _idx;
                auto size = max_bound - min_bound;
                auto _size = static_cast<int64_t>(size); //todo size_
                double new_idx = _idx % _size + idx_;
                return spline->prime(new_idx);
            }
        }
    };

    template<typename T>
    struct PiecewiseCubicHermitePolynomialBasedWithNoTrain {
        //using boost::math::interpolators::makima;
        std::optional<boost::math::interpolators::pchip<std::vector<double>>> spline = {};
        double min_bound = 0.0;
        double max_bound = 0.0;

        void loadData(const T & x, const T & y){
            std::vector<double> x_(x.size());
            std::vector<double> y_(x.size());
            for (auto i = 0; i < x.size(); i++){
                auto const & el = x[i];
                x_[i] = el;
            }
            for (auto i = 0; i < x.size(); i++){
                auto const & el = y[i];
                y_[i] = el;
            }
            double min = 9999999999999999999999.0; //todo
            double max = -9999999999999999999999.0; //todo
            for (auto const & el : x_){
                if (el < min){
                    min = el;
                }
                if (el > max){
                    max = el;
                }
            }
            min_bound = min;
            max_bound = max;
            spline = boost::math::interpolators::pchip
                <std::vector<double>>(std::move(x_), std::move(y_));
        }

        void loadData(const T & y){
            //IC(y.size());
            std::vector<double> x_(y.size());
            std::vector<double> y_(y.size());
            for (int i = 0; i < y.size(); i++){
                x_[i] = (double)i;
            }
            for (int i = 0; i < y.size(); i++){
                auto const & el = y[i];
                y_[i] = el;
            }
            //IC(y_.size());
            min_bound = 0.0;
            max_bound = y.size() - 1;

            spline = boost::math::interpolators::pchip
                <std::vector<double>>(std::move(x_), std::move(y_));
        }

        template<typename IdxT>
        double compute(IdxT idx){
            if (idx >= min_bound && idx <= max_bound){
                return (*spline)(idx);
            }
            else if (idx < min_bound){
                int64_t _idx = static_cast<int64_t>(idx);
                double idx_ = idx - _idx;
                auto size = max_bound - min_bound;
                int64_t _size = static_cast<int64_t>(size); //todo size_
                double new_idx = size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
                new_idx -= idx_;
                return (*spline)(new_idx);
            }
            else{
                int64_t _idx = static_cast<int64_t>(idx);
                double idx_ = idx - _idx;
                auto size = max_bound - min_bound;
                auto _size = static_cast<int64_t>(size); //todo size_
                double new_idx = _idx % _size + idx_;
                return (*spline)(new_idx);
            }
        }

        template<typename IdxT>
        double computeDerive(IdxT idx){
            if (idx >= min_bound && idx <= max_bound){
                return spline->prime(idx);
            }
            else if (idx < min_bound){
                int64_t _idx = static_cast<int64_t>(idx);
                double idx_ = idx - _idx;
                auto size = max_bound - min_bound;
                int64_t _size = static_cast<int64_t>(size); //todo size_
                double new_idx = size + (_idx - static_cast<int64_t>(min_bound)) % _size + min_bound;
                new_idx -= idx_;
                return spline->prime(new_idx);
            }
            else{
                int64_t _idx = static_cast<int64_t>(idx);
                double idx_ = idx - _idx;
                auto size = max_bound - min_bound;
                auto _size = static_cast<int64_t>(size); //todo size_
                double new_idx = _idx % _size + idx_;
                return spline->prime(new_idx);
            }
        }
    };
}
