module;

import npdsp_concepts;
import <vector>;
import <complex>;
import utility_math;
import <numbers>;

export module approximators;

namespace NP_DSP{
    namespace ONE_D{
        namespace APPROX{
            export
            template<Signal SignalT, typename LossFunc, typename StopPointFunc>
            struct FourerSeriesBased{
                using SignalType = SignalT;
                using DataType = SignalT::DataType;
                using Loss = LossFunc;
                using StopPoint = StopPointFunc;
                using IdxType = SignalT::IdxType;
                static constexpr bool is_signal_approximator = true;

                Loss loss;
                StopPoint stopPont;
                SignalType signal;
                bool is_actual = false;


                std::vector<std::complex<DataType>> fourier_series;
                std::vector<std::complex<DataType>> approximated_data;

                DataType max_value = 1000000;
                DataType max_error;

                FourerSeriesBased(Loss lossFn, StopPoint stopPointFn, SignalType signal_in, DataType max_error_in){
                    loss = lossFn;
                    stopPont = stopPointFn;
                    signal = signal_in;
                    fourier_series = new std::vector<std::complex<DataType>>(signal_in.size());
                    approximated_data = new std::vector<std::complex<DataType>>(signal_in.size());
                    max_error = max_error_in;
                }

                ~FourerSeriesBased(){
                    delete fourier_series;
                    delete approximated_data;
                }

                void computeData(){
                    is_actual = true;
                    UTILITY_MATH::ifftc2c(fourier_series, approximated_data);
                }

                inline std::complex<DataType> computeSample(IdxType idx){
                    std::complex<DataType> accum = {static_cast<DataType>(0), static_cast<DataType>(0)};
                    DataType w = std::numbers::pi * 2.0 * static_cast<DataType>(idx) / static_cast<DataType>(fourier_series.size());
                    for (auto i = 0; i < fourier_series.size(); i++){
                        accum+= fourier_series[i]*std::complex<DataType>{std::cos(w*i), std::sin(w*i)};
                    }
                    return accum;
                }

                DataType compute(IdxType idx){
                    if(idx % 0 == 0){
                        if(is_actual){
                            return approximated_data[idx].real;
                        }
                        else{
                            computeData();
                            return approximated_data[idx].real;
                        }
                    }
                    else{
                        return computeSample(idx).real;
                    }
                }

                void train() {
                    for (auto i = 0; i < fourier_series.size(); i++) {
                        //auto triginimetric_sample = ONE_D::UTILITY_MATH::convertFSampleC2T(fourier_series[i]);
                        std::pair<DataType, DataType> triginimetric_sample;
                        //optimize
                        // optimize theta
                        auto theta_min = static_cast<DataType>(-std::numbers::pi/2);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C({static_cast<DataType>(1.0), theta_min});
                        auto left_loss_pos = loss(*this);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C({static_cast<DataType>(-1.0), theta_min});
                        auto left_loss_neg = loss(*this);
                        auto theta_max = static_cast<DataType>(std::numbers::pi/2);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C({static_cast<DataType>(1.0), theta_max});
                        auto right_loss_pos = loss(*this);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C({static_cast<DataType>(-1.0), theta_max});
                        auto right_loss_neg = loss(*this);
                        auto theta_central = static_cast<DataType>(0.0);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C({static_cast<DataType>(1.0), theta_central});
                        auto central_loss_pos = loss(*this);
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C({static_cast<DataType>(-1.0), theta_central});
                        auto central_loss_neg = loss(*this);

                        DataType right_loss;
                        DataType left_loss;
                        DataType central_loss;

                        if (central_loss_neg < central_loss_pos || right_loss_neg < right_loss_pos || left_loss_neg < left_loss_pos) {
                            triginimetric_sample.first = - 1.0;
                            right_loss = right_loss_neg;
                            left_loss = left_loss_neg;
                            central_loss = central_loss_neg;
                        }
                        else {
                            triginimetric_sample.first = 1.0;
                            right_loss = right_loss_pos;
                            left_loss = left_loss_pos;
                            central_loss = central_loss_pos;
                        }


                        while (std::abs(right_loss - left_loss) > max_error) {
                            if (left_loss > central_loss && central_loss > right_loss) {
                                left_loss = central_loss;
                                theta_min = theta_central;
                                theta_central = (theta_min + theta_max) / 2;
                                triginimetric_sample.second = theta_central;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                central_loss = loss(*this);
                            }
                            else if (left_loss < central_loss && central_loss < right_loss) {
                                right_loss = central_loss;
                                theta_max = theta_central;
                                theta_central = (theta_min + theta_max) / 2;
                                triginimetric_sample.second = theta_central;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                central_loss = loss(*this);
                            }
                            else if (left_loss - central_loss > right_loss - central_loss) {
                                //left branch is higher
                                auto theta_left_avg = (theta_min + theta_central) / 2;
                                triginimetric_sample.second = theta_left_avg;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                auto new_left_loss = loss(*this);
                                if(new_left_loss < left_loss && new_left_loss > central_loss) {
                                    left_loss = new_left_loss;
                                    theta_min = theta_left_avg;
                                    theta_central = (theta_min + theta_max) / 2;
                                    triginimetric_sample.second = theta_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    central_loss = loss(*this);
                                }
                                else {
                                    right_loss = central_loss;
                                    theta_max = theta_central;
                                    theta_central = theta_left_avg;
                                    central_loss = new_left_loss;
                                }
                            }
                            else {
                                //right branch is higher
                                auto theta_right_avg = (theta_max + theta_central) / 2;
                                triginimetric_sample.second = theta_right_avg;
                                fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                auto new_right_loss = loss(*this);

                                if (new_right_loss < right_loss && new_right_loss > central_loss) {
                                    right_loss = new_right_loss;
                                    theta_max = theta_right_avg;
                                    theta_central = (theta_max + theta_min) / 2;
                                    triginimetric_sample.second = theta_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    central_loss = loss(*this);
                                }
                                else {
                                    left_loss = central_loss;
                                    theta_min = theta_central;
                                    central_loss = new_right_loss;
                                    theta_central = theta_right_avg;
                                }
                            }
                        }
                        triginimetric_sample.second = theta_central;

                        //optimize ampl
                        if (triginimetric_sample.first == 1.0) {
                            DataType ampl_left = 0.0;
                            DataType ampl_right = max_value;
                            DataType ampl_central = (ampl_left + ampl_right) / 2;
                            triginimetric_sample.first = ampl_left;
                            fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                            auto loss_left = loss(*this);
                            triginimetric_sample.first = ampl_right;
                            fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                            auto loss_right = loss(*this);
                            triginimetric_sample.first = ampl_central;
                            fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                            auto loss_central = loss(*this);

                            while(std::abs(loss_right - loss_left) > max_error) {
                                if(loss_left > loss_central && loss_central > loss_right) {
                                    loss_left = loss_central;
                                    ampl_left = ampl_central;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    triginimetric_sample.first = ampl_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    loss_central = loss(*this);
                                }
                                else if (loss_left < loss_central && loss_central < loss_right) {
                                    loss_right = loss_central;
                                    ampl_right = ampl_central;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    triginimetric_sample.first = ampl_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    loss_central = loss(*this);
                                }
                                else if (loss_left - loss_central > loss_right - loss_central) {
                                    //left branch is higher
                                    auto ampl_left_avg = (ampl_left + ampl_central) / 2;
                                    triginimetric_sample.first = ampl_left_avg;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    auto new_loss_left = loss(*this);
                                    if (new_loss_left < loss_left && new_loss_left > central_loss) {
                                        loss_left = new_loss_left;
                                        ampl_left = ampl_left_avg;
                                        ampl_central = (ampl_left + ampl_right) / 2;
                                        triginimetric_sample.first = ampl_central;
                                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                        loss_central = loss(*this);
                                    }
                                    else {
                                        loss_right = loss_central;
                                        ampl_right = ampl_central;
                                        ampl_central = ampl_left_avg;
                                        loss_central = new_loss_left;
                                    }
                                }
                                else {
                                    auto ampl_right_avg = (ampl_right + ampl_central) / 2;
                                    triginimetric_sample.first = ampl_right_avg;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    auto new_loss_right = loss(*this);

                                    if (new_loss_right < loss_right && new_loss_right > loss_central) {
                                        loss_right = new_loss_right;
                                        ampl_right = ampl_right_avg;
                                        ampl_central = (ampl_left + ampl_right) / 2;
                                        triginimetric_sample.first = ampl_central;
                                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                        loss_central = loss(*this);
                                    }
                                    else {
                                        loss_left = loss_central;
                                        ampl_left = ampl_central;
                                        loss_central = new_loss_right;
                                        ampl_central = ampl_right_avg;
                                    }
                                }
                            }
                            triginimetric_sample.first = ampl_central;
                        }
                        else {
                            DataType ampl_left = - max_value;
                            DataType ampl_right = 0;
                            DataType ampl_central = (ampl_left + ampl_right) / 2;
                            triginimetric_sample.first = ampl_left;
                            fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                            auto loss_left = loss(*this);
                            triginimetric_sample.first = ampl_right;
                            fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                            auto loss_right = loss(*this);
                            triginimetric_sample.first = ampl_central;
                            fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                            auto loss_central = loss(*this);

                            while(std::abs(loss_right - loss_left) > max_error) {
                                if(loss_left > loss_central && loss_central > loss_right) {
                                    loss_left = loss_central;
                                    ampl_left = ampl_central;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    triginimetric_sample.first = ampl_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    loss_central = loss(*this);
                                }
                                else if (loss_left < loss_central && loss_central < loss_right) {
                                    loss_right = loss_central;
                                    ampl_right = ampl_central;
                                    ampl_central = (ampl_left + ampl_right) / 2;
                                    triginimetric_sample.first = ampl_central;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    loss_central = loss(*this);
                                }
                                else if (loss_left - loss_central > loss_right - loss_central) {
                                    //left branch is higher
                                    auto ampl_left_avg = (ampl_left + ampl_central) / 2;
                                    triginimetric_sample.first = ampl_left_avg;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    auto new_loss_left = loss(*this);
                                    if (new_loss_left < loss_left && new_loss_left > central_loss) {
                                        loss_left = new_loss_left;
                                        ampl_left = ampl_left_avg;
                                        ampl_central = (ampl_left + ampl_right) / 2;
                                        triginimetric_sample.first = ampl_central;
                                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                        loss_central = loss(*this);
                                    }
                                    else {
                                        loss_right = loss_central;
                                        ampl_right = ampl_central;
                                        ampl_central = ampl_left_avg;
                                        loss_central = new_loss_left;
                                    }
                                }
                                else {
                                    auto ampl_right_avg = (ampl_right + ampl_central) / 2;
                                    triginimetric_sample.first = ampl_right_avg;
                                    fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                    auto new_loss_right = loss(*this);

                                    if (new_loss_right < loss_right && new_loss_right > loss_central) {
                                        loss_right = new_loss_right;
                                        ampl_right = ampl_right_avg;
                                        ampl_central = (ampl_left + ampl_right) / 2;
                                        triginimetric_sample.first = ampl_central;
                                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                                        loss_central = loss(*this);
                                    }
                                    else {
                                        loss_left = loss_central;
                                        ampl_left = ampl_central;
                                        loss_central = new_loss_right;
                                        ampl_central = ampl_right_avg;
                                    }
                                }
                            }
                            triginimetric_sample.first = ampl_central;
                        }
                        fourier_series[i] = UTILITY_MATH::convertFSampleT2C(triginimetric_sample);
                    }
                }
            };
        }
    }
}