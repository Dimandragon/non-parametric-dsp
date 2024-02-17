

module; 
#include <icecream.hpp>

import approximators;
import npdsp_concepts;
import npdsp_config;
import <utility>;
import <vector>;
import <cmath>;
import <numbers>;
import <complex>;
import utility_math;

export module phase_computers;

namespace NP_DSP{
    namespace ONE_D{
        namespace PHASE_COMPUTERS{
            export
            enum class InstFreqDerivativeBasedKind {Momental, TimeAverage, DeriveAverage, DeriveDouble};

            export 
            template <Signal DataT, Signal OutT>
            struct ExtremumsBasedNonOpt {
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = GENERAL::Nil;

                constexpr static bool is_phase_computer = true;

                void compute(const DataType & data, OutType & out, AdditionalDataType & nil){
                    std::vector<int> extremums;
                    extremums.push_back(0);
                    for (int i = 1; i < data.size()-1; i++){
                        if ((data[i] >= data[i-1] &&
                            data[i] > data[i+1]) ||
                            (data[i] > data[i-1] &&
                            data[i] >= data[i+1]) ||
                            (data[i] <= data[i-1] &&
                             data[i] < data[i+1]) ||
                            (data[i] < data[i-1] &&
                             data[i] <= data[i+1]))
                        {
                            extremums.push_back(i);
                        }
                    }
                    extremums.push_back(static_cast<int>(data.size()-1));

                    if(extremums.size() > 2){
                        extremums[0] = extremums[1] * 2 - extremums[2];
                        auto last = extremums.size() - 1;
                        extremums[last] = extremums[last - 1] * 2 - extremums[last - 2];
                        if (extremums[0] > 0){
                            extremums[0] = - extremums[1];
                        }
                        if (extremums[last] < extremums.size() - 1){
                            extremums[last] = data.size() + data.size() - extremums[last-1];
                        }
                    }
                    else{
                        auto last = extremums.size() - 1;
                        extremums[0] = - extremums[1];
                        extremums[last] = data.size() + data.size() - extremums[last-1];
                    }

                    std::vector<std::pair<int, double>> support;
                    for (auto i = 0; i < extremums.size(); i++){
                        support.push_back({extremums[i], i * std::numbers::pi});
                    }
                    
                    if constexpr (CONFIG::debug){
                        IC(support);
                    }

                    support[0] = {0, UTILITY_MATH::linearInterpolate<int, double>
                        (support[0], support[1], 0)};
                    auto pad = support[0].second;
                    for (auto i = 0; i < support.size(); i++){
                        support[i].second = support[i].second - pad;
                    }

                    support[support.size() - 1] = {data.size() - 1, UTILITY_MATH::linearInterpolate<int, double>
                        (support[support.size() - 2], support[support.size() - 1], data.size()-1)};

                    if constexpr (CONFIG::debug){
                        IC(support);
                    }

                    for(int i = 0; i < support.size() - 1; i++){
                        for (int j = support[i].first; j < support[i+1].first; j++){
                            out[j] = UTILITY_MATH::linearInterpolate<int, double>(support[i], support[i+1], j);
                        }
                    }
                    auto i = support.size() - 2;
                    auto j = support[i+1].first;
                    out[j] = UTILITY_MATH::linearInterpolate<int, double>(support[i], support[i+1], j);
                }
            };

            export
            template <Signal DataT, Signal OutT, Signal AdditionalDataT,
                Integrator IntegratorT, Derivator DerivatorT>
            struct ExtremumsBasedUsingFS{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = AdditionalDataT;
                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;
                int tile_size = 128;
                double approx_order_coeff = 1.0;

                constexpr static bool is_phase_computer = true;

                ExtremumsBasedUsingFS(IntegratorType integrator_in, DerivatorType derivator_in){
                    integrator = integrator_in;
                    derivator = derivator_in;
                }

                void compute(const DataType & data, OutType & out, AdditionalDataT & computer_buffer)
                {
                    //compute extremums
                    //compute "support" vector of pairs size_t and double, 
                    //where first is extremum position and second is value = idx * std::numbers::pi
                    //compute f(x) = integral(|arctg'(data')|)
                    //compute a(x) -> min for each pair in support 
                    // reduce f(position) * a(position) - support)^2 
                    // + for each x reduce
                    //if f(x) * a(x) > f(x-1) * a(x-1) then (f(x) * a(x) - f(x-1) * a(x-1))^2 else 0

                    std::vector<int> extremums;
                    extremums.push_back(0);
                    for (int i = 1; i < data.size()-1; i++){
                        if ((data[i] >= data[i-1] &&
                            data[i] > data[i+1]) ||
                            (data[i] > data[i-1] &&
                            data[i] >= data[i+1]) ||
                            (data[i] <= data[i-1] &&
                             data[i] < data[i+1]) ||
                            (data[i] < data[i-1] &&
                             data[i] <= data[i+1]))
                        {
                            extremums.push_back(i);
                        }
                    }
                    extremums.push_back(static_cast<int>(data.size()-1));
                    if constexpr (CONFIG::debug){
                        IC(extremums);
                    }
                    if(extremums.size() > 2){
                        extremums[0] = extremums[1] * 2 - extremums[2];
                        auto last = extremums.size() - 1;
                        extremums[last] = extremums[last - 1] * 2 - extremums[last - 2];
                        if (extremums[0] > 0){
                            extremums[0] = - extremums[1];
                        }
                        if (extremums[last] < extremums.size() - 1){
                            extremums[last] = data.size() + data.size() - extremums[last-1];
                        }
                    }
                    else{
                        auto last = extremums.size() - 1;
                        extremums[0] = - extremums[1];
                        extremums[last] = data.size() + data.size() - extremums[last-1];
                    }
                    if constexpr (CONFIG::debug){
                        IC(extremums);
                    }

                    std::vector<std::pair<size_t, float>> support;
                    for (auto i = 0; i < extremums.size(); i++){
                        support.push_back({extremums[i], i * std::numbers::pi});
                    }
                    if constexpr (CONFIG::debug){
                        IC(support);
                    }
                    
                    GENERAL::Nil nil;

                    derivator.compute(data, computer_buffer, nil);
                    for (int i = 0; i < data.size(); i++){
                        computer_buffer[i] = std::atan(computer_buffer[i]);
                    }
                    derivator.compute(computer_buffer, out, nil);
                    for (int i = 0; i < data.size(); i++){
                        out[i] = std::abs(computer_buffer[i]);
                    }
                    integrator.compute(out, computer_buffer, nil);
                    computer_buffer.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    std::vector<double> errors_cache;

                    auto loss = [&](auto & approximator){
                        double error = 0.0;
                        for (int i = 0; i < support.size(); i++){
                            auto position = support[i].first;
                            auto val = support[i].second;
                            error += 
                                (val - approximator.compute(position) * computer_buffer[position])
                                * (val - approximator.compute(position) * computer_buffer[position]) 
                                * computer_buffer.size() / support.size();
                        }
                        for (int i = 1; i < computer_buffer.size(); i++){
                            for (int j = 1; j <= i; j++){
                                if (computer_buffer[i]*approximator.compute(i) <
                                        computer_buffer[i-j]*approximator.compute(i-j)){
                                    error += (computer_buffer[i]*approximator.compute(i) - computer_buffer[i-j]*approximator.compute(i-j)) *
                                        (computer_buffer[i]*approximator.compute(i) - computer_buffer[i-j]*approximator.compute(i-j));
                                }
                            }
                        }
                        return error;
                    };

                    auto bySampleError = [&](auto & approximator, auto idx){
                        double error = 0.0;
                        for (int i = 0; i < support.size(); i++){
                            auto position = support[i].first;
                            auto val = support[i].second;
                            error += 
                                (val - approximator.compute(position) * computer_buffer[position])
                                * (val - approximator.compute(position) * computer_buffer[position]) 
                                * computer_buffer.size() / support.size();
                        }
                        
                        for (int i = 1; i < computer_buffer.size(); i++){
                            for (int j = 1; j <= i; j++){
                                if (computer_buffer[i]*approximator.compute(i) <
                                        computer_buffer[i-j]*approximator.compute(i-j)){
                                    error += (computer_buffer[i]*approximator.compute(i) - computer_buffer[i-j]*approximator.compute(i-j)) / i *
                                        (computer_buffer[i]*approximator.compute(i) - computer_buffer[i-j]*approximator.compute(i-j)) / computer_buffer[i];
                                }
                            }
                        }
                        //error += errors_cache[error_idx];
                        return error / approximator.tile_size;
                    };

                    auto stopPoint = [](auto losses_different, auto & approximator) {
                        if (losses_different > 0.00001){
                            return false;
                        }
                        else{
                            return true;
                        }
                    };

                    auto approximator = APPROX::FourierSeriesBased<DataType, decltype(loss), decltype(stopPoint), 
                        APPROX::FSApproxKind::Positive, decltype(bySampleError)>
                            (loss, out, stopPoint);
                    
                    approximator.bySampleLoss = &bySampleError;       
                    approximator.is_actual = false;
                    approximator.setApproxOrderRatio(approx_order_coeff);
                    approximator.tile_size = tile_size;
                    approximator.train();

                    approximator.is_actual = false;
                    for (auto i = 0; i < out.size(); i++){
                        out[i] = computer_buffer[i] * approximator.compute(i);
                    }
                    //todo test it and write example
                }
            };


            export
            template <Signal DataT, Signal OutT, Signal AdditionalDataT,
                Integrator IntegratorT, Derivator DerivatorT>
            struct ArctgScaledToExtremums{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = AdditionalDataT;
                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;
                int tile_size = 128;
                double approx_order_coeff = 1.0;

                constexpr static bool is_phase_computer = true;

                ArctgScaledToExtremums(IntegratorType integrator_in, DerivatorType derivator_in){
                    integrator = integrator_in;
                    derivator = derivator_in;
                }

                void compute(const DataType & data, OutType & out, AdditionalDataT & computer_buffer)
                {
                    //compute extremums
                    //compute "support" vector of pairs size_t and double, 
                    //where first is extremum position and second is value = idx * std::numbers::pi
                    //compute f(x) = integral(|arctg'(data')|)

                    std::vector<int> extremums;
                    extremums.push_back(0);
                    for (int i = 1; i < data.size()-1; i++){
                        if ((data[i] >= data[i-1] &&
                            data[i] > data[i+1]) ||
                            (data[i] > data[i-1] &&
                            data[i] >= data[i+1]) ||
                            (data[i] <= data[i-1] &&
                             data[i] < data[i+1]) ||
                            (data[i] < data[i-1] &&
                             data[i] <= data[i+1]))
                        {
                            extremums.push_back(i);
                        }
                    }
                    extremums.push_back(static_cast<int>(data.size()-1));

                    if(extremums.size() > 2){
                        extremums[0] = extremums[1] * 2 - extremums[2];
                        auto last = extremums.size() - 1;
                        extremums[last] = extremums[last - 1] * 2 - extremums[last - 2];
                        if (extremums[0] > 0){
                            extremums[0] = - extremums[1];
                        }
                        if (extremums[last] < extremums.size() - 1){
                            extremums[last] = data.size() + data.size() - extremums[last-1];
                        }
                    }
                    else{
                        auto last = extremums.size() - 1;
                        extremums[0] = - extremums[1];
                        extremums[last] = data.size() + data.size() - extremums[last-1];
                    }

                    std::vector<std::pair<int, double>> support;
                    for (auto i = 0; i < extremums.size(); i++){
                        support.push_back({extremums[i], i * std::numbers::pi});
                    }
                    
                    if constexpr (CONFIG::debug){
                        IC(support);
                    }

                    GENERAL::Nil nil;

                    derivator.compute(data, computer_buffer, nil);
                    for (int i = 0; i < data.size(); i++){
                        computer_buffer[i] = std::atan(computer_buffer[i]);
                    }
                    derivator.compute(computer_buffer, out, nil);
                    for (int i = 0; i < data.size(); i++){
                        computer_buffer[i] = std::abs(out[i]);
                    }
                    integrator.compute(computer_buffer, out, nil);
                    if constexpr (CONFIG::debug){
                        out.show(NP_DSP::ONE_D::PlottingKind::Simple);
                    }
                    
                    //norm 
                    support[0] = {0, UTILITY_MATH::linearInterpolate<int, double>
                        (support[0], support[1], 0)};
                    auto pad = support[0].second;
                    for (auto i = 0; i < support.size(); i++){
                        support[i].second = support[i].second - pad;
                    }

                    support[support.size() - 1] = {data.size() - 1, UTILITY_MATH::linearInterpolate<int, double>
                        (support[support.size() - 2], support[support.size() - 1], data.size()-1)};

                    if constexpr (CONFIG::debug){
                        IC(support);
                    }

                    std::vector<std::pair<int, double>> div_support;
                    
                    for (auto i = 0; i < support.size(); i = i + 1){
                        //IC(i, support.size());
                        if (i == 0){
                            div_support.push_back({0, 1});
                        }
                        else{
                            div_support.push_back({support[i].first, 
                                support[i].second / out[support[i].first]});
                        }
                    }
                    
                    div_support[0].second = div_support[1].second;

                    if constexpr (CONFIG::debug){
                        IC(div_support);
                    }

                    for(int i = 0; i < div_support.size() - 1; i++){
                        for (int j = div_support[i].first; j < div_support[i+1].first; j++){
                            out[j] = out[j] * 
                                UTILITY_MATH::linearInterpolate<int, double>(div_support[i], div_support[i+1], j);
                        }
                    }
                    auto i = div_support.size() - 2;
                    auto j = div_support[i+1].first;
                    out[j] = out[j] * 
                                UTILITY_MATH::linearInterpolate<int, double>(div_support[i], div_support[i+1], j);
                    /*
                    //todo opt by smooth using div_support
                    std::vector<double> der_div_support_points;
                    for (int i = 0; i < div_support.size() - 1; i++){
                        der_div_support_points.push_back((div_supports[i + 1].second - div_support[i].second) / 
                            (div_supports[i + 1].first - div_support[i].first));
                    }
                    */
                }
            };

            export
            template <Signal DataT, Signal OutT, Signal AdditionalDataT,
                Integrator IntegratorT, Derivator DerivatorT, InstFreqDerivativeBasedKind kind>
            struct ArctgScaledToExtremumsSquare{
                using DataType = DataT;
                using OutType = OutT;
                using AdditionalDataType = AdditionalDataT;
                using IntegratorType = IntegratorT;
                using DerivatorType = DerivatorT;


                IntegratorType integrator;
                DerivatorType derivator;

                constexpr static bool is_phase_computer = true;

                ArctgScaledToExtremumsSquare(IntegratorType integrator_in, DerivatorType derivator_in){
                    integrator = integrator_in;
                    derivator = derivator_in;
                }

                void compute(const DataType & data, OutType & out, AdditionalDataT & computer_buffer)
                {
                    //compute extremums
                    //compute "support" vector of pairs size_t and double,
                    //where first is extremum position and second is value = idx * std::numbers::pi
                    //compute f(x) = integral(|arctg'(data')|)

                    std::vector<int> extremums;
                    extremums.push_back(0);
                    for (int i = 1; i < data.size()-1; i++){
                        if ((data[i] >= data[i-1] &&
                            data[i] > data[i+1]) ||
                            (data[i] > data[i-1] &&
                            data[i] >= data[i+1]) ||
                            (data[i] <= data[i-1] &&
                             data[i] < data[i+1]) ||
                            (data[i] < data[i-1] &&
                             data[i] <= data[i+1]))
                        {
                            extremums.push_back(i);
                        }
                    }
                    extremums.push_back(static_cast<int>(data.size()-1));

                    if(extremums.size() > 2){
                        extremums[0] = extremums[1] * 2 - extremums[2];
                        auto last = extremums.size() - 1;
                        extremums[last] = extremums[last - 1] * 2 - extremums[last - 2];
                        if (extremums[0] > 0){
                            extremums[0] = - extremums[1];
                        }
                        if (extremums[last] < extremums.size() - 1){
                            extremums[last] = data.size() + data.size() - extremums[last-1];
                        }
                    }
                    else{
                        auto last = extremums.size() - 1;
                        extremums[0] = - extremums[1];
                        extremums[last] = data.size() + data.size() - extremums[last-1];
                    }

                    std::vector<std::pair<int, double>> support;
                    for (auto i = 0; i < extremums.size(); i++){
                        support.push_back({extremums[i], i * std::numbers::pi});
                    }

                    if constexpr (CONFIG::debug){
                        IC(support);
                    }

                    GENERAL::Nil nil;

                    derivator.compute(data, computer_buffer, nil);
                    for (int i = 0; i < data.size(); i++){
                        computer_buffer[i] = std::atan(computer_buffer[i]);
                    }
                    derivator.compute(computer_buffer, out, nil);
                    for (int i = 0; i < data.size(); i++){
                        out[i] = std::abs(out[i]);
                    }
                    integrator.compute(out, computer_buffer, nil);
                    if constexpr (CONFIG::debug){
                        computer_buffer.show(PlottingKind::Simple);
                    }

                    if constexpr (kind == InstFreqDerivativeBasedKind::TimeAverage) {
                        for (auto i = 0; i < data.size(); i++){
                            auto approx_answer = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer = approx_answer;
                            auto counter = 0;
                            while (approx_answer < 2.0 * std::numbers::pi){
                                counter++;
                                old_approx_answer = approx_answer;
                                approx_answer = computer_buffer.interpolate(i+counter) - computer_buffer.interpolate(i - counter);
                            }
                            auto left_loss = std::numbers::pi * 2.0 - old_approx_answer;
                            auto right_loss = approx_answer - std::numbers::pi * 2.0;
                            auto sum_loss = left_loss + right_loss;
                            auto period = (static_cast<OutType::SampleType>(counter) - right_loss/sum_loss) * 2;
                            out[i] = static_cast<OutType::SampleType>(1.0/period);
                        }

                        for (int i = 0; i < data.size(); i++) {
                            computer_buffer[i] = out[i] / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(computer_buffer, out, nil);
                    }
                    else if constexpr (kind == InstFreqDerivativeBasedKind::DeriveAverage) {
                        for (auto i = 0; i < data.size(); i++){
                            auto approx_answer_right = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer_right = approx_answer_right;
                            auto counter = 0;
                            while (approx_answer_right < std::numbers::pi){
                                counter++;
                                old_approx_answer_right = approx_answer_right;
                                approx_answer_right = computer_buffer.interpolate(i + counter) - computer_buffer.interpolate(i);
                            }
                            auto left_loss_right_edge = std::numbers::pi - old_approx_answer_right;
                            auto right_loss_right_edge = approx_answer_right - std::numbers::pi;
                            auto sum_loss_right_edge = left_loss_right_edge + right_loss_right_edge;
                            auto right_edge = static_cast<OutType::SampleType>(counter) - right_loss_right_edge/sum_loss_right_edge;

                            auto approx_answer_left = static_cast<OutType::SampleType>(0.0);
                            auto old_approx_answer_left = approx_answer_left;
                            counter = 0;
                            while (approx_answer_left < std::numbers::pi){
                                counter++;
                                old_approx_answer_left = approx_answer_left;
                                approx_answer_left = computer_buffer.interpolate(i)-computer_buffer.interpolate(i-counter);
                            }
                            auto left_loss_left_edge = std::numbers::pi - old_approx_answer_left;
                            auto right_loss_left_edge = approx_answer_left - std::numbers::pi;
                            auto sum_loss_left_edge = left_loss_left_edge + right_loss_left_edge;
                            auto left_edge = static_cast<OutType::SampleType>(counter) - right_loss_left_edge/sum_loss_left_edge;

                            auto period = right_edge + left_edge;
                            out[i] = static_cast<OutType::SampleType>(1.0/period);
                        }

                        for (int i = 0; i < data.size(); i++) {
                            computer_buffer[i] = out[i] / (std::numbers::pi * 2.0);
                        }
                        integrator.compute(computer_buffer, out, nil);
                    }
                    else if constexpr (kind == InstFreqDerivativeBasedKind::Momental) {
                        for (int i = 0; i < data.size(); i++) {
                            out[i] = computer_buffer[i];
                        }
                    }

                    //norm
                    support[0] = {0, UTILITY_MATH::linearInterpolate<int, double>
                        (support[0], support[1], 0)};
                    auto pad = support[0].second;
                    for (auto i = 0; i < support.size(); i++){
                        support[i].second = support[i].second - pad;
                    }

                    support[support.size() - 1] = {data.size() - 1, UTILITY_MATH::linearInterpolate<int, double>
                        (support[support.size() - 2], support[support.size() - 1], data.size()-1)};

                    if constexpr (CONFIG::debug){
                        IC(support);
                    }

                    std::vector<std::pair<int, double>> div_support;

                    for (auto i = 0; i < support.size(); i = i + 1){
                        //IC(i, support.size());
                        if (i == 0){
                            div_support.push_back({0, 1});
                        }
                        else{
                            div_support.push_back({support[i].first,
                                support[i].second / out[support[i].first]});
                        }
                    }

                    div_support[0].second = div_support[1].second;

                    if constexpr (CONFIG::debug){
                        IC(div_support);
                    }

                    for(int i = 0; i < div_support.size() - 1; i++){
                        for (int j = div_support[i].first; j < div_support[i+1].first; j++){
                            out[j] = out[j] *
                                UTILITY_MATH::linearInterpolate<int, double>(div_support[i], div_support[i+1], j);
                        }
                    }
                    auto i = div_support.size() - 2;
                    auto j = div_support[i+1].first;
                    out[j] = out[j] *
                                UTILITY_MATH::linearInterpolate<int, double>(div_support[i], div_support[i+1], j);
                    /*
                    //todo opt by smooth using div_support
                    std::vector<double> der_div_support_points;
                    for (int i = 0; i < div_support.size() - 1; i++){
                        der_div_support_points.push_back((div_supports[i + 1].second - div_support[i].second) /
                            (div_supports[i + 1].first - div_support[i].first));
                    }
                    */
                }
            };
        }
    }
}