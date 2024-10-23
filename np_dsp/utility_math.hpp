#pragma once

#include "pocketfft_hdronly.h"

#include <utility>
#include <complex>
#include <vector>

#include <npdsp_concepts.hpp>
#include <npdsp_config.hpp>

namespace NP_DSP::ONE_D::UTILITY_MATH {   
    struct SquarePolynome{
        double a;
        double b;
        double c;

        void solve(auto x1, auto x2, auto x3, auto y1, auto y2, auto y3){
            b = ((y1-y3)*(x1*x1 - x2*x2) - (y1 - y2)*(x1*x1 - x3*x3))/
            ((x1 - x2)*(x1 - x3)*(x2 - x3));
            a = (y1 - y2 - b * (x1 - x2)) / (x1*x1 - x2*x2);
            c = y3 - a*x3*x3 - b*x3;
        }

        double compute(auto x){
            return a*x*x + b*x + c;
        }

        double derive(auto x){
            return b + 2.0 * a;
        }
    };

    struct Linear{
        double k;
        double b;
    
        void solve(auto x1, auto x2, auto y1, auto y2){
            auto dx = x2 - x1;
            auto dy = y2 - y1;

            k = dy/dx;
            b = x1 - k*x1;
        }
    
        double compute(auto x){
            return k*x + b;
        }
        double derive(auto x){
            return k;
        }
    };

    template<typename T>
    T complexL2(std::complex<T> a, std::complex<T> b) {
        return (a.imag() - b.imag()) * (a.imag() - b.imag()) + (a.real() - b.real()) * (a.real() - b.real());
    }
    
    template<typename T>
    T pairL2(std::pair<T, T> a, std::pair<T, T> b) {
        return (a.first - b.first) * (a.first - b.first) + (a.second - b.second) * (a.second - b.second);
    }

    auto sec(auto z_r) {
        return 1 / std::cos(z_r);
    }

    
    template<typename T, typename data1T, typename data2T>
    T signalsL2Distance(const data1T& data1, const data2T& data2) {
        double error = 0.0;
        error = 0.0;
        for (int i = 0; i < data1.size(); i++) {
            error = error + std::sqrt((data1[i] - data2[i]) * (data1[i] - data2[i]));
        }
        return error;
    }
    
    template<typename T, typename data1T, typename data2T>
    T signalsL2NormedDistance(const data1T& data1, const data2T& data2) {
        double error = 0.0;
        for (int i = 0; i < data1.size(); i++) {
            error = error + std::sqrt((data1[i] - data2[i]) * (data1[i] - data2[i])) / data1.size();
        }
        return error;
    }

    
    template<typename T, typename data1T, typename data2T>
    T signalsL2DistanceDouble(const data1T& data1, const data2T& data2) {
        double error = 0.0;
        error = 0.0;
        for (int i = 0; i < data1.size(); i++) {
            error = error + std::sqrt((data1[i].first - data2[i].first) * (data1[i].first - data2[i].first) + 
                (data1[i].second - data2[i].second) * (data1[i].second - data2[i].second));
        }
        return error;
    }

    
    template<typename T, typename data1T, typename data2T>
    T signalsL2NormedDistanceDouble(const data1T& data1, const data2T& data2) {
        double error = 0.0;
        error = 0.0;
        for (int i = 0; i < data1.size(); i++) {
            error = error + std::sqrt((data1[i].first - data2[i].first) * (data1[i].first - data2[i].first) + 
                (data1[i].second - data2[i].second) * (data1[i].second - data2[i].second)) / data1.size();
        }
        return error;
    }

    
    template<typename xType, typename yType>
    yType linearInterpolate(std::pair<xType, yType> point1, std::pair<xType, yType> point2, xType x_in) {
        if (point2.first < point1.first){
            auto swap = point1;
            point1 = point2;
            point2 = swap;
        }
        auto dx = point2.first - point1.first;
        auto dy = point2.second - point1.second;
        if (dx == 0) {
            return (point1.second + point2.second) / 2;
        }
        return point1.second + dy * (x_in - point1.first) / dx;
    }

    
    template<typename xType, typename yType>
    xType backLinearInterpolate(std::pair<xType, yType> point1, std::pair<xType, yType> point2, yType y_in) {
        auto dx = point2.first - point1.first;
        auto dy = point2.second - point1.second;
        if (dy == 0) {
            return (point1.first + point2.first) / 2;
        }
        return point1.first + dx * (y_in - point1.second) / dy;
    }

    
    template<typename T>
    void fftc2c(std::vector<std::complex<T>> const& data_in, std::vector<std::complex<T>>& data_out) {
        auto len = data_in.size();
        pocketfft::shape_t shape{len};
        pocketfft::stride_t stridef(shape.size());
        size_t tmpf = sizeof(std::complex<T>);
        for (int i = shape.size() - 1; i >= 0; --i) {
            stridef[i] = tmpf;
            tmpf *= shape[i];
        }
        pocketfft::shape_t axes;
        for (size_t i = 0; i < shape.size(); ++i) {
            axes.push_back(i);
        }
        pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::FORWARD,
                       data_in.data(), data_out.data(), (T) 1.0 / (T) len);
    }

    
    template<typename T>
    void ifftc2c(std::vector<std::complex<T>> const& data_in, std::vector<std::complex<T>>& data_out) {
        auto len = data_in.size();
        pocketfft::shape_t shape{len};
        pocketfft::stride_t stridef(shape.size());
        size_t tmpf = sizeof(std::complex<T>);
        for (int i = shape.size() - 1; i >= 0; --i) {
            stridef[i] = tmpf;
            tmpf *= shape[i];
        }
        pocketfft::shape_t axes;
        for (size_t i = 0; i < shape.size(); ++i) {
            axes.push_back(i);
        }
        pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::BACKWARD,
                       data_in.data(), data_out.data(), static_cast<T>(1));
    }

    
    template<Signal DataT, Signal OutT, typename T>
    void fftc2c(const DataT& in, OutT& out) {
        //using T = typename OutT::SampleType;
        auto len = in.size();
        pocketfft::shape_t shape{len};
        pocketfft::stride_t stridef(shape.size());
        size_t tmpf = sizeof(std::complex<T>);
        for (int i = shape.size() - 1; i >= 0; --i) {
            stridef[i] = tmpf;
            tmpf *= shape[i];
        }
        pocketfft::shape_t axes;
        for (size_t i = 0; i < shape.size(); ++i) {
            axes.push_back(i);
        }
        std::vector<std::complex<T>> data_in;
        std::vector<std::complex<T>> data_out;
        for (int i = 0; i < len; i++) {
            data_in.push_back(in[i]);
            data_out.push_back({0., 0.});
        }
        pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::FORWARD,
                       data_in.data(), data_out.data(), 1.0 / len);

        for (int i = 0; i < len; i++) {
            out[i] = data_out[i];
        }
    }

    
    template<Signal DataT, Signal OutT, typename T>
    void ifftc2c(const DataT& in, OutT& out) {
        auto len = in.size();
        pocketfft::shape_t shape{len};
        pocketfft::stride_t stridef(shape.size());
        size_t tmpf = sizeof(std::complex<T>);
        for (int i = shape.size() - 1; i >= 0; --i) {
            stridef[i] = tmpf;
            tmpf *= shape[i];
        }
        pocketfft::shape_t axes;
        for (size_t i = 0; i < shape.size(); ++i) {
            axes.push_back(i);
        }
        std::vector<std::complex<T>> data_in;
        std::vector<std::complex<T>> data_out;
        for (int i = 0; i < len; i++) {
            data_in.push_back(in[i]);
            data_out.push_back({0., 0.});
        }
        pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::BACKWARD,
                       data_in.data(), data_out.data(), static_cast<T>(1));
        for (int i = 0; i < len; i++) {
            out[i] = data_out[i];
        }
    }

    template<typename DataT1, typename DataT2, typename OutT, typename T>
    void fastConvolution(const DataT1& in1, const DataT2& in2, OutT& out) {
        std::vector<std::complex<T>> data_in1;
        std::vector<std::complex<T>> data_in2;
        std::vector<std::complex<T>> sp1;
        std::vector<std::complex<T>> sp2;
        std::vector<std::complex<T>> data_out;

        size_t len = in1.size();
        if (in2.size() > len) {
            len = in2.size();
        }

        for (int i = 0; i < in1.size(); i++) {
            data_in1.push_back({in1[i], 0.0});
        }
        for (int i = 0; i < in2.size(); i++) {
            data_in2.push_back({in2[i], 0.0});
        }
        for (int i = in1.size(); i < len; i++) {
            data_in1.push_back({0.0, 0.0});
        }
        for (int i = in2.size(); i < len; i++) {
            data_in2.push_back({0.0, 0.0});
        }
        for (int i = 0; i < len; i++) {
            sp1.push_back({0.0, 0.0});
            sp2.push_back({0.0, 0.0});
            data_out.push_back({0.0, 0.0});
        }

        fftc2c(data_in1, sp1);
        fftc2c(data_in2, sp2);

        for (int i = 0; i < len; i++) {
            sp1[i] = sp1[i] * sp2[i] * static_cast<T>(len);
        }

        ifftc2c(sp1, data_out);

        for (int i = 0; i < len; i++) {
            out[i] = data_out[i].real();
        }
    }

    
    template<Signal DataT1, Signal DataT2, Signal OutT>
    void fastConvolution(const DataT1& in1, const DataT2& in2, OutT& out) {
        using T = typename OutT::SampleType;
        std::vector<std::complex<T>> data_in1;
        std::vector<std::complex<T>> data_in2;
        std::vector<std::complex<T>> sp1;
        std::vector<std::complex<T>> sp2;
        std::vector<std::complex<T>> data_out;

        size_t len = in1.size();
        if (in2.size() > len) {
            len = in2.size();
        }

        for (int i = 0; i < in1.size(); i++) {
            data_in1.push_back({in1[i], 0.0});
        }
        for (int i = 0; i < in2.size(); i++) {
            data_in2.push_back({in2[i], 0.0});
        }
        for (int i = in1.size(); i < len; i++) {
            data_in1.push_back({0.0, 0.0});
        }
        for (int i = in2.size(); i < len; i++) {
            data_in2.push_back({0.0, 0.0});
        }
        for (int i = 0; i < len; i++) {
            sp1.push_back({0.0, 0.0});
            sp2.push_back({0.0, 0.0});
            data_out.push_back({0.0, 0.0});
        }

        fftc2c(data_in1, sp1);
        fftc2c(data_in2, sp2);

        for (int i = 0; i < len; i++) {
            sp1[i] = sp1[i] * sp2[i] * static_cast<T>(len);
        }

        ifftc2c(sp1, data_out);

        for (int i = 0; i < len; i++) {
            out[i] = data_out[i].real();
        }
    }

    
    template<typename T>
    void fftc2c(std::vector<std::complex<T>> const& data_in, 
        std::vector<std::complex<T>>& data_out, size_t len, size_t pad) {
        pocketfft::shape_t shape{len};
        pocketfft::stride_t stridef(shape.size());
        size_t tmpf = sizeof(std::complex<T>);
        for (int i = shape.size() - 1; i >= 0; --i) {
            stridef[i] = tmpf;
            tmpf *= shape[i];
        }
        pocketfft::shape_t axes;
        for (size_t i = 0; i < shape.size(); ++i) {
            axes.push_back(i);
        }
        pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::FORWARD,
                       data_in.data() + pad, data_out.data() + pad, ((T) 1.0) / ((T) len));
    }

    
    template<typename T>
    void ifftc2c(std::vector<std::complex<T>> const& data_in, std::vector<std::complex<T>>& data_out, size_t len,
                 size_t pad) {
        pocketfft::shape_t shape{len};
        pocketfft::stride_t stridef(shape.size());
        size_t tmpf = sizeof(std::complex<T>);
        for (int i = shape.size() - 1; i >= 0; --i) {
            stridef[i] = tmpf;
            tmpf *= shape[i];
        }
        pocketfft::shape_t axes;
        for (size_t i = 0; i < shape.size(); ++i) {
            axes.push_back(i);
        }
        pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::BACKWARD,
                       data_in.data() + pad, data_out.data() + pad, static_cast<T>(1));
    }


    
    template<typename TIndex, SignalBase Base, typename TValue>
    std::pair<TIndex, TIndex> interpolationSearch(Base data, TIndex idx1, TIndex idx2, TValue value) {
        //todo interpolate or values limits
        if (idx1 > idx2) {
            auto buffer = idx1;
            idx1 = idx2;
            idx2 = buffer;
        }
        while (std::abs(idx2 - idx1) > 1) {
            auto dx = static_cast<double>(idx2 - idx1);
            auto dy = static_cast<double>(data[idx2] - data[idx1]);
            if (dy == 0) {
                break;
            }
            auto idx_new = static_cast<TIndex>(idx1 + dx * (value - data[idx1]) / dy);
            if (data[idx_new] > value) {
                idx2 = idx_new;
            } else if (data[idx_new] < value) {
                idx1 = idx_new;
            } else if (data[idx_new] == value) {
                idx1 = idx_new;
                idx2 = idx1 + 1;
            }
        }
        return {idx1, idx2};
    }

    
    template<typename TIndex, typename TValue, typename IdxLambdaT>
    std::pair<TIndex, TIndex> interpolationSearch(TIndex idx1, TIndex idx2, TValue value, IdxLambdaT idx_lambda) {
        //todo interpolate or values limits
        if (idx1 > idx2) {
            auto buffer = idx1;
            idx1 = idx2;
            idx2 = buffer;
        }
        while (std::abs(idx2 - idx1) > 1) {
            auto dx = static_cast<double>(idx2 - idx1);
            auto dy = static_cast<double>(idx_lambda(idx2) - idx_lambda(idx1));
            if (dy == 0) {
                idx2 = idx1 + 1;
                break;
            }
            auto idx_new = static_cast<TIndex>(idx1 + dx * (value - idx_lambda(idx1)) / dy);
            if (idx_lambda(idx_new) > value) {
                idx2 = idx_new;
            } else if (idx_lambda(idx_new) < value) {
                idx1 = idx_new;
            } else if (idx_lambda(idx_new) == value) {
                idx1 = idx_new;
                idx2 = idx1 + 1;
            }
        }
        return {idx1, idx2};
    }

    
    template<typename T>
    std::pair<T, T> convertFSampleC2T(std::complex<T> sample) {
        return std::pair{
            std::sqrt(sample.real() * sample.real() + sample.imag() * sample.imag()),
            std::atan2(sample.imag(), sample.real())
        };
    }

    
    template<typename T>
    std::complex<T> convertFSampleT2C(std::pair<T, T> sample) {
        T theta = sample.second;
        T ampl = sample.first;

        auto b1 = -ampl * std::tan(theta) / std::sqrt(sec(theta) * sec(theta));
        auto a1 = -ampl / std::sqrt(sec(theta) * sec(theta));
        auto b2 = ampl * std::tan(theta) / std::sqrt(sec(theta) * sec(theta));
        auto a2 = ampl / std::sqrt(sec(theta) * sec(theta));

        std::complex<T> result1{a1, b1};
        std::complex<T> result2{a1, b2};
        std::complex<T> result3{a2, b1};
        std::complex<T> result4{a2, b2};

        T error1 = pairL2<T>(sample, convertFSampleC2T(result1));
        T error2 = pairL2<T>(sample, convertFSampleC2T(result2));
        T error3 = pairL2<T>(sample, convertFSampleC2T(result3));
        T error4 = pairL2<T>(sample, convertFSampleC2T(result4));

        if (error1 <= error2 && error1 <= error3 && error1 <= error4) {
            return result1;
        } else if (error2 <= error1 && error2 <= error3 && error2 <= error4) {
            return result2;
        } else if (error3 <= error1 && error3 <= error2 && error3 <= error4) {
            return result3;
        } else {
            return result4;
        }
    }

     
    template<typename T>
    T powFact(T a, T pow){
        double result = 0.0;
        for (int i = 0; i < a; i++)
        {
            double temp = i + 1;
            result = result + std::pow(temp, a);//temp.powf(pow);
        }
        return result;
    }

    enum class HTKind{Add, AddRepl, Mull, MullRepl};

    size_t getUniqueSpecterSamplesCount(size_t signal_size){
        return signal_size / 2 + signal_size % 2;
    }

    double getFreqByIdx(int len, int idx) {
        return static_cast<double>(idx) / len;
    }

    int getIdxByFreq(int len, double freq){
        if (freq == 0.0){
            return 0;
        }
        else{
            return static_cast<int>(freq * len);
        }
    }
    double getIdxByFreqD(int len, double freq){
        if (freq == 0.0){
            return 0;
        }
        else{
            return freq * len;
        }
    }

    struct ftResamplingData{
        int new_idx;
        int new_size;
        int old_size;
    };

    ftResamplingData getResamplingSize(int len, double freq){
        ftResamplingData result;
        result.old_size = len;
        double old_idx = getIdxByFreqD(len, freq);
        auto old_idx_i = static_cast<int>(old_idx);
        if (old_idx == old_idx_i){
            result.new_size = len;
            result.new_idx = old_idx_i;
            return result;
        }
        auto idx_new = static_cast<double>(old_idx_i + 1);
        result.new_idx = static_cast<int>(idx_new);
        result.new_size = static_cast<int>(static_cast<double>(len) / old_idx * idx_new);
        return result;
    }

    template<Signal DataT, Signal OutT, HTKind kind>
    void hilbertTransform(DataT & data, OutT & out, 
        std::vector<std::complex<double>> & specter, 
            std::vector<std::complex<double>> & buffer)
    {
        for (int i = 0; i < data.size(); i++){
            buffer[i] = {data[i], 0.0};
        }
        
        fftc2c(buffer, specter);

        auto size = data.size();
        for (int i = 1; i < data.size() / 2 + data.size() % 2; i++){
            if constexpr (kind == HTKind::Mull || kind == HTKind::MullRepl){
                specter[i] = specter[i] * 2.0;
            }
            else if constexpr (kind == HTKind::Add || kind == HTKind::AddRepl){
                specter[i] += specter[size - i];
            }
            
            specter[size - i] = specter[size - i] * 0.0;
        }

        ifftc2c(specter, buffer);

        for (int i = 0; i < data.size(); i++){
            if constexpr (kind == HTKind::AddRepl || kind == HTKind::MullRepl){
                data[i] = buffer[i].real();
            }
            out[i] = buffer[i].imag();
        }
    }

    template<Signal DataT, Signal OutT, HTKind kind>
    void hilbertTransformConst(const DataT & data, OutT & out, 
        std::vector<std::complex<double>> & specter, 
            std::vector<std::complex<double>> & buffer)
    {
        for (int i = 0; i < data.size(); i++){
            buffer[i] = {data[i], 0.0};
        }
        
        fftc2c(buffer, specter);

        auto size = data.size();
        for (int i = 1; i < data.size() / 2 + data.size() % 2; i++){
            if constexpr (kind == HTKind::Mull || kind == HTKind::MullRepl){
                specter[i] = specter[i] * 2.0;
            }
            else if constexpr (kind == HTKind::Add || kind == HTKind::AddRepl){
                specter[i] += specter[size - i];
            }
            
            specter[size - i] = specter[size - i] * 0.0;
        }

        ifftc2c(specter, buffer);

        for (int i = 0; i < data.size(); i++){
            out[i] = buffer[i].imag();
        }
    }

    template<Signal DataT, Signal OutT, Signal Weights, HTKind kind>
    void WeightedHilbertTransform(DataT & data, OutT & out, 
        std::vector<std::complex<double>> & specter, 
            std::vector<std::complex<double>> & buffer,
            const Weights & weights)
    {
        for (int i = 0; i < data.size(); i++){
            buffer[i] = {data[i], 0.0};
        }
        
        fftc2c(buffer, specter);

        auto size = data.size();
        for (int i = 1; i < data.size() / 2 + data.size() % 2; i++){
            if constexpr (kind == HTKind::Mull || kind == HTKind::MullRepl){
                specter[i] = specter[i] * (1.0 + weights[i]);
            }
            else if constexpr (kind == HTKind::Add || kind == HTKind::AddRepl){
                specter[i] += specter[size - i] * weights[i];
            }
            
            specter[size - i] = specter[size - i] * (1.0 - weights[i]);
        }

        ifftc2c(specter, buffer);

        for (int i = 0; i < data.size(); i++){
            if constexpr (kind == HTKind::AddRepl || kind == HTKind::MullRepl){
                data[i] = buffer[i].real();
            }
            out[i] = buffer[i].imag();
        }

        //todo test
    }

    template<Signal DataT, Signal OutT, HTKind kind>
    void WeightedHilbertTransformConst(const DataT & data, OutT & out, 
        std::vector<std::complex<double>> & specter, 
            std::vector<std::complex<double>> & buffer,
            const std::vector<double> & weights)
    {
        for (int i = 0; i < data.size(); i++){
            buffer[i] = {data[i], 0.0};
        }
        
        fftc2c(buffer, specter);

        auto size = data.size();
        for (int i = 1; i < data.size() / 2 + data.size() % 2; i++){
            if constexpr (kind == HTKind::Mull || kind == HTKind::MullRepl){
                specter[i] = specter[i] * (1.0 + weights[i]);
            }
            else if constexpr (kind == HTKind::Add || kind == HTKind::AddRepl){
                specter[i] += specter[size - i] * weights[i];
            }
            
            specter[size - i] = specter[size - i] * (1.0 - weights[i]);
        }

        ifftc2c(specter, buffer);

        for (int i = 0; i < data.size(); i++){
            out[i] = buffer[i].imag();
        }

        //todo test
    }
    template<typename T, Signal SignalT>
    void resampling(SignalT & data, std::vector<T> & out, size_t target_size){
        auto size = data.size();
        out.clear();
        double step = static_cast<double>(size - 1) / static_cast<double>(target_size - 1);
        for (int i = 0; i < target_size; i++){
            out.push_back(data.interpolate(static_cast<double>(i) * step, SignalKind::Universal));
        }
    }

    struct circleExtendResult{
        int freq_idx;
        int pad;
    };

    template<typename T, Signal SignalT>
    circleExtendResult circleExtend(SignalT & data, std::vector<T> & out, double target_freq){
        out.clear();
        ftResamplingData res_data =
                getResamplingSize(data.size(), target_freq);
        int pad = (res_data.new_size - data.size()) / 2;
        for (int i = 0; i < pad; i++){
            out.push_back(data[pad - i - 1]);
        }
        for (int i = pad; i < data.size() + pad; i++){
            out.push_back(data[i - pad]);
        }
        for (int i = data.size() + pad; i < res_data.new_size; i++){
            out.push_back(data[data.size() - (i - data.size() - pad + 1)]);
        }
        circleExtendResult res;
        res.freq_idx = res_data.new_idx;
        res.pad = pad;

        return res;
    }

    enum class ExtremumsKind { Simple, DerArctg };

    template<typename T, typename U>
    void computeExtremums(const T & signal, std::vector<U> & extremums, ExtremumsKind kind){
        extremums.clear();
        extremums.push_back(0);
        if (kind == ExtremumsKind::Simple){
            for (int i = 1; i < signal.size() - 1; i++) {
                if ((signal[i] >= signal[i - 1] &&
                     signal[i] > signal[i + 1]) ||
                    (signal[i] > signal[i - 1] &&
                     signal[i] >= signal[i + 1]) ||
                    (signal[i] <= signal[i - 1] &&
                     signal[i] < signal[i + 1]) ||
                    (signal[i] < signal[i - 1] &&
                     signal[i] <= signal[i + 1])) {
                    extremums.push_back(i);
                }
            }
        }
        /*else if (kind == ExtremumsKind::DerArctg){
            std::vector<double> der;
            //der.has_ovnership = true;
            for (int i = 0; i < signal.size(); i++){
                der.push_back(0.0);
            }
            GENERAL::Nil nil;
            DERIVATORS::FinniteDifference<DERIVATORS::FinniteDifferenceType::Backward> derivator;
            derivator.compute(signal, der, &nil);
            for (int i = 1; i < der.size() - 1; i++) {
                if ((der[i] >= der[i - 1] &&
                     der[i] > der[i + 1]) ||
                    (der[i] > der[i - 1] &&
                     der[i] >= der[i + 1]) ||
                    (der[i] <= der[i - 1] &&
                     der[i] < der[i + 1]) ||
                    (der[i] < der[i - 1] &&
                     der[i] <= der[i + 1])) {
                    extremums.push_back(i);
                }
            }
        }*/
        
        extremums.push_back(static_cast<int>(signal.size() - 1));
    }

    template<typename T, typename U>
    bool compareSignals(const T & signal1, const T & signal2){
        if (signal1.size() != signal2.size()){
            return false;
        }
        else{
            for (auto i = 0; i < signal1.size(); i++){
                if (signal1[i] != signal2[i]){
                    return false;
                }
            }
        }
        return true;
    }

    template<typename T, typename U>
    void assignVectors(T & vec1, const U & vec2){
        vec1.clear();
        for (int i = 0; i < vec2.size(); i++){
            vec1.push_back(vec2[i]);
        }
    }
}
