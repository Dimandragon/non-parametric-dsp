module;

#include "pocketfft_hdronly.h"
#include "icecream.hpp"

export module utility_math;

import <utility>;
import <complex>;
import <vector>;

import npdsp_concepts;
import npdsp_config;

namespace NP_DSP::ONE_D::UTILITY_MATH {
    export
    template<typename T>
    T complexL2(std::complex<T> a, std::complex<T> b) {
        return (a.imag() - b.imag()) * (a.imag() - b.imag()) + (a.real() - b.real()) * (a.real() - b.real());
    }

    export
    template<typename T>
    T pairL2(std::pair<T, T> a, std::pair<T, T> b) {
        return (a.first - b.first) * (a.first - b.first) + (a.second - b.second) * (a.second - b.second);
    }

    auto sec(auto z_r) {
        return 1 / std::cos(z_r);
    }

    export
    template<typename T, typename data1T, typename data2T>
    T signalsL2Distance(const data1T& data1, const data2T& data2) {
        double error = 0.0;
        error = 0.0;
        for (int i = 0; i < data1.size(); i++) {
            error = error + std::sqrt((data1[i] - data2[i]) * (data1[i] - data2[i]));
        }
        return error;
    }

    export
    template<typename xType, typename yType>
    yType linearInterpolate(std::pair<xType, yType> point1, std::pair<xType, yType> point2, xType x_in) {
        auto dx = point2.first - point1.first;
        auto dy = point2.second - point1.second;
        if (dx == 0) {
            return (point1.second + point2.second) / 2.0;
        }
        return point1.second + dy * (x_in - point1.first) / dx;
    }

    export
    template<typename xType, typename yType>
    xType backLinearInterpolate(std::pair<xType, yType> point1, std::pair<xType, yType> point2, yType y_in) {
        auto dx = point2.first - point1.first;
        auto dy = point2.second - point1.second;
        if (dy == 0) {
            return (point1.first + point2.first) / 2;
        }
        return point1.first + dx * (y_in - point1.second) / dy;
    }

    export
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

    export
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

    export
    template<typename T>
    void fftc2c(const Signal<std::complex<T>>& in, Signal<std::complex<T>>& out) {
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

    export
    template<typename T>
    void ifftc2c(const Signal<std::complex<T>>& in, Signal<std::complex<T>>& out) {
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

    export
    template<typename T>
    void fastConvolution(Signal<T>& in1, Signal<T>& in2, Signal<T>& out) {
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
            sp1[i] = sp1[i] * sp2[i] * static_cast<double>(len);
        }

        ifftc2c(sp1, data_out);

        for (int i = 0; i < len; i++) {
            out[i] = data_out[i].real();
        }
    }

    export
    template<typename T>
    void fftc2c(std::vector<std::complex<T>> const& data_in, std::vector<std::complex<T>>& data_out, size_t len,
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
        pocketfft::c2c(shape, stridef, stridef, axes, pocketfft::FORWARD,
                       data_in.data() + pad, data_out.data() + pad, ((T) 1.0) / ((T) len));
    }

    export
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

    /*
    export
    template <typename TIndex, SignalBase Base, typename TValue>
    std::pair<TIndex, TIndex> interpoationSearch(Base data, TIndex idx1, TIndex idx2, TValue value) {
        //todo interpolate or values limits
        if (idx1 > idx2) {
            auto buffer = idx1;
            idx1 = idx2;
            idx2 = buffer;
        }
        while (std::abs(idx2 - idx1) > 1) {
            auto dx = static_cast<double>(idx2 - idx1);
            auto dy = static_cast<double>(data[idx2] - data[idx1]);
            if (dy == 0){
                break;
            }
            auto idx_new = static_cast<TIndex>(idx1 + dx * (value - data[idx1])/dy);
            if constexpr (CONFIG::debug){
                IC(idx1, idx2, dx, dy, idx_new);
            }
            if (data[idx_new] > value) {
                idx2 = idx_new;
            }
            else if(data[idx_new] < value){
                idx1 = idx_new;
            }
            else if (data[idx_new] == value){
                idx1 = idx_new;
                idx2 = idx1 + 1;
            }
            if constexpr (CONFIG::debug){
                IC(idx1, idx2);
            }
        }
        return {idx1, idx2};
    }*/

    /*export
    template <typename TIndex, typename TValue, typename IdxLambdaT>
    std::pair<TIndex, TIndex> interpoationSearch(TIndex idx1, TIndex idx2, TValue value, IdxLambdaT idx_lambda) {
        //todo interpolate or values limits
        if (idx1 > idx2) {
            auto buffer = idx1;
            idx1 = idx2;
            idx2 = buffer;
        }
        while (std::abs(idx2 - idx1) > 1) {
            auto dx = static_cast<double>(idx2 - idx1);
            auto dy = static_cast<double>(idx_lambda(idx2) - idx_lambda(idx1));
            if (dy == 0){
                idx2 = idx1 + 1;
                break;
            }
            auto idx_new = static_cast<TIndex>(idx1 + dx * (value - idx_lambda(idx1))/dy);
            if constexpr (CONFIG::debug){
                std::string mark = "creating idx_new in interpolation search";
                IC(idx1, idx2, dx, dy, idx_new);
            }
            if (idx_lambda(idx_new) > value) {
                idx2 = idx_new;
            }
            else if(idx_lambda(idx_new) < value){
                idx1 = idx_new;
            }
            else if (idx_lambda(idx_new) == value){
                idx1 = idx_new;
                idx2 = idx1 + 1;
            }
            if constexpr (CONFIG::debug){
                IC(idx1, idx2);
            }
        }
        return {idx1, idx2};
    }*/

    export
    template<typename T>
    std::pair<T, T> convertFSampleC2T(std::complex<T> sample) {
        return std::pair{
            std::sqrt(sample.real() * sample.real() + sample.imag() * sample.imag()),
            std::atan2(sample.imag(), sample.real())
        };
    }

    export
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
}
