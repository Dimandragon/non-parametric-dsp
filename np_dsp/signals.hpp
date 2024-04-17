#pragma once

#include "icecream.hpp"

#include <matplot/matplot.h>

#include <vector>
#include <cassert>
#include <cstddef>
#include <utility>
#include <type_traits>
#include <optional>
#include <numbers>
#include <complex>

#include <npdsp_concepts.hpp>
#include <utility_math.hpp>
#include <npdsp_config.hpp>

namespace NP_DSP::ONE_D {
    
    template<typename T>
    struct SimpleVecWrapper {
        bool has_ovnership = false;

        using SampleType = T;
        using IdxType = std::size_t;
        constexpr static bool is_signal_base = true;
        constexpr static bool is_writable = true;

        std::vector<T>* vec;

        SimpleVecWrapper() {
            has_ovnership = true;
            vec = new std::vector<T>;
        }

        ~SimpleVecWrapper() {
            if (has_ovnership) {
                delete vec;
            }
        }

        SimpleVecWrapper(std::vector<T>& vec_in) {
            vec = &vec_in;
        }

        inline T& operator[](std::size_t idx) {
            return (*vec)[idx];
        }

        inline const T& operator[](std::size_t idx) const {
            return (*vec)[idx];
        }

        inline std::size_t size() const {
            return vec->size();
        }
    };

    static_assert(is_signal_base_first<SimpleVecWrapper<int>>);

    
    template<typename T, typename IdxT, typename DataValExpr,
        typename DataRefExpr, typename SizeExpr, bool is_writeble_b>
    struct ExpressionWrapper {
        using SampleType = T;
        using IdxType = IdxT;
        using DataValueExpression = DataValExpr;
        using DataReferenceExpression = DataRefExpr;
        using SizeExpression = SizeExpr;

        constexpr static bool is_signal_base = true;
        constexpr static bool is_writable = is_writeble_b;

        DataValExpr* val_expression;
        DataRefExpr* ref_expression;
        SizeExpr* size_expression;

        ExpressionWrapper(DataValExpr& val_expr, DataRefExpr& ref_expr, SizeExpr& size_expr) {
            val_expression = &val_expr;
            ref_expression = &ref_expr;
            size_expression = &size_expr;
        }

        inline T& operator[](std::size_t idx) {
            if constexpr (!std::is_same_v<DataReferenceExpression, GENERAL::Nil>) {
                return ref_expression(idx);
            } else {
                std::unreachable();
            }
        }

        inline T operator[](std::size_t idx) const {
            if constexpr (!std::is_same_v<DataValueExpression, GENERAL::Nil>) {
                return (*val_expression)(idx);
            } else {
                std::unreachable();
            }
        }

        std::size_t size() const {
            if constexpr (!std::is_same_v<SizeExpression, GENERAL::Nil>) {
                return (*size_expression)();
            } else {
                std::unreachable();
            }
        }
    };

    
    template<typename T, typename IdxT, typename DataValExpr,
        typename SizeExpr>
    struct ExpressionWrapper<T, IdxT, DataValExpr,
                GENERAL::Nil, SizeExpr, false> {
        using SampleType = T;
        using IdxType = IdxT;
        using DataValueExpression = DataValExpr;
        using DataReferenceExpression = GENERAL::Nil;
        using SizeExpression = SizeExpr;

        constexpr static bool is_signal_base = true;
        constexpr static bool is_writable = false;

        DataValExpr* val_expression;
        SizeExpr* size_expression;

        ExpressionWrapper(DataValExpr& val_expr, SizeExpr& size_expr) {
            val_expression = &val_expr;
            size_expression = &size_expr;
        }

        inline T operator[](std::size_t idx) const {
            if constexpr (!std::is_same_v<DataValueExpression, GENERAL::Nil>) {
                return (*val_expression)(idx);
            } else {
                std::unreachable();
            }
        }

        std::size_t size() const {
            if constexpr (!std::is_same_v<SizeExpression, GENERAL::Nil>) {
                return (*size_expression)();
            } else {
                std::unreachable();
            }
        }
    };

    // enum class SignalKind {Monotone, Stohastic, Harmonic, Smooth};

    
    template<SignalBase BaseT, bool is_writable_b>
    struct GenericSignal {
        bool has_ovnership = false;
        constexpr static bool is_writable = is_writable_b;
        constexpr static bool is_signal = true;
        using Base = BaseT;
        using IdxType = typename Base::IdxType;
        using SampleType = typename Base::SampleType;
        Base* base;

        GenericSignal(Base& base_o) {
            base = &base_o;
        }

        GenericSignal() {
            base = new Base;
            has_ovnership = true;
        }

        ~GenericSignal() {
            if (has_ovnership) {
                delete base;
            }
        }

        inline SampleType& operator[](IdxType idx) {
            return (*base)[idx];
        }

        inline SampleType operator[](IdxType idx) const {
            return (*base)[idx];
        }

        inline size_t size() const {
            return base->size();
        }

        void show(PlottingKind kind) const {
            if (kind == PlottingKind::Interpolate) {
                std::vector<SampleType> plotting_data = {};
                int i = -size();
                while (i < static_cast<int>(size() * 2)) {
                    ++i;
                    plotting_data.push_back(static_cast<float>(interpolate<double>(static_cast<double>(i), SignalKind::Harmonic)));
                }
                matplot::plot(plotting_data);
                matplot::show();
            } else if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < base->size(); ++i) {
                    auto sample = (*base)[i];
                    plotting_data.push_back(sample);
                }
                matplot::plot(plotting_data);
                matplot::show();
            }
            else if (kind == PlottingKind::Spectre){
                std::vector<std::complex<double>> spectre;
                std::vector<std::complex<double>> fft_vec;

                for (int i = 0; i < base->size(); i++){
                    spectre.push_back({0.0, 0.0});
                    fft_vec.push_back((*base)[i]);
                }

                UTILITY_MATH::fftc2c(fft_vec, spectre);

                std::vector<double> plotting_vec_re;
                std::vector<double> plotting_vec_im;

                for (int i = 0; i < base->size(); i++){
                    plotting_vec_re.push_back(spectre[i].real());
                    plotting_vec_im.push_back(spectre[i].imag());
                }
                matplot::plot(plotting_vec_re);
                matplot::hold(true);
                matplot::plot(plotting_vec_im);
                matplot::hold(false);
                matplot::show();
            }
        }

        void show(PlottingKind kind, const std::string& filename, const std::string& format) const {
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < base->size(); ++i) {
                    plotting_data.push_back((*base)[i]);
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename, format);
            } else if (kind == PlottingKind::Interpolate) {
                std::vector<SampleType> plotting_data = {};
                int i = -size();
                while (i < static_cast<int>(size() * 2)) {
                    ++i;
                    plotting_data.push_back(interpolate<int>(i, SignalKind::Stohastic));
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename, format);
            }
        }

        void show(PlottingKind kind, const std::string& filename) const {
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < base->size(); ++i) {
                    plotting_data.push_back((*base)[i]);
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename);
            } else if (kind == PlottingKind::Interpolate) {
                std::vector<SampleType> plotting_data = {};
                int i = -size();
                while (i < static_cast<int>(size() * 2)) {
                    ++i;
                    plotting_data.push_back(interpolate<int>(i, SignalKind::Stohastic));
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename);
            }

        }




        //получение значения в неизвестной точке внутри диапазона определения (те в нашем случае по дробному индексу)
        template<typename Idx>
        SampleType interpolate(Idx idx, SignalKind kind) const {

            if (kind == SignalKind::Universal){
                int64_t size_l = static_cast<int64_t>(size());
                int64_t tile_number = static_cast<int64_t>(idx) / size_l;
                //[size] is begin of 1 tile
                //[-size] is end of -1 tile and [size] is begin of -2 tile
                if (idx <= -1){
                    if (static_cast<int64_t>(idx) % size_l != size_l){
                        tile_number -= 1;
                    }
                    Idx mirror_idx;
                    if (tile_number % 2 == 0){
                        mirror_idx = static_cast<Idx>(size_l - 1 - static_cast<int64_t>(idx + 1) % size_l);
                        //for example with size 1000 -1001 is 999, -2000 is 0
                    }
                    else{
                        //-1 for example
                        mirror_idx = static_cast<Idx>(-static_cast<int64_t>(idx + 1) % size_l);
                        //for example with size 1000 -3000 is 999; -2001 is 0; 
                    }
                    mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));
                    SampleType mirror_val = interpolate<Idx>(mirror_idx, SignalKind::Universal);
                    SampleType dfull = (*this)[size_l-1] - (*this)[0];
                    SampleType d_mirror = mirror_val - (*this)[0];
                    return (tile_number + 1) * dfull - d_mirror + (*this)[0];
                }
                else if (idx > size_l){
                    Idx mirror_idx;
                    if (tile_number % 2 == 0){
                        //for example 2
                        mirror_idx = static_cast<Idx>(static_cast<int64_t>(idx) % size_l);
                        //for example with size 1000 2000 is 0; 2001 is 1; 
                    }
                    else{
                        //for example 1
                        mirror_idx = static_cast<Idx>(size_l - 1 - static_cast<int64_t>(idx) % size_l);
                        //for example with size 1000 1000 is 999, 1001 is 998
                    }
                    mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));
                    SampleType mirror_val = interpolate<Idx>(mirror_idx, SignalKind::Universal);
                    SampleType dfull = (*this)[size_l-1] - (*this)[0];
                    SampleType d_mirror = (*this)[size_l - 1] - mirror_val;
                    return (*this)[0] + tile_number * dfull + d_mirror;
                }
                else{
                    if (idx <= size_l - 2){
                        return UTILITY_MATH::linearInterpolate<double, SampleType>(
                            {static_cast<int64_t>(idx), (*base)[static_cast<int64_t>(idx)]},
                            {static_cast<int64_t>(idx + 1), (*base)[static_cast<int64_t>(idx) + 1]}, static_cast<double>(idx));
                    }
                    else{
                        return UTILITY_MATH::linearInterpolate<double, SampleType>(
                            {static_cast<int64_t>(size() - 2), (*base)[size() - 2]},
                            {static_cast<int64_t>(size() - 1), (*base)[size() - 1]}, static_cast<double>(idx));
                    }
                }
            }
            else if (kind == SignalKind::Monotone) {
                if (idx >= 0 && idx < static_cast<Idx>(size() - 2)) {
                    if constexpr (CONFIG::debug) {
                        std::string mark = "1b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(idx), (*base)[static_cast<int>(idx)]},
                        {static_cast<double>(idx + 1), (*base)[idx + 1]}, static_cast<double>(idx));
                } else if (idx < 0) {
                    Idx idx_new = idx + ((0 - static_cast<int>(idx)) % static_cast<int>(size())) * 2;
                    if constexpr (CONFIG::debug) {
                        std::string mark = "2b";
                        IC(mark);
                        IC(idx_new);
                    }
                    if (idx_new == idx) {
                        idx_new++;
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(idx_new), value_new}, {0.0, (*base)[0]}, static_cast<double>(idx));
                } else if (idx > size() - 1) {
                    Idx idx_new = idx - ((static_cast<int>(idx) - static_cast<int>(size())) % static_cast<int>(size()))
                                  * 2 - 1;
                    if (idx_new == idx) {
                        idx_new--;
                    }
                    if (idx_new > size() - 1) {
                        idx_new = size() - 2;
                    }
                    if constexpr (CONFIG::debug) {
                        std::string mark = "3b";
                        IC(mark);
                        IC(idx_new);
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(size() - 1), (*base)[size() - 1]},
                        {static_cast<double>(idx_new), value_new}, static_cast<double>(idx));
                } else {
                    if constexpr (CONFIG::debug) {
                        std::string mark = "4b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(size() - 1), (*base)[size() - 1]},
                        {static_cast<double>(size() - 2), (*base)[size() - 2]}, static_cast<double>(idx));
                }
            } else if (kind == SignalKind::Harmonic) {
                std::vector<std::complex<double>> data;
                std::vector<std::complex<double>> spectr;
                //using DT = GenericSignal<SimpleVecWrapper<std::complex<double>>, true>;
                if (spectr.size() != size()) {
                    spectr.clear();
                    data.clear();
                    for (int i = 0; i < size(); i++) {
                        spectr.push_back({0.0, 0.0});
                        data.push_back({(*this)[i], 0.0});
                    }
                } else {
                    for (int i = 0; i < size(); i++) {
                        data[i] = {(*this)[i], 0.0};
                        spectr[i] = {0.0, 0.0};
                    }
                }
                UTILITY_MATH::fftc2c(data, spectr);
                //IC()
                //spectr.show(PlottingKind::Simple);


                std::complex<double> accum = {0.0, 0.0};
                for (int i = 0; i < size(); i++) {
                    std::complex b = {
                        std::cos(2.0 * std::numbers::pi * i * idx / size()),
                        std::sin(2.0 * std::numbers::pi * i * idx / size())
                    };
                    accum += spectr[i] * b;
                }
                return accum.real();
            } else if (kind == SignalKind::Stohastic) {
                if (idx >= 0 && idx < static_cast<Idx>(size() - 2)) {
                    if constexpr (NP_DSP::CONFIG::debug) {
                        std::string mark = "1b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(idx), (*base)[static_cast<int>(idx)]},
                        {static_cast<double>(idx + 1), (*base)[idx + 1]}, static_cast<double>(idx));
                } else if (idx < 0) {
                    Idx idx_new = 0.0 - idx;

                    if (idx_new == idx) {
                        ++idx_new;
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return value_new;
                } else if (idx > size() - 1) {
                    Idx idx_new = 2.0 * size() - 2.0 - idx;

                    SampleType value_new = interpolate(idx_new, kind);
                    return value_new;
                } else {
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(size() - 1), (*base)[size() - 1]},
                        {static_cast<double>(size() - 2), (*base)[size() - 2]}, static_cast<double>(idx));
                }
            } else if (kind == SignalKind::Smooth) {
                //todo Teilors Series
                std::unreachable();
            }
        }

        template<typename Idx>
        IdxType findMonotone(SampleType value, std::optional<Idx> idx1, std::optional<Idx> idx2) const {
            if (!idx1) {
                if (*idx2 != 0) {
                    idx1 = {0};
                } else {
                    idx1 = {-size()};
                }
            }
            if (!idx2) {
                if (*idx2 != 0) {
                    idx1 = {2 * size()};
                } else {
                    idx2 = {size()};
                }
            }
            auto idx_lambda = [&](int idx) {
                return interpolate(idx, SignalKind::Monotone);
            };
            std::pair<int, int> idxes = ONE_D::UTILITY_MATH::interpolationSearch(*idx1, *idx2, value, idx_lambda);
            if constexpr (CONFIG::debug) {
                std::string mark = "find monotone";
                IC(mark, idxes.first, idxes.second, value);
            }
            return UTILITY_MATH::backLinearInterpolate<Idx, SampleType>(
                {*idx1, (*base)[*idx1]}, {*idx2, (*base)[*idx2]}, value);
        }
    };

    
    template<SignalBase BaseT>
    struct GenericSignal<BaseT, false> {
        bool has_ovnership = false;
        constexpr static bool is_writable = false;
        constexpr static bool is_signal = true;
        using Base = BaseT;
        using IdxType = typename Base::IdxType;
        using SampleType = typename Base::SampleType;
        Base* base;

        GenericSignal(Base& base_o) {
            base = &base_o;
        }

        GenericSignal() {
            base = new Base;
            has_ovnership = true;
        }

        ~GenericSignal() {
            if (has_ovnership) {
                delete base;
            }
        }

        inline SampleType operator[](IdxType idx) const {
            return (*base)[idx];
        }

        inline size_t size() const {
            return base->size();
        }

        void show(PlottingKind kind) {
            if (kind == PlottingKind::Interpolate) {
                std::vector<SampleType> plotting_data = {};
                int i = -size();
                while (i < static_cast<int>(size() * 2)) {
                    i++;
                    plotting_data.push_back(static_cast<float>(interpolate<double>(static_cast<double>(i), SignalKind::Stohastic)));
                }
                matplot::plot(plotting_data);
                matplot::show();
            } else if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < base->size(); i++) {
                    auto sample = (*base)[i];
                    plotting_data.push_back(sample);
                }
                matplot::plot(plotting_data);
                matplot::show();
            }
        }

        void show(PlottingKind kind, const std::string& filename, const std::string& format) const {
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < base->size(); i++) {
                    plotting_data.push_back((*base)[i]);
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename, format);
            } else if (kind == PlottingKind::Interpolate) {
                std::vector<SampleType> plotting_data = {};
                int i = -size();
                while (i < static_cast<int>(size() * 2)) {
                    i++;
                    plotting_data.push_back(interpolate<int>(i, SignalKind::Stohastic));
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename, format);
            }
        }

        void show(PlottingKind kind, const std::string& filename) const {
            if (kind == PlottingKind::Simple) {
                std::vector<SampleType> plotting_data = {};
                for (auto i = 0; i < base->size(); i++) {
                    plotting_data.push_back((*base)[i]);
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename);
            } else if (kind == PlottingKind::Interpolate) {
                std::vector<SampleType> plotting_data = {};
                int i = -size();
                while (i < static_cast<int>(size() * 2)) {
                    i++;
                    plotting_data.push_back(interpolate<int>(i, SignalKind::Stohastic));
                }
                matplot::plot(plotting_data);
                matplot::show();
                matplot::save(filename);
            }
        }


        //получение значения в неизвестной точке внутри диапазона определения (те в нашем случае по дробному индексу)
        template<typename Idx>
        SampleType interpolate(Idx idx, SignalKind kind) const {
            if (kind == SignalKind::Monotone) {
                if (idx >= 0 && idx < static_cast<Idx>(size() - 2)) {
                    if constexpr (CONFIG::debug) {
                        std::string mark = "1b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(idx), (*base)[static_cast<int>(idx)]},
                        {static_cast<double>(idx + 1), (*base)[idx + 1]}, static_cast<double>(idx));
                } else if (idx < 0) {
                    Idx idx_new = idx + ((0 - static_cast<int>(idx)) % static_cast<int>(size())) * 2;
                    if constexpr (CONFIG::debug) {
                        std::string mark = "2b";
                        IC(mark);
                        IC(idx_new);
                    }
                    if (idx_new == idx) {
                        idx_new++;
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(idx_new), value_new}, {0.0, (*base)[0]}, static_cast<double>(idx));
                } else if (idx > size() - 1) {
                    Idx idx_new = idx - ((static_cast<int>(idx) - static_cast<int>(size())) % static_cast<int>(size()))
                                  * 2 - 1;
                    if (idx_new == idx) {
                        idx_new--;
                    }
                    if (idx_new > size() - 1) {
                        idx_new = size() - 2;
                    }
                    if constexpr (CONFIG::debug) {
                        std::string mark = "3b";
                        IC(mark);
                        IC(idx_new);
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(size() - 1), (*base)[size() - 1]},
                        {static_cast<double>(idx_new), value_new}, static_cast<double>(idx));
                } else {
                    if constexpr (CONFIG::debug) {
                        std::string mark = "4b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(size() - 1), (*base)[size() - 1]},
                        {static_cast<double>(size() - 2), (*base)[size() - 2]}, static_cast<double>(idx));
                }
            } else if (kind == SignalKind::Harmonic) {
                std::vector<std::complex<double>> data;
                std::vector<std::complex<double>> spectr;
                //using DT = GenericSignal<SimpleVecWrapper<std::complex<double>>, true>;
                if (spectr.size() != size()) {
                    spectr.clear();
                    data.clear();
                    for (int i = 0; i < size(); i++) {
                        spectr.push_back({0.0, 0.0});
                        data.push_back({(*this)[i], 0.0});
                    }
                } else {
                    for (int i = 0; i < size(); i++) {
                        data[i] = {(*this)[i], 0.0};
                        spectr[i] = {0.0, 0.0};
                    }
                }
                UTILITY_MATH::fftc2c(data, spectr);
                //IC()
                //spectr.show(PlottingKind::Simple);


                std::complex<double> accum = {0.0, 0.0};
                for (int i = 0; i < size(); i++) {
                    std::complex b = {
                        std::cos(2.0 * std::numbers::pi * i * idx / size()),
                        std::sin(2.0 * std::numbers::pi * i * idx / size())
                    };
                    accum += spectr[i] * b;
                }
                return accum.real();
            } else if (kind == SignalKind::Stohastic) {
                if (idx >= 0 && idx < static_cast<Idx>(size() - 2)) {
                    if constexpr (NP_DSP::CONFIG::debug) {
                        std::string mark = "1b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(idx), (*base)[static_cast<int>(idx)]},
                        {static_cast<double>(idx + 1), (*base)[idx + 1]}, static_cast<double>(idx));
                } else if (idx < 0) {
                    Idx idx_new = 0.0 - idx;

                    if (idx_new == idx) {
                        ++idx_new;
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return value_new;
                } else if (idx > size() - 1) {
                    Idx idx_new = 2.0 * size() - 2.0 - idx;

                    SampleType value_new = interpolate(idx_new, kind);
                    return value_new;
                } else {
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(size() - 1), (*base)[size() - 1]},
                        {static_cast<double>(size() - 2), (*base)[size() - 2]}, static_cast<double>(idx));
                }
            } else if (kind == SignalKind::Smooth) {
                //todo Teilors Series
                std::unreachable();
            }
        }

        template<typename Idx>
        IdxType findMonotone(SampleType value, std::optional<Idx> idx1, std::optional<Idx> idx2) const {
            if (!idx1) {
                if (*idx2 != 0) {
                    idx1 = {0};
                } else {
                    idx1 = {-size()};
                }
            }
            if (!idx2) {
                if (*idx2 != 0) {
                    idx1 = {2 * size()};
                } else {
                    idx2 = {size()};
                }
            }
            auto idx_lambda = [&](int idx) {
                return interpolate(idx, SignalKind::Monotone);
            };
            std::pair<int, int> idxes = UTILITY_MATH::interpolationSearch(*idx1, *idx2, value, idx_lambda);
            if constexpr (CONFIG::debug) {
                std::string mark = "find monotone";
                IC(mark, idxes.first, idxes.second, value);
            }
            return UTILITY_MATH::backLinearInterpolate<Idx, SampleType>(
                {*idx1, (*base)[*idx1]}, {*idx2, (*base)[*idx2]}, value);
        }
    };

    static_assert(ONE_D::is_signal<GenericSignal<SimpleVecWrapper<int>, true>>);
}
