#pragma once

#include <cstdint>

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

        inline T& operator[](IdxT idx) {
            if constexpr (!std::is_same_v<DataReferenceExpression, GENERAL::Nil>) {
                return ref_expression(idx);
            } else {
                /*std::unreachable();*/
            }
        }

        inline T operator[](IdxT idx) const {
            if constexpr (!std::is_same_v<DataValueExpression, GENERAL::Nil>) {
                return (*val_expression)(idx);
            } else {
                /*std::unreachable();*/
            }
        }

        std::size_t size() const {
            if constexpr (!std::is_same_v<SizeExpression, GENERAL::Nil>) {
                return (*size_expression)();
            } else {
                /*std::unreachable();*/
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

        inline T operator[](IdxT idx) const {
            if constexpr (!std::is_same_v<DataValueExpression, GENERAL::Nil>) {
                return (*val_expression)(idx);
            } else {
                /*std::unreachable();*/
            }
        }

        std::size_t size() const {
            if constexpr (!std::is_same_v<SizeExpression, GENERAL::Nil>) {
                return (*size_expression)();
            } else {
                /*std::unreachable();*/
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
        }

        void show(PlottingKind kind, const std::string& filename, const std::string& format) const {
        }

        void show(PlottingKind kind, const std::string& filename) const {
        }




        //получение значения в неизвестной точке внутри диапазона определения (те в нашем случае по дробному индексу)
        template<typename Idx>
        SampleType interpolate(Idx idx, SignalKind kind) const {
            if (kind == SignalKind::Universal){
                int64_t integral_val = round(idx);
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
                        mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));
                        //for example with size 1000 2000 is 0; 2001 is 1; 
                    }
                    else{
                        //for example 1
                        mirror_idx = static_cast<Idx>(size_l - 1 - static_cast<int64_t>(idx) % size_l);
                        mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));
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
                            {(double)static_cast<int64_t>(idx), (*base)[static_cast<int64_t>(idx)]},
                            {(double)static_cast<int64_t>(idx + 1), (*base)[static_cast<int64_t>(idx) + 1]}, static_cast<double>(idx));
                    }
                    else{
                        return UTILITY_MATH::linearInterpolate<double, SampleType>(
                            {(double)static_cast<int64_t>(size() - 2), (*base)[size() - 2]},
                            {(double)static_cast<int64_t>(size() - 1), (*base)[size() - 1]}, static_cast<double>(idx));
                    }
                }
            }
            else if (kind == SignalKind::Monotone) {
                double left = 0.0;
                double right = size() - 1; 
                double dx = right - left;
                int64_t chunk_num;
                if (idx >= left){
                    chunk_num = static_cast<int64_t>((idx - left) / dx);
                }
                else{
                    chunk_num = static_cast<int64_t>((idx - left) / dx) - 1;
                }
                if (chunk_num == 0 || idx == right){
                    if (idx <= right - 1){
                        auto out = UTILITY_MATH::linearInterpolate<double, SampleType>(
                            {(double)static_cast<int64_t>(idx), (*base)[static_cast<int64_t>(idx)]},
                            {(double)static_cast<int64_t>(idx + 1), (*base)[static_cast<int64_t>(idx + 1)]},
                            static_cast<double>(idx)
                        );
                        return out;
                    }
                    else{
                        auto out = UTILITY_MATH::linearInterpolate<double, SampleType>(
                            {(double)static_cast<int64_t>(idx - 1), (*base)[static_cast<int64_t>(idx - 1)]},
                            {(double)static_cast<int64_t>(idx), (*base)[static_cast<int64_t>(idx)]},
                            static_cast<double>(idx)
                        );
                        return out;
                    }
                }

                auto mod = [](double a, double b) -> double{
                    return a - static_cast<double>(static_cast<int64_t>(a / b)) * b;
                };

                auto div = [](double a, double b) -> double{
                    return static_cast<double>(static_cast<int64_t>(a / b));
                };

                double top = interpolate(right, kind);
                double bot = interpolate(left, kind);
                double dy = top - bot;

                double _dx = idx - left;
                double dx_ = right - idx;

                if (chunk_num >= 0){
                    if (chunk_num % 2 == 0){
                        return bot + dy * div(_dx, dx) + interpolate(mod(_dx, dx), kind) - bot; 
                    }
                    else{
                        return bot + dy * div(_dx, dx) + dy - interpolate(right - mod(_dx, dx), kind) + bot;
                    }
                }
                else{
                    if (chunk_num % 2 == -1){
                        return bot - dy * div(std::abs(_dx), dx) - interpolate(mod(std::abs(_dx), dx), kind) + bot;
                    }
                    else{
                        return bot - dy * div(std::abs(_dx), dx) - dy + interpolate(right - mod(std::abs(_dx), dx), kind) - bot;
                    }
                }
            } else if (kind == SignalKind::Harmonic) {
                std::vector<std::complex<double>> data;
                std::vector<std::complex<double>> spectr;
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
                /*std::unreachable();*/
            }
        }

        double findMonotone(SampleType value, std::optional<double> idx1, std::optional<double> idx2, std::optional<double> start_idx, SampleType delta) const {
            if (!idx1) {
                idx1 = {0};
            }
            if (!idx2) {
                idx2 = {size() - 1};
            }
            if (*idx1 > *idx2){
                auto temp = *idx1;
                *idx1 = *idx2;
                *idx2 = temp; 
            }
            auto idx_lambda = [&](int idx) {
                return interpolate(idx, SignalKind::Monotone);
            };

            double approx_answer_x = 0.0;
            double approx_answer_y = 0.0;
            double x_left = *idx1;
            double x_right = *idx2;
            double y_left = interpolate<double>(x_left, SignalKind::Monotone);
            double y_right = interpolate<double>(x_right, SignalKind::Monotone);   

            while(std::abs(*idx1 - *idx2) > delta){
                approx_answer_x = UTILITY_MATH::backLinearInterpolate<double, double>
                    ({x_left, y_left}, {x_right, y_right}, 
                    value);
                approx_answer_y = interpolate<double>(approx_answer_x, SignalKind::Monotone);
                
                if (approx_answer_x > x_right){
                    y_left = y_right;
                    y_right = approx_answer_y;
                    x_left = x_right;
                    x_right = approx_answer_x;
                }
                else if (approx_answer_x < x_left){
                    y_right = y_left;
                    x_right = x_left;
                    y_left = approx_answer_y;
                    x_left = approx_answer_x;
                }
                else if (approx_answer_y > value){
                    y_right = approx_answer_y;
                    x_right = approx_answer_x;
                }
                else if (approx_answer_y < value){
                    y_left = approx_answer_y;
                    x_left = approx_answer_x;
                }
                if (std::abs(approx_answer_y - value) <= delta){
                    /*x_left = static_cast<int64_t>(approx_answer_x);
                    x_right = x_left + 1;
                    y_left = findPeriodDistanceAverage(phase, i, x_left);
                    y_right = findPeriodDistanceAverage(phase, i, x_right);*/
                    break;
                }
                if (y_left == y_right){
                    x_right += x_right * 0.1 + (x_right - x_left) * 0.1 + 10.0;
                    y_right = interpolate<double>(x_right, SignalKind::Monotone);
                }
            }

            return approx_answer_x;
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
        }

        void show(PlottingKind kind, const std::string& filename, const std::string& format) const {
        }

        void show(PlottingKind kind, const std::string& filename) const {
        }


        //получение значения в неизвестной точке внутри диапазона определения (те в нашем случае по дробному индексу)
        template<typename Idx>
        SampleType interpolate(Idx idx, SignalKind kind) const {
            if (kind == SignalKind::Universal){
                int64_t integral_val = round(idx);
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
                        mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));
                        //for example with size 1000 2000 is 0; 2001 is 1; 
                    }
                    else{
                        //for example 1
                        mirror_idx = static_cast<Idx>(size_l - 1 - static_cast<int64_t>(idx) % size_l);
                        mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));
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
                            {(double)static_cast<int64_t>(idx), (*base)[static_cast<int64_t>(idx)]},
                            {(double)static_cast<int64_t>(idx + 1), (*base)[static_cast<int64_t>(idx) + 1]}, static_cast<double>(idx));
                    }
                    else{
                        return UTILITY_MATH::linearInterpolate<double, SampleType>(
                            {(double)static_cast<int64_t>(size() - 2), (*base)[size() - 2]},
                            {(double)static_cast<int64_t>(size() - 1), (*base)[size() - 1]}, static_cast<double>(idx));
                    }
                }
            }
            else if (kind == SignalKind::Monotone) {
                double left = 0.0;
                double right = size() - 1; 
                double dx = right - left;
                int64_t chunk_num;
                if (idx >= left){
                    chunk_num = static_cast<int64_t>((idx - left) / dx);
                }
                else{
                    chunk_num = static_cast<int64_t>((idx - left) / dx) - 1;
                }
                if (chunk_num == 0 || idx == right){
                    if (idx <= right - 1){
                        auto out = UTILITY_MATH::linearInterpolate<double, SampleType>(
                            {(double)static_cast<int64_t>(idx), (*base)[static_cast<int64_t>(idx)]},
                            {(double)static_cast<int64_t>(idx + 1), (*base)[static_cast<int64_t>(idx + 1)]},
                            static_cast<double>(idx)
                        );
                        return out;
                    }
                    else{
                        auto out = UTILITY_MATH::linearInterpolate<double, SampleType>(
                            {(double)static_cast<int64_t>(idx - 1), (*base)[static_cast<int64_t>(idx - 1)]},
                            {(double)static_cast<int64_t>(idx), (*base)[static_cast<int64_t>(idx)]},
                            static_cast<double>(idx)
                        );
                        return out;
                    }
                }

                auto mod = [](double a, double b) -> double{
                    return a - static_cast<double>(static_cast<int64_t>(a / b)) * b;
                };

                auto div = [](double a, double b) -> double{
                    return static_cast<double>(static_cast<int64_t>(a / b));
                };

                double top = interpolate(right, kind);
                double bot = interpolate(left, kind);
                double dy = top - bot;

                double _dx = idx - left;
                double dx_ = right - idx;

                if (chunk_num >= 0){
                    if (chunk_num % 2 == 0){
                        return bot + dy * div(_dx, dx) + interpolate(mod(_dx, dx), kind) - bot; 
                    }
                    else{
                        return bot + dy * div(_dx, dx) + dy - interpolate(right - mod(_dx, dx), kind) + bot;
                    }
                }
                else{
                    if (chunk_num % 2 == -1){
                        return bot - dy * div(std::abs(_dx), dx) - interpolate(mod(std::abs(_dx), dx), kind) + bot;
                    }
                    else{
                        return bot - dy * div(std::abs(_dx), dx) - dy + interpolate(right - mod(std::abs(_dx), dx), kind) - bot;
                    }
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
                /*std::unreachable();*/
            }
        }

        double findMonotone(SampleType value, std::optional<double> idx1, std::optional<double> idx2, std::optional<double> start_idx, SampleType delta) const {
            if (!idx1) {
                idx1 = {0};
            }
            if (!idx2) {
                idx2 = {size() - 1};
            }
            if (*idx1 > *idx2){
                auto temp = *idx1;
                *idx1 = *idx2;
                *idx2 = temp; 
            }
            auto idx_lambda = [&](int idx) {
                return interpolate(idx, SignalKind::Monotone);
            };

            double approx_answer_x = 0.0;
            double approx_answer_y = 0.0;
            double x_left = *idx1;
            double x_right = *idx2;
            double y_left = interpolate<double>(x_left, SignalKind::Monotone);
            double y_right = interpolate<double>(x_right, SignalKind::Monotone);  
            while(std::abs(*idx1 - *idx2) > delta){
                approx_answer_x = UTILITY_MATH::backLinearInterpolate<double, double>
                    ({x_left, y_left}, {x_right, y_right}, 
                    value);
                approx_answer_y = interpolate<double>(approx_answer_x, SignalKind::Monotone);
                
                if (approx_answer_x > x_right){
                    y_left = y_right;
                    y_right = approx_answer_y;
                    x_left = x_right;
                    x_right = approx_answer_x;
                }
                else if (approx_answer_x < x_left){
                    y_right = y_left;
                    x_right = x_left;
                    y_left = approx_answer_y;
                    x_left = approx_answer_x;
                }
                else if (approx_answer_y > value){
                    y_right = approx_answer_y;
                    x_right = approx_answer_x;
                }
                else if (approx_answer_y < value){
                    y_left = approx_answer_y;
                    x_left = approx_answer_x;
                }
                if (std::abs(approx_answer_y - value) <= delta){
                    /*x_left = static_cast<int64_t>(approx_answer_x);
                    x_right = x_left + 1;
                    y_left = findPeriodDistanceAverage(phase, i, x_left);
                    y_right = findPeriodDistanceAverage(phase, i, x_right);*/
                    break;
                }
                if (y_left == y_right){
                    x_right += x_right * 0.1 + (x_right - x_left) * 0.1 + 10.0;
                    y_right = interpolate<double>(x_right, SignalKind::Monotone);
                }
            }

            return approx_answer_x;
        }
    };

    template<SignalBase BaseT>
    struct GenericSignalRExpr {
        bool has_ovnership = false;
        constexpr static bool is_writable = false;
        constexpr static bool is_signal = true;
        using Base = BaseT;
        using IdxType = typename Base::IdxType;
        using SampleType = typename Base::SampleType;
        Base* base;

        GenericSignalRExpr(Base& base_o) {
            base = &base_o;
        }

        GenericSignalRExpr() {
            base = new Base;
            has_ovnership = true;
        }

        ~GenericSignalRExpr() {
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
        }

        void show(PlottingKind kind, const std::string& filename, const std::string& format) const {
        }

        void show(PlottingKind kind, const std::string& filename) const {
        }


        //получение значения в неизвестной точке внутри диапазона определения (те в нашем случае по дробному индексу)
        template<typename Idx>
        SampleType interpolate(Idx idx, SignalKind kind) const {
            if (kind == SignalKind::Universal){
                int64_t integral_val = round(idx);
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
                        mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));
                        //for example with size 1000 2000 is 0; 2001 is 1; 
                    }
                    else{
                        //for example 1
                        mirror_idx = static_cast<Idx>(size_l - 1 - static_cast<int64_t>(idx) % size_l);
                        mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));
                        //for example with size 1000 1000 is 999, 1001 is 998
                    }
                    
                    mirror_idx = mirror_idx + idx - static_cast<Idx>(static_cast<int64_t>(idx));

                    SampleType mirror_val = interpolate<Idx>(mirror_idx, SignalKind::Universal);
                    SampleType dfull = (*this)[size_l-1] - (*this)[0];
                    SampleType d_mirror = (*this)[size_l - 1] - mirror_val;
                    return (*this)[0] + tile_number * dfull + d_mirror;
                }
                else{
                    return (*base)[idx];
                }
            }
            else if (kind == SignalKind::Monotone) {
                double left = 0.0;
                double right = size() - 1; 
                double dx = right - left;
                int64_t chunk_num;
                if (idx >= left){
                    chunk_num = static_cast<int64_t>((idx - left) / dx);
                }
                else{
                    chunk_num = static_cast<int64_t>((idx - left) / dx) - 1;
                }
                if (chunk_num == 0 || idx == right){
                    return (*base)[idx];
                }

                auto mod = [](double a, double b) -> double{
                    return a - static_cast<double>(static_cast<int64_t>(a / b)) * b;
                };

                auto div = [](double a, double b) -> double{
                    return static_cast<double>(static_cast<int64_t>(a / b));
                };

                double top = interpolate(right, kind);
                double bot = interpolate(left, kind);
                double dy = top - bot;

                double _dx = idx - left;
                double dx_ = right - idx;

                if (chunk_num >= 0){
                    if (chunk_num % 2 == 0){
                        return bot + dy * div(_dx, dx) + interpolate(mod(_dx, dx), kind) - bot; 
                    }
                    else{
                        return bot + dy * div(_dx, dx) + dy - interpolate(right - mod(_dx, dx), kind) + bot;
                    }
                }
                else{
                    if (chunk_num % 2 == -1){
                        return bot - dy * div(std::abs(_dx), dx) - interpolate(mod(std::abs(_dx), dx), kind) + bot;
                    }
                    else{
                        return bot - dy * div(std::abs(_dx), dx) - dy + interpolate(right - mod(std::abs(_dx), dx), kind) - bot;
                    }
                }
            } else if (kind == SignalKind::Harmonic) {
                std::vector<std::complex<double>> data;
                std::vector<std::complex<double>> spectr;
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
                /*std::unreachable();*/
            }
        }

        double findMonotone(SampleType value, std::optional<double> idx1, std::optional<double> idx2, std::optional<double> start_idx, SampleType delta) const {
            if (!idx1) {
                idx1 = {0};
            }
            if (!idx2) {
                idx2 = {size() - 1};
            }
            if (*idx1 > *idx2){
                auto temp = *idx1;
                *idx1 = *idx2;
                *idx2 = temp; 
            }
            auto idx_lambda = [&](double idx) {
                return interpolate(idx, SignalKind::Monotone);
            };

            double approx_answer_x = 0.0;
            double approx_answer_y = 0.0;
            double x_left = *idx1;
            double x_right = *idx2;
            double y_left = interpolate<double>(x_left, SignalKind::Monotone);
            double y_right = interpolate<double>(x_right, SignalKind::Monotone);   

            while(std::abs(*idx1 - *idx2) > delta){
                approx_answer_x = UTILITY_MATH::backLinearInterpolate<double, double>
                    ({x_left, y_left}, {x_right, y_right}, 
                    value);
                approx_answer_y = interpolate<double>(approx_answer_x, SignalKind::Monotone);
                
                if (approx_answer_x > x_right){
                    y_left = y_right;
                    y_right = approx_answer_y;
                    x_left = x_right;
                    x_right = approx_answer_x;
                }
                else if (approx_answer_x < x_left){
                    y_right = y_left;
                    x_right = x_left;
                    y_left = approx_answer_y;
                    x_left = approx_answer_x;
                }
                else if (approx_answer_y > value){
                    y_right = approx_answer_y;
                    x_right = approx_answer_x;
                }
                else if (approx_answer_y < value){
                    y_left = approx_answer_y;
                    x_left = approx_answer_x;
                }
                if (std::abs(approx_answer_y - value) <= delta){
                    /*x_left = static_cast<int64_t>(approx_answer_x);
                    x_right = x_left + 1;
                    y_left = findPeriodDistanceAverage(phase, i, x_left);
                    y_right = findPeriodDistanceAverage(phase, i, x_right);*/
                    break;
                }
                if (y_left == y_right || x_right == x_left){
                    x_right += std::abs(x_right * 0.1 + (x_right - x_left)) * 0.1 + 10.0;
                    //x_left -= std::abs(x_left * 0.1 + (x_right - x_left) * 0.1 + 10.0);
                    y_right = interpolate<double>(x_right, SignalKind::Monotone);
                    //y_left = interpolate<double>(x_left, SignalKind::Monotone);
                }
            }

            return approx_answer_x;
        }
    };


    static_assert(is_signal<GenericSignal<SimpleVecWrapper<int>, true>>);
}
