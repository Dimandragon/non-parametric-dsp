module;

#include "icecream.hpp"

export module signals;

import <matplot/matplot.h>;

import <vector>;
import <tuple>;
import <cassert>;
import <cstddef>;
import <utility>;
import <type_traits>;
import <optional>;
import <memory>;
import <iostream>;
import <complex>;

import npdsp_concepts;
import utility_math;
import npdsp_config;

namespace NP_DSP
{
    namespace ONE_D{
        export
        template<typename T>
        struct SimpleVecWrapper{
            bool has_ovnership = false;

            using SampleType = T;
            using IdxType = std::size_t;
            constexpr static bool is_signal_base = true;
            constexpr static bool is_writable = true;

            std::vector<T> * vec;

            SimpleVecWrapper(){
                has_ovnership = true;
                vec = new std::vector<T>;
            }

            ~SimpleVecWrapper(){
                if (has_ovnership){
                    delete vec;
                }
            }

            SimpleVecWrapper(std::vector<T> & vec_in){
                vec = &vec_in;
            }

            inline T& operator[](std::size_t idx){
                return (*vec)[idx];
            }

            inline const T& operator[](std::size_t idx) const{
                return (*vec)[idx];
            }

            inline std::size_t size() const {
                return vec->size();
            }
        };

        static_assert(is_signal_base_first<SimpleVecWrapper<int>>);

        export
        template<typename T, typename IdxT, typename DataValExpr, 
            typename DataRefExpr, typename SizeExpr, bool is_writeble_b>
        struct ExpressionWrapper{
            using SampleType = T;
            using IdxType = IdxT; 
            using DataValueExpression = DataValExpr;
            using DataReferenceExpression = DataRefExpr;
            using SizeExpression = SizeExpr;

            constexpr static bool is_signal_base = true;
            constexpr static bool is_writable = is_writeble_b;

            DataValExpr * val_expression;
            DataRefExpr * ref_expression;
            SizeExpr * size_expression;

            ExpressionWrapper(DataValExpr & val_expr, DataRefExpr & ref_expr, SizeExpr & size_expr) {
                val_expression = &val_expr;
                ref_expression = &ref_expr;
                size_expression = &size_expr;
            }

            inline T& operator[](std::size_t idx){
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return ref_expression(idx);
                }
                else{
                    std::unreachable();
                }
            }
            inline T operator[](std::size_t idx) const{
                if constexpr (! std::is_same_v<DataValueExpression, GENERAL::Nil>){
                    return (*val_expression)(idx);
                }
                else{
                    std::unreachable();
                }
            }

            std::size_t size() const {
                if constexpr (! std::is_same_v<SizeExpression, GENERAL::Nil>){
                    return (*size_expression)();
                }
                else{
                    std::unreachable();
                }
            }
        };

        export
        template<typename T, typename IdxT, typename DataValExpr, 
            typename SizeExpr>
        struct ExpressionWrapper<T, IdxT, DataValExpr, 
            GENERAL::Nil, SizeExpr, false>
        {
            using SampleType = T;
            using IdxType = IdxT; 
            using DataValueExpression = DataValExpr;
            using DataReferenceExpression = GENERAL::Nil;
            using SizeExpression = SizeExpr;

            constexpr static bool is_signal_base = true;
            constexpr static bool is_writable = false;

            DataValExpr * val_expression;
            SizeExpr * size_expression;

            ExpressionWrapper(DataValExpr & val_expr, SizeExpr & size_expr) {
                val_expression = &val_expr;
                size_expression = &size_expr;
            }

            inline T operator[](std::size_t idx) const{
                if constexpr (! std::is_same_v<DataValueExpression, GENERAL::Nil>){
                    return (*val_expression)(idx);
                }
                else{
                    std::unreachable();
                }
            }

            std::size_t size() const {
                if constexpr (! std::is_same_v<SizeExpression, GENERAL::Nil>){
                    return (*size_expression)();
                }
                else{
                    std::unreachable();
                }
            }
        };

        export enum class SignalKind {Monotone, Stohastic, Harmonic, Smooth};
        
        export
        template<SignalBase BaseT, bool is_writable_b>
        struct GenericSignal{
            bool has_ovnership = false;
            constexpr static bool is_writable = is_writable_b;
            constexpr static bool is_signal = true;
            using Base = BaseT;
            using IdxType = Base::IdxType;
            using SampleType = Base::SampleType;
            Base * base;

            SignalKind kind = SignalKind::Stohastic;

            GenericSignal(Base & base_o){
                base = &base_o;
            }

            GenericSignal(){
                base = new Base;
                has_ovnership = true;
            }

            ~GenericSignal(){
                if (has_ovnership){
                    delete base;
                }
            }

            inline SampleType & operator[](IdxType idx){
                return (*base)[idx];
            }

            inline SampleType operator[](IdxType idx) const {
                return (*base)[idx];
            }

            inline size_t size() const {
                return base->size();
            }

            void show(PlottingKind kind){
                if (kind == PlottingKind::Interpolate){
                    std::vector<SampleType> plotting_data = {};
                    int i = -size();
                    while(i < static_cast<int>(size() * 2)){
                        ++i;
                        plotting_data.push_back(static_cast<float>(interpolate<double>(static_cast<double>(i))));
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                }
                else if (kind == PlottingKind::Simple){
                    std::vector<SampleType> plotting_data = {};
                    for (auto i = 0; i < base->size(); ++i){
                        auto sample = (*base)[i];
                        plotting_data.push_back(sample);
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                }
            }

            void show(PlottingKind kind, const std::string & filename, const std::string & format) const {
                if (kind == PlottingKind::Simple){
                    std::vector<SampleType> plotting_data = {};
                    for (auto i = 0; i < base->size(); ++i){
                        plotting_data.push_back((*base)[i]);
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                    matplot::save(filename, format);
                }
                else if (kind == PlottingKind::Interpolate){
                    std::vector<SampleType> plotting_data = {};
                    int i = -size();
                    while(i < static_cast<int>(size() * 2)){
                        i++;
                        plotting_data.push_back(interpolate<int>(i));
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                    matplot::save(filename, format);
                }
            }

            void show(PlottingKind kind, const std::string & filename) const {
                if (kind == PlottingKind::Simple){
                    std::vector<SampleType> plotting_data = {};
                    for (auto i = 0; i < base->size(); ++i){
                        plotting_data.push_back((*base)[i]);
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                    matplot::save(filename);
                }
                else if (kind == PlottingKind::Interpolate){
                    std::vector<SampleType> plotting_data = {};
                    int i = -size();
                    while(i < static_cast<int>(size() * 2)){
                        ++i;
                        plotting_data.push_back(interpolate<int>(i));
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                    matplot::save(filename);
                }
            }


            //получение значения в неизвестной точке внутри диапазона определения (те в нашем случае по дробному индексу)
            template<typename Idx>
            SampleType interpolate(Idx idx) const {
                if (kind == SignalKind::Monotone) {
                    if (idx>=0 && idx < static_cast<Idx>(size() - 2)){
                        if constexpr (NP_DSP::CONFIG::debug){
                            std::string mark = "1b";
                            IC(mark);
                        }
                        return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(idx), (*base)[static_cast<int>(idx)]}, {static_cast<double>(idx + 1), (*base)[idx + 1]}, static_cast<double>(idx));
                    }
                    else if (idx<0){
                        Idx idx_new = idx + ((0 - static_cast<int>(idx)) % static_cast<int>(size())) * 2;
                        if constexpr (NP_DSP::CONFIG::debug){
                            std::string mark = "2b";
                            IC(mark);
                            IC(idx_new);
                        }
                        if (idx_new == idx){
                            idx_new++;
                        }
                        SampleType value_new = interpolate(idx_new);
                        return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(idx_new), value_new}, {0.0, (*base)[0]}, static_cast<double>(idx));
                    }
                    else if (idx > size() - 1){
                        Idx idx_new = idx - ((static_cast<int>(idx) - static_cast<int>(size())) % static_cast<int>(size())) * 2 - 1;
                        if (idx_new == idx){
                            idx_new--;
                        }
                        if (idx_new > size() - 1){
                            idx_new = size() - 2;
                        }
                        if constexpr (CONFIG::debug){
                            std::string mark = "3b";
                            IC(mark);
                            IC(idx_new);
                        }
                        SampleType value_new = interpolate(idx_new);
                        return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(size()-1), (*base)[size()-1]}, {static_cast<double>(idx_new), value_new}, static_cast<double>(idx));
                    }
                    else{
                        if constexpr (CONFIG::debug){
                            std::string mark = "4b";
                            IC(mark);
                        }
                        return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(size() - 1), (*base)[size() - 1]}, {static_cast<double>(size() - 2), (*base)[size() - 2]}, static_cast<double>(idx));
                    }
                }
                else if (kind == SignalKind::Harmonic) {
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
                    }
                    else {
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
                        std::complex b =  {std::cos(2.0 * std::numbers::pi * i * idx / size()),
                            std::sin(2.0 * std::numbers::pi * i * idx / size())};
                        accum += spectr[i] * b;
                    }
                    return accum.real();
                }
                else if (kind == SignalKind::Stohastic) {
                    if (idx>=0 && idx < static_cast<Idx>(size() - 2)){
                        if constexpr (NP_DSP::CONFIG::debug){
                            std::string mark = "1b";
                            IC(mark);
                        }
                        return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(idx), (*base)[static_cast<int>(idx)]}, {static_cast<double>(idx + 1), (*base)[idx + 1]}, static_cast<double>(idx));
                    }
                    else if (idx<0){
                        Idx idx_new = 0.0 - idx;

                        if (idx_new == idx){
                            ++idx_new;
                        }
                        SampleType value_new = interpolate(idx_new);
                        return value_new;
                    }
                    else if (idx > size() - 1){
                        Idx idx_new = 2.0 * size() - 2.0 - idx;

                        SampleType value_new = interpolate(idx_new);
                        return value_new;
                    }
                    else{
                        return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(size() - 1), (*base)[size() - 1]}, {static_cast<double>(size() - 2), (*base)[size() - 2]}, static_cast<double>(idx));
                    }
                }
                else if (kind == SignalKind::Smooth) {
                    //todo Teilors Series
                    std::unreachable();
                }
            }

            template<typename Idx>
            IdxType findMonotone(SampleType value, std::optional<Idx> idx1, std::optional<Idx> idx2) const {
                if (!idx1){
                    if (*idx2!=0){
                        idx1 = {0};
                    }
                    else{
                        idx1 = {-size()};
                    }
                }
                if (!idx2){
                    if (*idx2!=0){
                        idx1 = {2 * size()};
                    }
                    else{
                        idx2 = {size()};
                    }
                }
                auto idx_lambda = [&](int idx){
                    return interpolate(idx);
                };
                std::pair<int, int> idxes = ONE_D::UTILITY_MATH::interpoationSearch(*idx1, *idx2, value, idx_lambda);
                if constexpr (CONFIG::debug){
                    std::string mark = "find monotone";
                    IC(mark, idxes.first, idxes.second, value);
                }
                return UTILITY_MATH::backLinearInterpolate<Idx, SampleType>({*idx1, (*base)[*idx1]}, {*idx2, (*base)[*idx2]}, value);
            }
        };

        export
        template<SignalBase BaseT>
        struct GenericSignal<BaseT, false>{
            bool has_ovnership = false;
            constexpr static bool is_writable = false;
            constexpr static bool is_signal = true;
            using Base = BaseT;
            using IdxType = Base::IdxType;
            using SampleType = Base::SampleType;
            Base * base;

            GenericSignal(Base & base_o){
                base = &base_o;
            }

            GenericSignal(){
                base = new Base;
                has_ovnership = true;
            }

            ~GenericSignal(){
                if (has_ovnership){
                    delete base;
                }
            }

            inline SampleType operator[](IdxType idx) const {
                return (*base)[idx];
            }

            inline size_t size() const {
                return base->size();
            }

            void show(PlottingKind kind){
                if (kind == PlottingKind::Interpolate){
                    std::vector<SampleType> plotting_data = {};
                    int i = -size();
                    while(i < static_cast<int>(size() * 2)){
                        i++;
                        plotting_data.push_back(static_cast<float>(interpolate<double>(static_cast<double>(i))));
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                }
                else if (kind == PlottingKind::Simple){
                    std::vector<SampleType> plotting_data = {};
                    for (auto i = 0; i < base->size(); i++){
                        auto sample = (*base)[i];
                        plotting_data.push_back(sample);
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                }
            }

            void show(PlottingKind kind, const std::string & filename, const std::string & format) const {
                if (kind == PlottingKind::Simple){
                    std::vector<SampleType> plotting_data = {};
                    for (auto i = 0; i < base->size(); i++){
                        plotting_data.push_back((*base)[i]);
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                    matplot::save(filename, format);
                }
                else if (kind == PlottingKind::Interpolate){
                    std::vector<SampleType> plotting_data = {};
                    int i = -size();
                    while(i < static_cast<int>(size() * 2)){
                        i++;
                        plotting_data.push_back(interpolate<int>(i));
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                    matplot::save(filename, format);
                }
            }

            void show(PlottingKind kind, const std::string & filename) const {
                if (kind == PlottingKind::Simple){
                    std::vector<SampleType> plotting_data = {};
                    for (auto i = 0; i < base->size(); i++){
                        plotting_data.push_back((*base)[i]);
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                    matplot::save(filename);
                }
                else if (kind == PlottingKind::Interpolate){
                    std::vector<SampleType> plotting_data = {};
                    int i = -size();
                    while(i < static_cast<int>(size() * 2)){
                        i++;
                        plotting_data.push_back(interpolate<int>(i));
                    }
                    matplot::plot(plotting_data);
                    matplot::show();
                    matplot::save(filename);
                }
            }


            //получение значения в неизвестной точке внутри диапазона определения (те в нашем случае по дробному индексу)
            template<typename Idx>
            SampleType interpolate(Idx idx) const {
                if (idx>=0 && idx < static_cast<Idx>(size() - 2)){
                    return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(idx), (*base)[static_cast<int>(idx)]}, {static_cast<double>(idx + 1), (*base)[idx + 1]}, static_cast<double>(idx));
                }
                else if (idx<0){
                    Idx idx_new = idx + ((0 - static_cast<int>(idx)) % static_cast<int>(size())) * 2;
                    if (idx_new == idx){
                        idx_new++;
                    }
                    SampleType value_new = interpolate(idx_new);
                    return ONE_D::UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(idx_new), value_new}, {0.0, (*base)[0]}, static_cast<double>(idx));
                }
                else if (idx > size() - 1){
                    Idx idx_new = idx - ((static_cast<int>(idx) - static_cast<int>(size())) % static_cast<int>(size())) * 2 - 1;
                    if (idx_new == idx){
                        idx_new--;
                    }
                    SampleType value_new = interpolate(idx_new);
                    return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(size()-1), (*base)[size()-1]}, {static_cast<double>(idx_new), value_new}, static_cast<double>(idx));
                }
                else{
                    return UTILITY_MATH::linearInterpolate<double, SampleType>({static_cast<double>(size() - 1), (*base)[size() - 1]}, {static_cast<double>(size() - 2), (*base)[size() - 2]}, static_cast<double>(idx));
                }
            }

            template<typename Idx>
            IdxType findMonotone(SampleType value, std::optional<Idx> idx1, std::optional<Idx> idx2) const {
                if (!idx1){
                    if (*idx2!=0){
                        idx1 = {0};
                    }
                    else{
                        idx1 = {-size()};
                    }
                }
                if (!idx2){
                    if (*idx2!=0){
                        idx1 = {2 * size()};
                    }
                    else{
                        idx2 = {size()};
                    }
                }
                auto idx_lambda = [&](int idx){
                    return interpolate(idx);
                };
                std::pair<int, int> idxes = UTILITY_MATH::interpoationSearch(*idx1, *idx2, value, idx_lambda);
                if constexpr (CONFIG::debug){
                    std::string mark = "find monotone";
                    IC(mark, idxes.first, idxes.second, value);
                }
                return UTILITY_MATH::backLinearInterpolate<Idx, SampleType>({*idx1, (*base)[*idx1]}, {*idx2, (*base)[*idx2]}, value);
            }
        };

        static_assert(ONE_D::is_signal<GenericSignal<SimpleVecWrapper<int>, true>>);
    }

    namespace GENERAL{
        export
        template<typename T, typename IdxT, typename DataValExpr, typename DataRefExpr, typename SizeExpr, bool is_writeble_b>
        struct ExpressionWrapper{
            using SampleType = T;
            using IdxType = IdxT;
            using DataValueExpression = DataValExpr;
            using DataReferenceExpression = DataRefExpr;
            using SizeExpression = SizeExpr;
            constexpr static std::size_t dims_count = std::tuple_size_v<IdxType>;

            constexpr static bool is_signal_base = true;
            constexpr static bool is_writable = is_writeble_b;

            DataValExpr val_expression;
            DataRefExpr ref_expression;
            SizeExpr dim_size_expression;

            ExpressionWrapper(DataValExpr val_expr, DataRefExpr ref_expr, SizeExpr size_expr) {
                val_expression = val_expr;
                ref_expression = ref_expr;
                dim_size_expression = size_expr;
            }

            inline T& operator[](std::size_t idx){
                if constexpr (! std::is_same_v<DataReferenceExpression, Nil>){
                    return ref_expression(idx);
                }
                else{
                    std::unreachable();
                }
            }
            inline const T& operator[](std::size_t idx) const{
                if constexpr (! std::is_same_v<DataReferenceExpression, Nil>){
                    return val_expression(idx);
                }
                else{
                    std::unreachable();
                }
            }

            std::size_t dimSize(std::size_t dim_number){
                if constexpr (! std::is_same_v<DataReferenceExpression, Nil>){
                    return dim_size_expression(dim_number);
                }
                else{
                    std::unreachable();
                }
            }
        };
    }
}