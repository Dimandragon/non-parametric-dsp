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
import <functional>;

import npdsp_concepts;
import utility_math;
import npdsp_config;

namespace NP_DSP::ONE_D {
    export
    template<typename T>
    struct SimpleVecWrapper : SignalBase<T> {
        bool has_ovnership = false;

        constexpr static bool is_signal_base = true;
        constexpr static bool is_writable = true;

        std::vector<T>* vec;

        SimpleVecWrapper() {
            has_ovnership = true;
            vec = new std::vector<T>;
            this->ref_fn = [&](size_t idx)-> T& {
                return (*vec)[idx];
            };
            this->val_fn = [&](size_t idx)-> T const {
                return (*vec)[idx];
            };
            this->size_fn = [&]()-> size_t const {
                return vec->size();
            };
        }

        ~SimpleVecWrapper() override {
            if (has_ovnership) {
                delete vec;
            }
        }

        SimpleVecWrapper(std::vector<T>& vec_in) {
            vec = &vec_in;
            this->ref_fn = [&](size_t idx)-> T& {
                return (*vec)[idx];
            };
            this->val_fn = [&](size_t idx)-> T const {
                return (*vec)[idx];
            };
            this->size_fn = [&]()-> size_t const {
                return vec->size();
            };
        }

        /*
        T& operator[](std::size_t idx) override {
            return (*vec)[idx];
        }

        T operator[](std::size_t idx) const override {
            return (*vec)[idx];
        }

        size_t size() const override {
            return vec->size();
        }
        */
    };

    //static_assert(is_signal_base_first<SimpleVecWrapper<int>>);

    export
    template<typename T, typename IdxT, typename DataValExpr,
        typename DataRefExpr, typename SizeExpr, bool is_writeble_b>
    struct ExpressionWrapper : SignalBase<T> {
        using SampleType = T;
        using IdxType = IdxT;
        using DataValueExpression = DataValExpr;
        using DataReferenceExpression = DataRefExpr;
        using SizeExpression = SizeExpr;

        constexpr static bool is_signal_base = true;
        constexpr static bool is_writable = is_writeble_b;

        ExpressionWrapper(DataValExpr& val_expr, DataRefExpr& ref_expr, SizeExpr& size_expr) {
            this->val_fn = val_expr;
            this->ref_fn = ref_expr;
            this->size_fn = size_expr;
        }
    };

    export
    template<typename T, typename IdxT, typename DataValExpr,
        typename SizeExpr>
    struct ExpressionWrapper<T, IdxT, DataValExpr, GENERAL::Nil, SizeExpr, false> : SignalBase<T> {
        using SampleType = T;
        using IdxType = IdxT;
        using DataValueExpression = DataValExpr;
        using SizeExpression = SizeExpr;

        constexpr static bool is_signal_base = true;
        constexpr static bool is_writable = false;

        ExpressionWrapper(DataValExpr& val_expr, SizeExpr& size_expr) {
            //IC(&val_expr, &size_expr);
            this->val_fn = val_expr;
            this->size_fn = size_expr;
        }
    };


    export
    template<typename T, bool is_writable_b>
    struct GenericSignal : Signal<T> {
        bool has_ovnership = false;
        constexpr static bool is_writable = is_writable_b;
        constexpr static bool is_signal = true;
        using Base = SignalBase<T>;
        using IdxType = typename Base::IdxType;
        using SampleType = typename Base::SampleType;

        GenericSignal(Base& base_o) {
            this->base = &base_o;
        }


        GenericSignal() {
            this->base = nullptr;
        }

        template<typename BaseT>
        GenericSignal(GENERAL::Tag<BaseT> tag) {
            this->base = new BaseT;
            has_ovnership = true;
        }

        ~GenericSignal() override final {
            if (has_ovnership) {
                delete this->base;
            }
        }

        //получение значения в неизвестной точке внутри диапазона определения (те в нашем случае по дробному индексу)

        SampleType interpolate(double idx, SignalKind kind) const override final {
            using Idx = double;
            if (kind == SignalKind::Monotone) {
                if (idx >= 0 && idx < static_cast<Idx>(this->size() - 2)) {
                    if constexpr (CONFIG::debug) {
                        std::string mark = "1b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(idx), (*this->base)[static_cast<int>(idx)]},
                        {static_cast<double>(idx + 1), (*this->base)[idx + 1]}, static_cast<double>(idx));
                } else if (idx < 0) {
                    Idx idx_new = idx + ((0 - static_cast<int>(idx)) % static_cast<int>(this->size())) * 2;
                    if constexpr (CONFIG::debug) {
                        std::string mark = "2b";
                        IC(mark);
                        IC(idx_new);
                    }
                    if (idx_new == idx) {
                        ++idx_new;
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {idx_new, value_new}, {0.0, (*this->base)[0]}, static_cast<double>(idx));
                } else if (idx > this->size() - 1) {
                    Idx idx_new = idx - ((static_cast<int>(idx) - static_cast<int>(this->size())) % static_cast<int>(
                                             this->size())) * 2 - 1;
                    if (idx_new == idx) {
                        --idx_new;
                    }
                    if (idx_new > this->size() - 1) {
                        idx_new = this->size() - 2;
                    }
                    if constexpr (CONFIG::debug) {
                        std::string mark = "3b";
                        IC(mark);
                        IC(idx_new);
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(this->size() - 1), (*this->base)[this->size() - 1]}, {idx_new, value_new},
                        idx);
                } else {
                    if constexpr (CONFIG::debug) {
                        std::string mark = "4b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(this->size() - 1), (*this->base)[this->size() - 1]}, {
                            static_cast<double>(this->size() - 2), (*this->base)[this->size() - 2]
                        }, idx);
                }
            } else if (kind == SignalKind::Harmonic) {
                if constexpr (GENERAL::is_complex<T>) {
                    std::vector<std::complex<double>> data;
                    std::vector<std::complex<double>> spectr;
                    if (spectr.size() != this->size()) {
                        spectr.clear();
                        data.clear();
                        for (int i = 0; i < this->size(); i++) {
                            spectr.push_back({0.0, 0.0});
                            data.push_back((*this)[i]);
                        }
                    } else {
                        for (int i = 0; i < this->size(); i++) {
                            data.push_back((*this)[i]);
                            spectr[i] = {0.0, 0.0};
                        }
                    }
                    UTILITY_MATH::fftc2c(data, spectr);
                    //IC()
                    //spectr.show(PlottingKind::Simple);


                    std::complex accum = {0.0, 0.0};
                    for (int i = 0; i < this->size(); i++) {
                        std::complex b = {
                            std::cos(2.0 * std::numbers::pi * i * idx / this->size()),
                            std::sin(2.0 * std::numbers::pi * i * idx / this->size())
                        };
                        accum += spectr[i] * b;
                    }
                    return accum;
                } else {
                    std::vector<std::complex<double>> data;
                    std::vector<std::complex<double>> spectr;
                    if (spectr.size() != this->size()) {
                        spectr.clear();
                        data.clear();
                        for (int i = 0; i < this->size(); i++) {
                            spectr.push_back({0.0, 0.0});
                            data.push_back({(*this)[i], 0.0});
                        }
                    } else {
                        for (int i = 0; i < this->size(); i++) {
                            data[i] = {(*this)[i], 0.0};
                            spectr[i] = {0.0, 0.0};
                        }
                    }
                    UTILITY_MATH::fftc2c(data, spectr);
                    //IC()
                    //spectr.show(PlottingKind::Simple);


                    std::complex accum = {0.0, 0.0};
                    for (int i = 0; i < this->size(); i++) {
                        std::complex b = {
                            std::cos(2.0 * std::numbers::pi * i * idx / this->size()),
                            std::sin(2.0 * std::numbers::pi * i * idx / this->size())
                        };
                        accum += spectr[i] * b;
                    }
                    return accum.real();
                }
            } else if (kind == SignalKind::Stohastic) {
                if (idx >= 0 && idx < static_cast<Idx>(this->size() - 2)) {
                    if constexpr (CONFIG::debug) {
                        std::string mark = "1b";
                        IC(mark);
                    }
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {idx, (*this->base)[static_cast<int>(idx)]}, {idx + 1, (*this->base)[idx + 1]}, idx);
                } else if (idx < 0) {
                    Idx idx_new = 0.0 - idx;

                    if (idx_new == idx) {
                        ++idx_new;
                    }
                    SampleType value_new = interpolate(idx_new, kind);
                    return value_new;
                } else if (idx > this->size() - 1) {
                    Idx idx_new = 2.0 * this->size() - 2.0 - idx;

                    SampleType value_new = interpolate(idx_new, kind);
                    return value_new;
                } else {
                    return UTILITY_MATH::linearInterpolate<double, SampleType>(
                        {static_cast<double>(this->size() - 1), (*this->base)[this->size() - 1]}, {
                            static_cast<double>(this->size() - 2), (*this->base)[this->size() - 2]
                        }, idx);
                }
            } else if (kind == SignalKind::Smooth) {
                //todo Teilors Series
                std::unreachable();
            }
        }
    };
}
