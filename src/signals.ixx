module;

import <vector>;
import <tuple>;
import <cassert>;
import <cstddef>;
import <utility>;
import <type_traits>;
import <optional>;

import npdsp_concepts;

export module signals;


namespace NP_DSP
{
    namespace ONE_D{
        export
        template<typename T>
        struct SimpleVecWrapper{
            using SampleType = T;
            using IdxType = std::size_t;
            constexpr static bool is_signal_base = true;
            constexpr static bool is_writable = true;

            std::vector<T> & vec;
            inline T& operator[](std::size_t idx){
                return vec[idx];
            }

            inline const T& operator[](std::size_t idx) const{
                return vec[idx];
            }

            inline std::size_t size(){
                return vec.size();
            }
        };

        static_assert(is_signal_base<SimpleVecWrapper<int>>);

        export
        template<typename T, typename IdxT, typename DataValExpr, typename DataRefExpr, typename SizeExpr, bool is_writeble_b>
        struct ExpressionWrapper{
            using SampleType = T;
            using IdxType = IdxT;
            using DataValueExpression = DataValExpr;
            using DataReferenceExpression = DataRefExpr;
            using SizeExpression = SizeExpr;

            constexpr static bool is_signal_base = true;
            constexpr static bool is_writable = is_writeble_b;

            DataValExpr val_expression;
            DataRefExpr ref_expression;
            SizeExpr size_expression;

            ExpressionWrapper(DataValExpr val_expr, DataRefExpr ref_expr, SizeExpr size_expr) {
                val_expression = val_expr;
                ref_expression = ref_expr;
                size_expression = size_expr;
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
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return val_expression(idx);
                }
                else{
                    std::unreachable();
                }
            }

            std::size_t size(){
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return size_expression();
                }
                else{
                    std::unreachable();
                }
            }
        };

        export
        template<SignalBase BaseT>
        struct GenericSignal{
            constexpr static bool is_signal = true;
            using Base = BaseT;
            using IdxType = Base::IdxType;
            using SampleType = Base::SapleType;
            Base base;

            explicit GenericSignal(Base & base_o){
                base = base_o;
            }
            inline SampleType & operator[](IdxType idx){
                return base[idx];
            }

            inline SampleType operator[](IdxType idx) const {
                return base[idx];
            }

            inline size_t size(){
                return base.size();
            }

            template<typename Idx>
            SampleType interpolate(Idx idx){
                //todo
            }
            template<typename Idx>
            IdxType findInterpolate(SampleType value, std::optional<Idx>, std::optional<Idx>){
                //todo
            }
            template<typename Idx>
            IdxType findIncr(SampleType value, std::optional<Idx>, std::optional<Idx>){
                //todo
            }
            template<typename Idx>
            IdxType findDecr(SampleType value, std::optional<Idx>, std::optional<Idx>){
                //todo
            }
        };

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
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return ref_expression(idx);
                }
                else{
                    std::unreachable();
                }
            }
            inline const T& operator[](std::size_t idx) const{
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return val_expression(idx);
                }
                else{
                    std::unreachable();
                }
            }

            std::size_t dimSize(std::size_t dim_number){
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return dim_size_expression(dim_number);
                }
                else{
                    std::unreachable();
                }
            }
        };
    }
}

//todo padding