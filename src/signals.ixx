module;

import <vector>;
import <tuple>;
import <cassert>;
import <cstddef>;
import <utility>;
import <type_traits>;

import npdsp_concepts;

export module signals;


namespace NP_DSP
{
    namespace GENERAL{
        struct Nil{};
    }

    namespace ONE_D{
        export
        template<typename T>
        struct SimpleVecWrapper{
            using DataType = T;
            using IdxType = size_t;
            constexpr static bool is_signal = true;
            constexpr static bool is_writable = true;

            std::vector<T> & vec;
            inline T& getRefByIdx(size_t idx){
                return vec[idx];
            }

            inline T getByIdx(size_t idx){
                return vec[idx];
            }

            size_t getSize(){
                return vec.size();
            }
        };

        static_assert(is_signal<SimpleVecWrapper<int>>);

        export
        template<typename T, typename DataValExpr, typename DataRefExpr, typename SizeExpr, bool is_writeble_b>
        struct ExpressionWrapper{
            using DataType = T;
            using IdxType = size_t;
            using DataValueExpression = DataValExpr;
            using DataReferenceExpression = DataRefExpr;
            using SizeExpression = SizeExpr;

            constexpr static bool is_signal = true;
            constexpr static bool is_writable = is_writeble_b;

            DataValExpr val_expression;
            DataRefExpr ref_expression;
            SizeExpr size_expression;

            ExpressionWrapper(DataValExpr val_expr, DataRefExpr ref_expr, SizeExpr size_expr) {
                val_expression = val_expr;
                ref_expression = ref_expr;
                size_expression = size_expr;
            }

            inline T& getRefByIdx(size_t idx){
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return ref_expression(idx);
                }
                else{
                    std::unreachable();
                }
            }

            inline T getByIdx(size_t idx){
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return val_expression(idx);
                }
                else{
                    std::unreachable();
                }
            }

            size_t getSize(){
                if constexpr (! std::is_same_v<DataReferenceExpression, GENERAL::Nil>){
                    return size_expression();
                }
                else{
                    std::unreachable();
                }
            }
        };
    }
}