#pragma once

#include "Global.h"
#include "expression_traits.h"
#include "expressions.h"
#include "expressions_lambda.h"

#include <functional>

namespace EH
{
    namespace Matrix
    {
        namespace functional

        {
            template < typename T2 >
            struct Multiply_Unary
            {
                T2 s;

                Multiply_Unary( T2 _s ) :
                    s( _s )
                {
                }

                template < typename T1 >
                typename std::common_type< T1 , T2 >::type
                _ehm_inline
                operator()( T1 x ) const
                {
                    return x * s;
                }
            };
            template < typename TA , typename TB >
            struct Multiply_Binary
            {
                typename std::common_type< TA , TB >::type
                _ehm_inline
                operator() ( TA a , TB b ) const
                {
                    return a * b;
                }
            };
        };
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator * ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto x )
                    {
                        return x * scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator * ( typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto x )
                    {
                        return x * scalar;
                    }
                );
        }

        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator / ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto x )
                    {
                        return x / scalar;
                    }
                );
        }

        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_expression< TA >::value && is_expression< TB >::value &&
                                                       expression_traits< TA >::cols == expression_traits< TB >::cols &&
                                                       expression_traits< TA >::rows == expression_traits< TB >::rows
                                                     >::type
                 >
        _ehm_inline
        auto operator + ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                        []( auto a , auto b )
                        {
                            return a + b;
                        }
                    );
        }
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_expression< TA >::value && is_expression< TB >::value &&
                                                       expression_traits< TA >::cols == expression_traits< TB >::cols &&
                                                       expression_traits< TA >::rows == expression_traits< TB >::rows
                                                     >::type
                 >
        _ehm_inline
        auto operator - ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                        []( auto a , auto b )
                        {
                            return a - b;
                        }
                    );
        }

        template < typename TA , typename TB , typename = void , typename = void ,
                   typename = typename std::enable_if< is_expression< TA >::value && is_expression< TB >::value &&
                                                       expression_traits< TA >::cols == expression_traits< TB >::rows &&
                                                       ( is_row_vector< TA >::value == false || is_column_vector< TB >::value == false )
                                                     >::type
                 >
        _ehm_inline
        auto operator * ( TA&& e1 , TB&& e2 )
        {
            return Expressions::MatMatMult< TA , TB >( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) );
        }
        template < typename TA , typename TB , typename = void , typename = void , typename = void ,
                   typename = typename std::enable_if< is_expression< TA >::value && is_expression< TB >::value &&
                                                       expression_traits< TA >::cols == expression_traits< TB >::rows &&
                                                       is_row_vector< TA >::value && is_column_vector< TB >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator * ( TA&& e1 , TB&& e2 )
        {
            constexpr const IndexType N = expression_traits< TA >::cols;
            typename std::common_type< typename expression_traits< TA >::result_type ,
                                       typename expression_traits< TB >::result_type >::type sum( 0 );
            for( IndexType i=0; i<N; ++i )
            {
                sum += GetBy( std::forward< TA >( e1 ) , i , 0 ) * GetBy( std::forward< TB >( e2 ) , 0 , i );
            }
            return sum;
        }



        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value &&
                                                       is_vector< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator + ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto x )
                    {
                        return x + scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value &&
                                                       is_vector< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator + ( typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto x )
                    {
                        return x + scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value &&
                                                       is_vector< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator - ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto x )
                    {
                        return x - scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value &&
                                                       is_vector< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator - ( typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto x )
                    {
                        return scalar - x;
                    }
                );
        }

        template < typename TA , typename TB , typename = void ,
                   typename = typename std::enable_if< is_expression< TA >::value && is_expression< TB >::value &&
                                                       is_vector< TA >::value && is_vector< TB >::value &&

                                                       expression_traits< TA >::cols == expression_traits< TB >::cols &&
                                                       expression_traits< TA >::rows == expression_traits< TB >::rows
                                                     >::type
                 >
        _ehm_inline
        auto operator * ( TA&& e1 , TB&& e2 )
        {
            //return Expressions::Multiply< TA , TB >( e1 , e2 );
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                    []( auto a , auto b )
                    {
                        return a * b;
                    }
                );
        }
        template < typename TA , typename TB , typename = void ,
                   typename = typename std::enable_if< is_expression< TA >::value && is_expression< TB >::value &&
                                                       is_vector< TA >::value && is_vector< TB >::value &&

                                                       expression_traits< TA >::cols == expression_traits< TB >::cols &&
                                                       expression_traits< TA >::rows == expression_traits< TB >::rows
                                                     >::type
                 >
        _ehm_inline
        auto operator / ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                    []( auto a , auto b )
                    {
                        return a / b;
                    }
                );
        }
    };
};
