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
        template < typename TA >
        _ehm_inline
        auto operator * ( const Expression< TA >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto x )
                    {
                        return x * scalar;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto  operator * ( typename expression_traits< TA >::result_type scalar , const Expression< TA >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto x )
                    {
                        return x * scalar;
                    }
                );
        }

        template < typename TA >
        auto operator / ( const Expression< TA >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary(  exp ,
                    [ scalar ]( auto x )
                    {
                        return x / scalar;
                    }
                );
        }
        template < typename TA >
        auto operator / ( typename expression_traits< TA >::result_type scalar , const Expression< TA >& exp )
        {
            return Expressions::make_unary(  exp ,
                    [ scalar ]( auto x )
                    {
                        return scalar / x;
                    }
                );
        }

        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_same_size< TA , TB >::value >::type
                 >
        _ehm_inline
        auto operator + ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                        []( auto a , auto b )
                        {
                            return a + b;
                        }
                    );
        }
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_same_size< TA , TB >::value >::type
                 >
        _ehm_inline
        auto operator - ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                        []( auto a , auto b )
                        {
                            return a - b;
                        }
                    );
        }

        template < typename TA , typename TB , typename = void , typename = void ,
                   typename = typename std::enable_if< expression_traits< TA >::cols == expression_traits< TB >::rows &&
                                                       ( is_row_vector< TA >::value == false || is_column_vector< TB >::value == false )
                                                     >::type
                 >
        _ehm_inline
        auto operator * ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::MatMatMult< TA , TB >( e1 , e2 );
        }
        template < typename TA , typename TB , typename = void , typename = void , typename = void ,
                   typename = typename std::enable_if< expression_traits< TA >::cols == expression_traits< TB >::rows &&
                                                       is_row_vector< TA >::value && is_column_vector< TB >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator * ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            constexpr const IndexType N = expression_traits< TA >::cols;
            typename std::common_type< typename expression_traits< TA >::result_type ,
                                       typename expression_traits< TB >::result_type >::type sum( 0 );
            for( IndexType i=0; i<N; ++i )
            {
                sum += GetBy( e1 , i , 0 ) * GetBy( e2 , 0 , i );
            }
            return sum;
        }



        template < typename TA >
        _ehm_inline
        auto operator + ( const expression_vector_size_type< TA , 0 >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto x )
                    {
                        return x + scalar;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator - ( const expression_vector_size_type< TA , 0 >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto x )
                    {
                        return x - scalar;
                    }
                );
        }

        template < typename TA >
        _ehm_inline
        auto operator + ( typename expression_traits< TA >::result_type scalar , const expression_vector_size_type< TA , 0 >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto x )
                    {
                        return x + scalar;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator - ( typename expression_traits< TA >::result_type scalar , const expression_vector_size_type< TA , 0 >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto x )
                    {
                        return scalar - x;
                    }
                );
        }


        template < typename TA , typename TB , typename = void ,
                   typename = typename std::enable_if< is_vector< TA >::value && is_vector< TB >::value &&
                                                       is_same_size< TA , TB >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator * ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )
                    {
                        return a * b;
                    }
                );
        }
        template < typename TA , typename TB , typename = void ,
                   typename = typename std::enable_if< is_vector< TA >::value && is_vector< TB >::value &&
                                                       is_same_size< TA , TB >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator / ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )
                    {
                        return a / b;
                    }
                );
        }

        template < typename TA , typename TB , typename = void ,
                   typename = typename std::enable_if< is_square< TA >::value && is_vector< TB >::value &&
                                                       expression_traits< TA >::rows == vector_size< TB >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator + ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , Expressions::Vector2Square< TB >( e2 ) ,
                    []( auto a , auto b )
                    {
                        return a + b;
                    }
                );
        }
        template < typename TA , typename TB , typename = void , typename = void ,
                   typename = typename std::enable_if< is_square< TB >::value && is_vector< TA >::value &&
                                                       expression_traits< TB >::rows == vector_size< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator + ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e2 , Expressions::Vector2Square< TA >( e1 ) ,
                    []( auto a , auto b )
                    {
                        return a + b;
                    }
                );
        }
        template < typename TA , typename TB , typename = void ,
                   typename = typename std::enable_if< is_square< TA >::value && is_vector< TB >::value &&
                                                       expression_traits< TA >::rows == vector_size< TB >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator - ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , Expressions::Vector2Square< TB >( e2 ) ,
                    []( auto a , auto b )
                    {
                        return a - b;
                    }
                );
        }
        template < typename TA , typename TB , typename = void , typename = void ,
                   typename = typename std::enable_if< is_square< TB >::value && is_vector< TA >::value &&
                                                       expression_traits< TB >::rows == vector_size< TA >::value
                                                     >::type
                 >
        _ehm_inline
        auto operator - ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e2 , Expressions::Vector2Square< TA >( e1 ) ,
                    []( auto a , auto b )
                    {
                        return b - a;
                    }
                );
        }
    };
};
