#pragma once

#include "Global.h"
#include "expressions_lambda.h"

namespace EH
{
    namespace Matrix
    {
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_expression< TA >::value && is_expression< TB >::value &&
                                              expression_traits< TA >::cols == expression_traits< TB >::cols &&
                                              expression_traits< TA >::rows == expression_traits< TB >::rows
                                            >::type
                 >
        _ehm_inline
        auto operator == ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                    []( auto a , auto b )->bool
                    {
                        return a == b;
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
        auto operator != ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                    []( auto a , auto b )->bool
                    {
                        return a != b;
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
        auto operator >= ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                    []( auto a , auto b )->bool
                    {
                        return a >= b;
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
        auto operator <= ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                    []( auto a , auto b )->bool
                    {
                        return a <= b;
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
        auto operator < ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                    []( auto a , auto b )->bool
                    {
                        return a < b;
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
        auto operator > ( TA&& e1 , TB&& e2 )
        {
            return Expressions::make_binary( std::forward< TA >( e1 ) , std::forward< TB >( e2 ) ,
                    []( auto a , auto b )->bool
                    {
                        return a > b;
                    }
                );
        }






        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator == ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return a == scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator != ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return a != scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator >= ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return a >= scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator <= ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return a <= scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator > ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return a > scalar;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator < ( TA&& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return a < scalar;
                    }
                );
        }





        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator == (typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar == a;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator != (typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar != a;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator <= (typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar <= a;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator >= (typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar >= a;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator < (typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar < a;
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if< is_expression< TA >::value >::type
                 >
        _ehm_inline
        auto operator > (typename expression_traits< TA >::result_type scalar , TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar > a;
                    }
                );
        }
    };
};
