#pragma once

#include "Global.h"
#include "expressions_lambda.h"

namespace EH
{
    namespace Matrix
    {
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_same_size< TA , TB >::value >::type
                 >
        _ehm_inline
        auto operator == ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )->bool
                    {
                        return a == b;
                    }
                );
        }
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_same_size< TA , TB >::value >::type
                 >
        _ehm_inline
        auto operator != ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )->bool
                    {
                        return a != b;
                    }
                );
        }
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_same_size< TA , TB >::value >::type
                 >
        _ehm_inline
        auto operator >= ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )->bool
                    {
                        return a >= b;
                    }
                );
        }
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_same_size< TA , TB >::value >::type
                 >
        _ehm_inline
        auto operator <= ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )->bool
                    {
                        return a <= b;
                    }
                );
        }
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_same_size< TA , TB >::value >::type
                 >
        _ehm_inline
        auto operator < ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )->bool
                    {
                        return a < b;
                    }
                );
        }
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< is_same_size< TA , TB >::value >::type
                 >
        _ehm_inline
        auto operator > ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )->bool
                    {
                        return a > b;
                    }
                );
        }






        template < typename TA >
        _ehm_inline
        auto operator == ( const Expression< TA >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return a == scalar;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator != ( const Expression< TA >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return a != scalar;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator >= ( const Expression< TA >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return a >= scalar;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator <= ( const Expression< TA >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return a <= scalar;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator > ( const Expression< TA >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return a > scalar;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator < ( const Expression< TA >& exp , typename expression_traits< TA >::result_type scalar )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return a < scalar;
                    }
                );
        }





        template < typename TA >
        _ehm_inline
        auto operator == ( typename expression_traits< TA >::result_type scalar , const Expression< TA >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar == a;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator != ( typename expression_traits< TA >::result_type scalar , const Expression< TA >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar != a;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator <= ( typename expression_traits< TA >::result_type scalar , const Expression< TA >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar <= a;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator >= ( typename expression_traits< TA >::result_type scalar , const Expression< TA >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar >= a;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator < ( typename expression_traits< TA >::result_type scalar , const Expression< TA >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar < a;
                    }
                );
        }
        template < typename TA >
        _ehm_inline
        auto operator > ( typename expression_traits< TA >::result_type scalar , const Expression< TA >& exp )
        {
            return Expressions::make_unary( exp ,
                    [ scalar ]( auto a )->bool
                    {
                        return scalar > a;
                    }
                );
        }

        template < typename TA , typename TB ,
                   typename = typename std::enable_if< 
                       is_same_size< TA , TB >::value &&
                       std::is_same< bool , typename expression_traits< TA >::result_type >::value &&
                       std::is_same< bool , typename expression_traits< TB >::result_type >::value
                   >::type >
        _ehm_inline
        auto operator && ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( bool a , bool b )->bool
                    {
                        return a && b;
                    }
                );
        }
        template < typename TA , typename TB ,
                   typename = typename std::enable_if< 
                       is_same_size< TA , TB >::value &&
                       std::is_same< bool , typename expression_traits< TA >::result_type >::value &&
                       std::is_same< bool , typename expression_traits< TB >::result_type >::value
                   >::type >
        _ehm_inline
        auto operator || ( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( bool a , bool b )->bool
                    {
                        return a || b;
                    }
                );
        }
    };
};
