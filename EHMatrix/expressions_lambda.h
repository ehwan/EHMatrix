#pragma once

#include "Global.h"
#include "expression_traits.h"

namespace EH
{
    namespace Matrix
    {
        namespace Expressions
        {
            template < typename TA , typename FUNC , int OPADD >
            struct Unary :
                Expression< Unary< TA , FUNC , OPADD > >
            {
                using typename Expression< Unary< TA , FUNC , OPADD > >::result_type;

                auto_reference< TA > a;
                FUNC func;

                template < typename FFUNC >
                Unary( auto_reference< TA > _a , FFUNC&& _func ) :
                    a( _a ) ,
                    func( _func )
                {
                }

                result_type
                constexpr _ehm_inline
                Get( IndexType i ) const
                {
                    return func( GetBy( a , i ) );
                }
                result_type
                constexpr _ehm_inline
                Get( IndexType x , IndexType y ) const
                {
                    return func( GetBy( a , x , y ) );
                }
            };
            template < int OPADD = 1 , typename FUNC , typename TA >
            Unary< TA ,
                   FUNC ,
                   OPADD >
            _ehm_inline
            make_unary( TA&& exp , FUNC&& func )
            {
                return Unary< TA ,
                              FUNC ,
                              OPADD >
                                ( std::forward< TA >( exp ) , std::forward< FUNC >( func ) );
            }
            template < typename TA , typename TB , typename FUNC , int OPADD >
            struct Binary :
                Expression< Binary< TA , TB , FUNC , OPADD > >
            {
                using typename Expression< Binary< TA , TB , FUNC , OPADD > >::result_type;

                auto_reference< TA > a;
                auto_reference< TB > b;
                FUNC func;

                template < typename FFUNC >
                Binary( auto_reference< TA > _a , auto_reference< TB > _b , FFUNC&& _func ) :
                    a( _a ) , b( _b ) , func( _func )
                {
                }

                result_type
                _ehm_inline
                Get( IndexType i ) const
                {
                    return func( GetBy( a , i ) , GetBy( b , i ) );
                }
                result_type
                _ehm_inline
                Get( IndexType x , IndexType y ) const
                {
                    return func( GetBy( a , x , y ) , GetBy( b , x , y ) );
                }
            };
            template < int OPADD = 1 , typename FUNC , typename TA , typename TB >
            Binary< TA ,
                    TB ,
                    FUNC ,
                    OPADD >
            _ehm_inline
            make_binary( TA&& exp1 , TB&& exp2 , FUNC&& func )
            {
                return Binary< TA ,
                               TB ,
                               FUNC ,
                               OPADD >
                                   ( std::forward< TA >( exp1 ) , std::forward< TB >( exp2 ) , std::forward< FUNC >( func ) );
            }
        };  // namespace Expressions

        template < typename TA , typename FUNC , int OPADD >
        struct expression_traits< Expressions::Unary< TA , FUNC , OPADD > > : expression_traits< TA >
        {
            using result_type = typename std::result_of< FUNC( typename expression_traits< TA >::result_type ) >::type;
            _ehm_const int operations = expression_traits< TA >::operations + OPADD;
            _ehm_const bool catch_reference = false;
        };
        template < typename TA , typename TB , typename FUNC , int OPADD >
        struct expression_traits< Expressions::Binary< TA , TB , FUNC , OPADD > > : expression_traits< TA , TB >
        {
            using result_type = typename std::result_of< FUNC( typename expression_traits< TA >::result_type ,
                                                               typename expression_traits< TB >::result_type ) >::type;
            _ehm_const int operations = expression_traits< TA >::operations + expression_traits< TB >::operations + OPADD;

            _ehm_const bool catch_reference = false;
        };
    };  // namespace Matrix
};  // namespace EH
