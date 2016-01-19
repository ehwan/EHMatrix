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

                auto_reference< const TA > a;
                FUNC func;

                template < typename FFUNC >
                Unary( auto_reference< const TA > _a , FFUNC&& _func ) :
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
            make_unary( const Expression< TA >& exp , FUNC&& func )
            {
                return Unary< TA , FUNC , OPADD >
                                ( exp , std::forward< FUNC >( func ) );
            }
            template < typename TA , typename TB , typename FUNC , int OPADD >
            struct Binary :
                Expression< Binary< TA , TB , FUNC , OPADD > >
            {
                using typename Expression< Binary< TA , TB , FUNC , OPADD > >::result_type;

                auto_reference< const TA > a;
                auto_reference< const TB > b;
                FUNC func;

                template < typename FFUNC >
                Binary( auto_reference< const TA > _a , auto_reference< const TB > _b , FFUNC&& _func ) :
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
            make_binary( const Expression< TA >& exp1 , const Expression< TB >& exp2 , FUNC&& func )
            {
                return Binary< TA , TB , FUNC , OPADD >
                                   ( exp1 , exp2 , std::forward< FUNC >( func ) );
            }
            template < typename TA , bool SINGLE , bool RESTRICT , int OPADD , IndexType M , IndexType N , typename FUNC >
            struct IndexFilter :
                WritableExpression< IndexFilter< TA , SINGLE , RESTRICT , OPADD , M , N , FUNC > >
            {
                using typename WritableExpression< IndexFilter< TA , SINGLE , RESTRICT , OPADD , M , N , FUNC > >::result_type;

                auto_reference< TA > a;
                FUNC func;

                template < typename FFUNC >
                IndexFilter( auto_reference< TA > _a , FFUNC&& _func ) :
                    a( _a ) , func( _func )
                {
                }

                _ehm_inline
                result_type Get( IndexType i ) const
                {
                    func( i );
                    return GetBy( a , i );
                }
                _ehm_inline
                result_type Get( IndexType x , IndexType y ) const
                {
                    func( x , y );
                    return GetBy( a , x , y );
                }

                _ehm_inline
                result_type& Ref( IndexType i )
                {
                    func( i );
                    return RefBy( a , i );
                }
                _ehm_inline
                result_type Ref( IndexType x , IndexType y )
                {
                    func( x , y );
                    return RefBy( a , x , y );
                }
            };
            template < bool SINGLE , bool RESTRICT , IndexType M , IndexType N , int OPADD=1 , typename FUNC , typename TA >
            auto make_index_filter( Expression< TA >& exp , FUNC&& func )
            {
                return IndexFilter< TA , SINGLE , RESTRICT , OPADD , M , N , FUNC >( exp , std::forward< FUNC >( func ) );
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
        template < typename TA , bool SINGLE , bool RESTRICT , int OPADD , IndexType M , IndexType N , typename FUNC >
        struct expression_traits< Expressions::IndexFilter< TA , SINGLE , RESTRICT , OPADD , M , N , FUNC > > : expression_traits< TA >
        {
            _ehm_const IndexType cols = N;
            _ehm_const IndexType rows = M;
            _ehm_const int operations = expression_traits< TA >::operations + OPADD;
            _ehm_const bool is_single_index = SINGLE;
            _ehm_const bool is_restrict = RESTRICT;
            _ehm_const bool catch_reference = false;
        };
    };  // namespace Matrix
};  // namespace EH
