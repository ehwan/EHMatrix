#pragma once

#include "Global.h"
#include "expression_interface.h"
#include "expression_traits.h"
#include "aliased_container.h"

#include <type_traits>
#include <initializer_list>

namespace EH
{
    namespace Matrix
    {
        template < typename T , IndexType M , IndexType N >
        struct Matrix :
            aliased_container< T , M , N > ,
            //Expression< Matrix< T , M , N > >
            WritableExpression< Matrix< T , M , N > >
        {
            using parent = WritableExpression< Matrix< T , M , N > >;
            using typename parent::result_type;

            using parent::operator=;

            constexpr _ehm_inline T Get( IndexType i ) const
            {
                return aliased_container< T , M , N >::s[ i ];
            }
            constexpr _ehm_inline T Get( IndexType x , IndexType y ) const
            {
                return aliased_container< T , M , N >::s[ y + x*M ];
            }
            constexpr _ehm_inline T& Ref( IndexType i )
            {
                return aliased_container< T , M , N >::s[ i ];
            }
            constexpr _ehm_inline T& Ref( IndexType x , IndexType y )
            {
                return aliased_container< T , M , N >::s[ y + x*M ];
            }

            constexpr Matrix(){}

            template < typename ... Ts ,
                       typename = typename std::enable_if<
                           EH::static_sequence< std::size_t , expression_traits< Ts >::rows*expression_traits< Ts >::cols ... >::sum::value == M*N
                        >::type >
            constexpr Matrix( Ts&& ... args )
            {
                parent::FillAggressive( std::forward< Ts >( args )... );
            }

            template < typename ... Ts , typename = void ,
                       typename = typename std::enable_if<
                           M == N &&
                           EH::static_sequence< std::size_t , expression_traits< Ts >::rows*expression_traits< Ts >::cols ... >::sum::value == M
                        >::type >
            constexpr Matrix( Ts&& ... args )
            {
                parent::Fill( 0 );
                parent::Diagonal().FillAggressive( std::forward< Ts >( args )... );
            }

            template < IndexType M0 = M , IndexType N0 = N ,
                       typename = typename std::enable_if< N0 == M0 >::type >
            constexpr Matrix( const T a )
            {
                parent::Fill( 0 );
                parent::Diagonal().Fill( a );
            }
            template < IndexType M0 = M , IndexType N0 = N , typename = void ,
                       typename = typename std::enable_if< M0 == 1 || N0 == 1 >::type >
            constexpr Matrix( const T a )
            {
                parent::Fill( a );
            }
        };

        template < typename T , IndexType M , IndexType N >
        struct expression_traits< Matrix< T , M , N > >
        {
            using result_type = T;
            using root_type   = Matrix< T , M , N >;

            _ehm_const IndexType cols            = N;
            _ehm_const IndexType rows            = M;

            _ehm_const bool      is_single_index = true;
            _ehm_const bool      is_restrict     = true;
            _ehm_const bool      catch_reference = true;

            _ehm_const bool      is_expression   = true;

            _ehm_const int       operations      = 0;
        };
    };  // namespace Matrix
};  //namespace EH
