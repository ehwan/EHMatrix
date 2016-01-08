#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_head.h"

#include <initializer_list>

namespace EH
{
    namespace Matrix
    {
        template < typename THIS , IndexType N >
        struct Matrix_Interface< THIS , N , N >
        {
            typedef Expression::ret_type< THIS > T;

            _ehm_inline void operator *= ( const T a )
            {
                static_cast< THIS& >( *this ).Multiply( a );
            }
            _ehm_inline void operator /= ( const T a )
            {
                static_cast< THIS& >( *this ).Divide( a );
            }

            template < typename CLS >
            void operator *= ( const Expression::expression_size_type< CLS , N , N >& m )
            {
                static_cast< THIS& >( *this ).FillAggressive( Expression::mat_type< THIS >( static_cast< const THIS& >( *this ) * m ) );
            }

            template < typename CLS >
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , N , N >& exp )
            {
                static_cast< THIS& >( *this ).PlusAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( exp ) );
            }
            template < typename CLS >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , N , N >& exp )
            {
                static_cast< THIS& >( *this ).MinusAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( exp ) );
            }
            _ehm_inline
            void
            operator += ( std::initializer_list< T > lst )
            {
                static_cast< THIS& >( *this ).Plus( lst.begin() , lst.end() );
            }
            _ehm_inline
            void
            operator -= ( std::initializer_list< T > lst )
            {
                static_cast< THIS& >( *this ).Minus( lst.begin() , lst.end() );
            }

            template < typename CLS , typename = void >
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , N , 1 >& v )
            {
                static_cast< THIS& >( *this ).Diagonal().PlusAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( v ) );
            }
            template < typename CLS , typename = void >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , N , 1 >& v )
            {
                static_cast< THIS& >( *this ).Diagonal().MinusAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( v ) );
            }
            _ehm_inline void operator += ( const T a )
            {
                static_cast< THIS& >( *this ).Diagonal().Plus( a );
            }
            _ehm_inline void operator -= ( const T a )
            {
                static_cast< THIS& >( *this ).Diagonal().Minus( a );
            }


            template < typename EXP_CLS >
            _ehm_inline void operator = ( const Expression::expression_size_type< EXP_CLS , N , N >& m )
            {
                static_cast< THIS& >( *this ).FillAggressive( typename Expression::AssignShouldMakeTemp< THIS , EXP_CLS >::type( m ) );
            }
            _ehm_inline
            void
            operator = ( std::initializer_list< T > lst )
            {
                static_cast< THIS& >( *this ).Fill( lst.begin() , lst.end() );
            }
            template < typename EXP_CLS , typename = void >
            _ehm_inline void operator = ( const Expression::expression_size_type< EXP_CLS , N , 1 >& m )
            {
                static_cast< THIS& >( *this ).Fill( 0 );
                static_cast< THIS& >( *this ).Diagonal().FillAggressive( m );
            }
            _ehm_inline
            void
            operator = ( const T a )
            {
                static_cast< THIS& >( *this ).Fill( 0 );
                static_cast< THIS& >( *this ).Diagonal().Fill( a );
            }

            void Log() const
            {
                using Expression::GetBy;
                std::cout << "Matrix Log" << '\n';
                for( IndexType m=0; m<N; ++m )
                {
                    std::cout << "(\t" ;
                    for( IndexType n=0; n<N; ++n )
                    {
                        std::cout << GetBy( static_cast< const THIS& >( *this ) , n , m );
                        if( n!=N-1 ){ std::cout << " , "; }
                    }
                    std::cout << "\t))\n";
                }
            }
        };

        namespace Expression
        {
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator + ( const expression_size_type< TA , M , M >& m ,
                         const expression_size_type< TB , M , 1 >& v )
            {
                return SquarePlus< TA , TB >( m , v );
            }
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator - ( const expression_size_type< TA , M , M >& m ,
                         const expression_size_type< TB , M , 1 >& v )
            {
                return SquarePlus< TA , decltype( v.Negative() ) >( m , v.Negative() );
            }
            template < typename TA , typename TB , typename = void ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator + ( const expression_size_type< TB , M , 1 >& v ,
                         const expression_size_type< TA , M , M >& m )
            {
                return SquarePlus< TA , TB >( m , v );
            }
            template < typename TA , typename TB , typename = void ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator - ( const expression_size_type< TB , M , 1 >& v ,
                         const expression_size_type< TA , M , M >& m )
            {
                return SquarePlus< decltype( m.Negative() ) , TB >( m.Negative() , v );
            }

            template < typename TA ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator + ( const expression_size_type< TA , M , M >& m ,
                         const ret_type< TA > s )
            {
                return SquarePlus< TA , ret_type< TA > >( m , s );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator - ( const expression_size_type< TA , M , M >& m ,
                         const ret_type< TA > s )
            {
                return SquarePlus< TA , ret_type< TA > >( m , -s );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator + ( const ret_type< TA > s ,
                         const expression_size_type< TA , M , M >& m )
            {
                return SquarePlus< TA , ret_type< TA > >( m , s );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator - ( const ret_type< TA > s ,
                         const expression_size_type< TA , M , M >& m )
            {
                return SquarePlus< decltype( m.Negative() ) , ret_type< TA > >( m.Negative() , s );
            }
        };  // namespace Expression

        template < typename TA >
        typename std::enable_if< Expression::Traits< TA >::rows == 2 && Expression::Traits< TA >::cols == 2 ,
                 typename Expression::Traits< TA >::return_type >::type
        Det( const Expression::Expression< TA >& exp )
        {
            using Expression::GetBy;
            return GetBy( exp , 0 , 0 )*GetBy( exp , 1 , 1 ) - GetBy( exp , 0 , 1 )*GetBy( exp , 1 , 0 );
        }

        template < typename TA >
        typename std::enable_if< Expression::Traits< TA >::rows == 2 && Expression::Traits< TA >::cols == 2 ,
                 Expression::mat_type< TA > >::type
        Inverse( const Expression::Expression< TA >& exp )
        {
            using Expression::GetBy;
            typename Expression::ShouldMakeTemp< TA , 2 >::type temp( exp );
            const decltype( Det( exp ) ) invD = 1 / Det( temp );

            return Expression::mat_type< TA >
            {
            {
                invD * GetBy( temp , 1 , 1 ) , -invD * GetBy( temp , 0 , 1 ) , -invD * GetBy( temp , 1 , 0 ) , invD * GetBy( temp , 0 , 0 )
            }
            };
        }

    };  // namespace Matrix

};  //namespace E
