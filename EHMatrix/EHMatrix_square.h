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
                static_cast< THIS& >( *this ).Fill_Safe( static_cast< const THIS& >(*this) * m );
            }

            template < typename CLS >
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , N , N >& exp )
            {
                static_cast< THIS& >( *this ).Plus_Safe( exp );
            }
            template < typename CLS >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , N , N >& exp )
            {
                static_cast< THIS& >( *this ).Minus_Safe( exp );
            }
            template < typename LST_TYPE >
            _ehm_inline
            void
            operator += ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Plus( lst.begin() , lst.end() );
            }
            template < typename LST_TYPE >
            _ehm_inline
            void
            operator -= ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Minus( lst.begin() , lst.end() );
            }

            template < typename CLS , typename = void >
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , N , 1 >& v )
            {
                static_cast< THIS& >( *this ).Diagonal(). Plus_Safe( v );
            }
            template < typename CLS , typename = void >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , N , 1 >& v )
            {
                static_cast< THIS& >( *this ).Diagonal().Minus_Safe( v );
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
                static_cast< THIS& >( *this ).Fill_Safe( m );
            }
            template < typename LST_TYPE >
            _ehm_inline
            void
            operator = ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Fill( lst.begin() , lst.end() );
            }
            template < typename EXP_CLS , typename = void >
            _ehm_inline void operator = ( const Expression::expression_size_type< EXP_CLS , N , 1 >& m )
            {
                static_cast< THIS& >( *this ).Fill( 0 );
                static_cast< THIS& >( *this ).Diagonal() = m;
            }
            template < typename SCALAR_TYPE >
            _ehm_inline
            typename std::enable_if< Expression::is_scalar< SCALAR_TYPE >::value >::type
            operator = ( const SCALAR_TYPE a )
            {
                static_cast< THIS& >( *this ).Fill( 0 );
                static_cast< THIS& >( *this ).Diagonal() = a;
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
                return SquarePlus< TA , NegativeTo< TB > >( m , v.Negative() );
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
                return SquarePlus< NegativeTo< TA > , TB >( m.Negative() , v );
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
                return SquarePlus< NegativeTo< TA > , ret_type< TA > >( m.Negative() , s );
                //return SquareMinus< TA , ret_type< TA > >( m , s );
            }
        };  // namespace Expression

        namespace Expression
        {
            // TA is square matrix;
            // TB is vector, or scalar
            template < typename TA , typename TB >
            struct SquarePlus :
                Expression< SquarePlus< TA , TB > >
            {
                constexpr const static int log2 = std::log2( Traits< TA >::cols + 1 );

                auto_creference< TA > a;
                auto_creference< TB > b;

                SquarePlus( auto_creference< TA > _a , auto_creference< TB > _b ) :
                    a( _a ) , b( _b )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr ) || HasSameRoot( b , ptr );
                }

                _ehm_inline ret_type< TA > operator [] ( IndexType i ) const
                {
                    return GetBy( a , i ) + ( (i&Traits< TA >::cols) == 0 )*GetBy( b , i >> log2 );
                }
                _ehm_inline ret_type< TA > Get ( IndexType x , IndexType y ) const
                {
                    return GetBy( a , x , y ) + ( x == y )*GetBy( b , 0 , x );
                }
            };
            template < typename TA , typename TB >
            struct Traits< SquarePlus< TA , TB > > : TraitsCombine< TA , TB >
            {
                // can use index if ( n + 1 ) is power of two
                constexpr const static bool is_gettable =
                    Traits< TA >::is_gettable || Traits< TB >::is_gettable || (Traits< TA >::cols&(Traits< TA >::cols+1))!=0;
                constexpr const static bool catch_reference = false;
            };
        };  // namespace Expression

    };  // namespace Matrix

};  //namespace E
