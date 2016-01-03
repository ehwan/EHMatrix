#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_head.h"

#include <initializer_list>
#include <cmath>

namespace EH
{
    namespace Matrix
    {
        template < typename THIS , IndexType N >
        struct Matrix_Interface< THIS , N , 1 >
        {
            constexpr const static IndexType M = Expression::Traits< THIS >::rows;
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
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , M , 1 >& exp )
            {
                static_cast< THIS& >( *this ).Plus( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( exp ) );
            }
            template < typename CLS >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , M , 1 >& exp )
            {
                static_cast< THIS& >( *this ).Minus( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( exp ) );
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

            _ehm_inline void operator += ( const T a )
            {
                static_cast< THIS& >( *this ).Plus( a );
            }
            _ehm_inline void operator -= ( const T a )
            {
                static_cast< THIS& >( *this ).Minus( a );
            }

            template < typename CLS >
            _ehm_inline void operator *= ( const Expression::expression_size_type< CLS , M , 1 >& v )
            {
                static_cast< THIS& >( *this ).Multiply( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( v ) );
            }
            template < typename CLS >
            _ehm_inline void operator /= ( const Expression::expression_size_type< CLS , M , 1 >& v )
            {
                static_cast< THIS& >( *this ).Divide( typename Expression::AssignShouldMakeTemp< THIS , CLS >::tpye( v ) );
            }
            template < typename LST_TYPE >
            _ehm_inline
            void
            operator *= ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Multiply( lst.begin() , lst.end() );
            }
            template < typename LST_TYPE >
            _ehm_inline
            void
            operator /= ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Divide( lst.begin() , lst.end() );
            }

            template < typename EXP_CLS >
            _ehm_inline void operator = ( const Expression::Expression< EXP_CLS >& m )
            {
                static_cast< THIS& >( *this ).Fill( typename Expression::AssignShouldMakeTemp< THIS , EXP_CLS >::type( m ) );
            }
            template < typename LST_TYPE >
            _ehm_inline
            void
            operator = ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Fill( lst.begin() , lst.end() );
            }
            template < typename SCALAR_TYPE >
            _ehm_inline
            typename std::enable_if< Expression::is_scalar< SCALAR_TYPE >::value >::type
            operator = ( const SCALAR_TYPE a )
            {
                static_cast< THIS& >( *this ).Fill( a );
            }

            template < typename SFINE = THIS >
            typename std::enable_if< Expression::Traits< SFINE >::operations <= 1 , T >::type
            LengthSquared() const
            {
                using Expression::GetBy;
                T sum = T( 0 );
                for( IndexType i=0; i<M; ++i )
                {
                    sum += GetBy( static_cast< const THIS& >( *this ) , 0 , i ) * GetBy( static_cast< const THIS& >( *this ) , 0 , i );
                }
                return sum;
            }

            template < typename SFINE = THIS >
            typename std::enable_if< (Expression::Traits< SFINE >::operations > 1) , T >::type
            LengthSquared() const
            {
                using Expression::GetBy;
                T sum = T( 0 );
                for( IndexType i=0; i<M; ++i )
                {
                    const T t = GetBy( static_cast< const THIS& >( *this ) , 0 , i );
                    sum += t * t;
                }
                return sum;
            }
            auto Length() const
            {
                return std::sqrt( LengthSquared() );
            }
            auto Normalize()
            {
                using Expression::GetByRef;
                const auto l = Length();
                const auto invl = decltype( l )( 1 ) / l;
                for( IndexType i=0; i<M; ++i )
                {
                    GetByRef( static_cast< THIS& >( *this ) , 0 , i ) *= invl;
                }
                return l;
            }

            void Log() const
            {
                using Expression::GetBy;
                std::cout << "Vector Log\n";
                for( IndexType m=0; m<M; ++m )
                {
                    std::cout << "(\t";
                    std::cout << GetBy( static_cast< const THIS& >( *this ) , 0 , m );
                    std::cout << "\t)\n";
                }
            }
        };

        namespace Expression
        {
            // zero extra typename=void is for square-scalar operations
            template < typename TA , typename = void >
            auto
            operator + ( const expression_size_type< TA , 0 , 1 >& v ,
                         const ret_type< TA > s )
            {
                return Plus< TA , ret_type< TA > >( v , s );
            }
            template < typename TA , typename = void >
            auto
            operator - ( const expression_size_type< TA , 0 , 1 >& v ,
                         const ret_type< TA > s )
            {
                return Plus< TA , ret_type< TA > >( v , -s );
            }
            template < typename TA , typename = void >
            auto
            operator + ( const ret_type< TA > s ,
                         const expression_size_type< TA , 0 , 1 >& v )
            {
                return Plus< ret_type< TA > , TA >( s , v );
            }
            template < typename TA , typename = void >
            auto
            operator - ( const ret_type< TA > s ,
                         const expression_size_type< TA , 0 , 1 >& v )
            {
                return Plus< ret_type< TA > , NegativeTo< TA > >( s , v.Negative() );
            }

            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator * ( const expression_size_type< TA , M , 1 >& v1 ,
                         const expression_size_type< TB , M , 1 >& v2 )
            {
                return Multiply< TA , TB >( v1 , v2 );
            }
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator / ( const expression_size_type< TA , M , 1 >& v1 ,
                         const expression_size_type< TB , M , 1 >& v2 )
            {
                return Divide< TA , TB >( v1 , v2 );
            }
        };  // namespace Expression


        template < typename TA , typename TB , typename = void >
        typename std::common_type< Expression::ret_type< TA > , Expression::ret_type< TB > >::type
        Cross( const Expression::expression_size_type< TA , 2 , 1 >& v1 ,
               const Expression::expression_size_type< TB , 2 , 1 >& v2 )
        {
            using Expression::GetBy;
            return GetBy( v1 , 0 , 0 )*GetBy( v2 , 0 , 1 ) - GetBy( v1 , 0 , 1 )*GetBy( v2 , 0 , 0 );
        }
        template < typename TA >
        auto
        Cross( const Expression::expression_size_type< TA , 2 , 1 >& v ,
               const Expression::ret_type< TA > a )
        {
            return Expression::Multiply<
                        Expression::OuterProduct< 2 , TA > , Expression::ret_type< TA >
                        >( Expression::OuterProduct< 2 , TA >( v ) , a );
        }
        template < typename TA >
        auto
        Cross( const Expression::ret_type< TA > a ,
               const Expression::expression_size_type< TA , 2 , 1 >& v )
        {
            return Expression::Multiply<
                        Expression::OuterProduct< 2 , TA > , Expression::ret_type< TA >
                        >( Expression::OuterProduct< 2 , TA >( v ) , -a );
        }
        template < typename TA >
        auto
        Cross( const Expression::expression_size_type< TA , 2 , 1 >& v )
        {
            return Expression::OuterProduct< 2 , TA >( v );
        }
        template < typename TA , typename TB >
        auto
        Cross( const Expression::expression_size_type< TA , 3 , 1 >& _a ,
               const Expression::expression_size_type< TB , 3 , 1 >& _b )
        {
            using Expression::GetBy;
            typename Expression::ShouldMakeTemp< TA , 2 >::type a( _a );
            typename Expression::ShouldMakeTemp< TB , 2 >::type b( _b );
            return Vector< typename std::common_type< Expression::ret_type< TA > , Expression::ret_type< TB > >::type , 3 >
            {
                {
                    GetBy( a , 0 , 1 )*GetBy( b , 0 , 2 ) - GetBy( a , 0 , 2 )*GetBy( b , 0 , 1 ) ,
                    GetBy( a , 0 , 2 )*GetBy( b , 0 , 0 ) - GetBy( a , 0 , 0 )*GetBy( b , 0 , 2 ) ,
                    GetBy( a , 0 , 0 )*GetBy( b , 0 , 1 ) - GetBy( a , 0 , 1 )*GetBy( b , 0 , 0 )
                }
            };
        }

        namespace Expression
        {
            template < typename T >
            struct OuterProduct< 2 , T > :
                Expression< OuterProduct< 2 , T > >
            {
                auto_creference< T > a;

                OuterProduct( auto_creference< T > _a ) :
                    a( _a )
                {
                }

                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr );
                }
                //    0       1
                // { v[1] , -v[0] }
                _ehm_inline ret_type< T > operator [] ( IndexType i ) const
                {
                    return -GetBy( a , i==0 ) * MAKE_SIGNED( i );
                }
                _ehm_inline ret_type< T > Get ( IndexType x , IndexType i ) const
                {
                    return -GetBy( a , 0 , i==0 ) * MAKE_SIGNED( i );
                }
            };
            template < typename T >
            struct Traits< OuterProduct< 2 , T > > : Traits< T >
            {
                // vecctor type is always not-gettable
                constexpr const static IndexType cols = 1;
                constexpr const static IndexType rows = 2;
                constexpr const static bool is_restrict = false;
                constexpr const static int operations = Traits< T >::operations + 1;
                constexpr const static bool catch_reference = false;
            };
        }; // namespace Expression
    };  // namespace Matrix
};  // namespace EH
