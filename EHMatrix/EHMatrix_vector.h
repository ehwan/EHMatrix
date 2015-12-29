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
                static_cast< THIS& >( *this ).template Multiply< T >( a );
            }
            _ehm_inline void operator /= ( const T a )
            {
                static_cast< THIS& >( *this ).template Divide< T >( a );
            }

            template < typename CLS >
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , M , 1 >& exp )
            {
                static_cast< THIS& >( *this ).Plus_Safe( exp );
            }
            template < typename CLS >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , M , 1 >& exp )
            {
                static_cast< THIS& >( *this ).Minus_Safe( exp );
            }
            template < typename LST_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< LST_TYPE , T >::value >::type
            operator += ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Plus( lst.begin() , lst.end() );
            }
            template < typename LST_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< LST_TYPE , T >::value >::type
            operator -= ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Minus( lst.begin() , lst.end() );
            }

            _ehm_inline void operator += ( const T a )
            {
                static_cast< THIS& >( *this ).template Plus< T >( a );
            }
            _ehm_inline void operator -= ( const T a )
            {
                static_cast< THIS& >( *this ).template Minus< T >( a );
            }

            template < typename CLS >
            _ehm_inline void operator *= ( const Expression::expression_size_type< CLS , M , 1 >& v )
            {
                static_cast< THIS& >( *this ).Multiply_Safe( v );
            }
            template < typename CLS >
            _ehm_inline void operator /= ( const Expression::expression_size_type< CLS , M , 1 >& v )
            {
                static_cast< THIS& >( *this ).Divide_Safe( v );
            }
            template < typename LST_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< LST_TYPE , T >::value >::type
            operator *= ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Multiply( lst.begin() , lst.end() );
            }
            template < typename LST_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< LST_TYPE , T >::value >::type
            operator /= ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Divide( lst.begin() , lst.end() );
            }

            template < typename EXP_CLS >
            _ehm_inline void operator = ( const Expression::expression_size_type< EXP_CLS , M , 1 >& m )
            {
                static_cast< THIS& >( *this ).Fill_Safe( m );
            }
            template < typename LST_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< LST_TYPE , T >::value >::type
            operator = ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Fill( lst.begin() , lst.end() );
            }
            template < typename SCALAR_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< SCALAR_TYPE , T >::value >::type
            operator = ( const SCALAR_TYPE a )
            {
                static_cast< THIS& >( *this ).template Fill< T >( a );
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
                /*
                using Expression::GetBy;
                using EH::LOG::LOG;
                using EH::LOG::LOGR;
                LOG( "Vector Log" );
                for( IndexType m=0; m<M; ++m )
                {
                    LOGR( "(\t" );
                    LOGR( GetBy( static_cast< const THIS& >( *this ) , 0 , m ) );
                    LOGR( "\t)\n" );
                }
                */
            }
        };
        template < typename T , IndexType N >
        struct Vector_aliased_container
        {
            union
            {
                T s[ N ];
                struct{ T x; T y; T z; T w; };
                struct{ T r; T g; T b; T a; };
            };
        };
        template < typename T >
        struct Vector_aliased_container< T , 2 >
        {
            union
            {
                T s[ 2 ];
                struct{ T x; T y; };
                struct{ T r; T g; };
            };
        };
        template < typename T >
        struct Vector_aliased_container< T , 3 >
        {
            union
            {
                T s[ 3 ];
                struct{ T x; T y; T z; };
                struct{ T r; T g; T b; };
            };
        };
        template < typename T , IndexType M >
        struct Matrix< T , M , 1 > :
            Expression::Expression< Matrix< T , M , 1 > > ,
            Vector_aliased_container< T , M >
        {
            typedef Expression::Expression< Matrix< T , M , 1 > > parent;

            //T s[ M ];

            _ehm_inline T& operator [] ( IndexType i )
            {
                return Vector_aliased_container< T , M >::s[ i ];
            }
            _ehm_inline T& Get( IndexType x , IndexType y )
            {
                return Vector_aliased_container< T , M >::s[ y ];
            }
            _ehm_inline T operator [] ( IndexType i ) const
            {
                return Vector_aliased_container< T , M >::s[ i ];
            }
            _ehm_inline T Get( IndexType x , IndexType y ) const
            {
                return Vector_aliased_container< T , M >::s[ y ];
            }

            Matrix(){}
            template < typename CLS >
            Matrix( const Expression::expression_size_type< CLS , M , 1 >& m )
            {
                parent::template Fill< CLS >( m );
            }
            template < typename IterType >
            explicit Matrix( IterType begin , IterType end )
            {
                parent::Fill( begin , end );
            }
            Matrix( std::initializer_list< T > lst )
            {
                parent::Fill( lst.begin() , lst.end() );
            }
            Matrix( const T a )
            {
                parent::template Fill< T >( a );
            }


            EXPRESSION_ASSIGN_OPERATOR( parent )


            _ehm_inline T* data()
            {
                return Vector_aliased_container< T , M >::s;
            }
            _ehm_inline const T* data() const
            {
                return Vector_aliased_container< T , M >::s;
            }
            _ehm_inline T* begin()
            {
                return Vector_aliased_container< T , M >::s;
            }
            _ehm_inline T* end()
            {
                return Vector_aliased_container< T , M >::s + M;
            }
            _ehm_inline T* begin() const
            {
                return Vector_aliased_container< T , M >::s;
            }
            _ehm_inline T* end() const
            {
                return Vector_aliased_container< T , M >::s + M;
            }
            _ehm_inline T* cbegin() const
            {
                return Vector_aliased_container< T , M >::s;
            }
            _ehm_inline T* cend() const
            {
                return Vector_aliased_container< T , M >::s + M;
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
