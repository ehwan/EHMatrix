#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_head.h"
#include <cmath>

namespace EH
{
    namespace Matrix
    {
        namespace Complex
        {
            template < typename T = float , typename ARG >
            auto constexpr Complex( const ARG theta )
            {
                return Matrix< T , 2 , 1 >
                {
                    std::cos( theta ) , std::sin( theta )
                };
            }

            template < typename TA , typename TB >
            auto Multiply( const Expression::expression_size_type< TA , 2 , 1 >& _c1 ,
                           const Expression::expression_size_type< TB , 2 , 1 >& _c2 )
            {
                return Expression::ComplexMultiply< TA , TB >( _c1 , _c2 );
            }

            template < typename T >
            auto Conjugate( const Expression::expression_size_type< T , 2 , 1 >& c )
            {
                return Expression::ComplexConjugate< T >( c );
            }

            template < typename T >
            auto Arg( const Expression::expression_size_type< T , 2 , 1 >& c )
            {
                return std::atan2( Expression::GetBy( c , 0 , 1 ) , Expression::GetBy( c , 0 , 0 ) );
            }
        };  // namespae Complex

        namespace Quaternion
        {
            template < typename T >
            auto Quaternion( const Expression::expression_size_type< T , 3 , 1 >& v , const Expression::ret_type< T > w )
            {
                return Vector< Expression::ret_type< T > , 4 >
                {
                    { GetBy( v , 0 , 0 ) , GetBy( v , 0 , 1 ) , GetBy( v , 0 , 2 ) , w }
                };
            }
            template < typename T >
            auto Quaternion( const Vector< T , 3 >& axis )
            {
                const Expression::ret_type< T > l = axis.Length();
                const Expression::ret_type< T > sinc =
                    l==Expression::ret_type< T >(0) ?
                        Expression::ret_type< T >(0) :
                        std::sin( l*Expression::ret_type< T >(0.5) )/l;
                return Quaternion( axis * sinc , std::cos( l*Expression::ret_type< T >(0.5) ) );
            }


            template < typename TA , typename TB >
            auto Multiply( const Expression::expression_size_type< TA , 4 , 1 >& q1 ,
                           const Expression::expression_size_type< TB , 4 , 1 >& q2 )
            {
                return Expression::QuatMult< TA , TB >( q1 , q2 );
            }
            template < typename TA , typename TB , typename = void >
            auto Multiply( const Expression::expression_size_type< TA , 4 , 1 >& _q ,
                           const Expression::expression_size_type< TB , 3 , 1 >& _v )
            {
#define Q0 Expression::GetBy(q,0,0)
#define Q1 Expression::GetBy(q,0,1)
#define Q2 Expression::GetBy(q,0,2)
#define Q3 Expression::GetBy(q,0,3)
#define V0 Expression::GetBy(v,0,0)
#define V1 Expression::GetBy(v,0,1)
#define V2 Expression::GetBy(v,0,2)
                typename Expression::ShouldMakeTemp< TA , 8 >::type q( _q );
                typename Expression::ShouldMakeTemp< TB , 3 >::type v( _v );

                return Vector< typename std::common_type< Expression::ret_type< TA > , Expression::ret_type< TB > >::type , 3 >
                {
                    V0 * ( Q0*Q0 - Q1*Q1 - Q2*Q2 + Q3*Q3 ) +
                    V1 * 2*( Q0*Q1 + Q2*Q3 ) +
                    V2 * 2*( Q0*Q2 - Q1*Q3 ) ,

                    V0 * 2*( Q0*Q1 - Q2*Q3 ) +
                    V1 * ( -Q0*Q0 + Q1*Q1 - Q2*Q2 + Q3*Q3 ) +
                    V2 * 2*( Q0*Q3 + Q1*Q2 ) ,

                    V0 * 2*( Q0*Q2 + Q1*Q3 ) +
                    V1 * 2*( -Q0*Q3 + Q1*Q2 ) +
                    V2 * ( -Q0*Q0 - Q1*Q1 + Q2*Q2 + Q3*Q3 )
                };
#undef Q0
#undef Q1
#undef Q2
#undef Q3
#undef V0
#undef V1
#undef V2
            }
            template < typename TA , typename TB >
            auto InverseMultiply( const Expression::expression_size_type< TA , 4 , 1 >& _q ,
                                  const Expression::expression_size_type< TB , 3 , 1 >& _v )
            {
#define Q0 Expression::GetBy(q,0,0)
#define Q1 Expression::GetBy(q,0,1)
#define Q2 Expression::GetBy(q,0,2)
#define Q3 Expression::GetBy(q,0,3)
#define V0 Expression::GetBy(v,0,0)
#define V1 Expression::GetBy(v,0,1)
#define V2 Expression::GetBy(v,0,2)
                typename Expression::ShouldMakeTemp< TA , 8 >::type q( _q );
                typename Expression::ShouldMakeTemp< TB , 3 >::type v( _v );

                return Vector< typename std::common_type< Expression::ret_type< TA > , Expression::ret_type< TB > >::type , 3 >
                {
                    V0 * ( Q0*Q0 - Q1*Q1 - Q2*Q2 + Q3*Q3 ) +
                    V1 * 2*( Q0*Q1 - Q2*Q3 ) +
                    V2 * 2*( Q0*Q2 + Q1*Q3 ) ,

                    V0 * 2*( Q0*Q1 + Q2*Q3 ) +
                    V1 * ( -Q0*Q0 + Q1*Q1 - Q2*Q2 + Q3*Q3 ) +
                    V2 * 2*( -Q0*Q3 + Q1*Q2 ) ,

                    V0 * 2*( Q0*Q2 - Q1*Q3 ) +
                    V1 * 2*( Q0*Q3 + Q1*Q2 ) +
                    V2 * ( -Q0*Q0 - Q1*Q1 + Q2*Q2 + Q3*Q3 )
                };
#undef Q0
#undef Q1
#undef Q2
#undef Q3
#undef V0
#undef V1
#undef V2
            }
            template < typename T >
            auto Conjugate( const Expression::expression_size_type< T , 4 , 1 >& q )
            {
                return Expression::QuatConjugate< T >( q );
            }
            template < typename CLS >
            auto Matrix( const Expression::expression_size_type< CLS , 4 , 1 >& _q )
            {
#define Q0 Expression::GetBy(q,0,0)
#define Q1 Expression::GetBy(q,0,1)
#define Q2 Expression::GetBy(q,0,2)
#define Q3 Expression::GetBy(q,0,3)
                typename Expression::ShouldMakeTemp< CLS , 8 >::type q( _q );

                return EH::Matrix::Matrix< Expression::ret_type< CLS > , 3 , 3 >
                {
                    ( Q0*Q0 - Q1*Q1 - Q2*Q2 + Q3*Q3 ) ,
                    2*( Q0*Q1 - Q2*Q3 ) ,
                    2*( Q0*Q2 + Q1*Q3 ) ,

                    2*( Q0*Q1 + Q2*Q3 ) ,
                    ( -Q0*Q0 + Q1*Q1 - Q2*Q2 + Q3*Q3 ) ,
                    2*( -Q0*Q3 + Q1*Q2 ) ,

                    2*( Q0*Q2 - Q1*Q3 ) ,
                    2*( Q0*Q3 + Q1*Q2 ) ,
                    ( -Q0*Q0 - Q1*Q1 + Q2*Q2 + Q3*Q3 )
                };
#undef Q0
#undef Q1
#undef Q2
#undef Q3
            }

        };  // namespace Quaternion
        namespace Expression
        {
            template < typename TA >
            struct ComplexConjugate :
                Expression< ComplexConjugate< TA > >
            {
                auto_creference< TA > a;

                ComplexConjugate( auto_creference< TA > _a ) :
                    a( _a )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr );
                }

                inline ret_type< TA > operator [] ( IndexType i ) const
                {
                    return -GetBy( a , i ) * MAKE_SIGNED( i );
                }
                inline ret_type< TA > Get( IndexType x , IndexType i ) const
                {
                    return -GetBy( a , 0 , i ) * MAKE_SIGNED( i );
                }
            };
            template < typename TA >
            struct Traits< ComplexConjugate< TA > > : Traits< TA >
            {
                constexpr const static IndexType cols = 1;
                constexpr const static IndexType rows = 2;
                constexpr const static int operations = Traits< TA >::operations + 1;
                constexpr const static bool catch_reference = false;
            };
            template < typename TA , typename TB >
            struct ComplexMultiply :
                Expression< ComplexMultiply< TA , TB > >
            {
                typename ShouldMakeTemp< TA , 2 >::type a;
                typename ShouldMakeTemp< TB , 2 >::type b;

                ComplexMultiply( auto_creference< TA > _a , auto_creference< TB > _b ) :
                    a( _a ) , b( _b )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr ) || HasSameRoot( b , ptr );
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >* ptr ) const
                {
                    return a.has_same_root( ptr ) || b.has_same_root( ptr );
                }

                // a0b0 - a1b1
                // a0b1 + a1b0
                inline typename std::common_type< ret_type< TA > , ret_type< TB > >::type operator [] ( IndexType i ) const
                {
                    return GetBy( a , 0 , 0 )*GetBy( b , 0 , i ) + MAKE_SIGNED( i )*GetBy( a , 0 , 1 )*GetBy( b , 0 , i==0 );
                }
                inline typename std::common_type< ret_type< TA > , ret_type< TB > >::type Get( IndexType x , IndexType i ) const
                {
                    return GetBy( a , 0 , 0 )*GetBy( b , 0 , i ) + MAKE_SIGNED( i )*GetBy( a , 0 , 1 )*GetBy( b , 0 , i==0 );
                }

            };
            template < typename TA , typename TB >
            struct Traits< ComplexMultiply< TA , TB > > : TraitsCombine< TA , TB >
            {
                constexpr const static bool is_gettable =
                    ShouldMakeTemp< TA , 2 >::is_gettable || ShouldMakeTemp< TB , 2 >::is_gettable;
                constexpr const static IndexType cols = 1;
                constexpr const static IndexType rows = 2;
                typedef typename std::common_type< ret_type< TA > , ret_type< TB > >::type return_type;
                constexpr const static bool is_restrict =
                    ShouldMakeTemp< TA , 2 >::value && ShouldMakeTemp< TB , 2 >::value;
                constexpr const static int operations =
                    ( ShouldMakeTemp< TA , 2 >::operations + ShouldMakeTemp< TB , 2 >::operations ) * 2 + 2;
                constexpr const static bool catch_reference =
                    ShouldMakeTemp< TA , 2 >::value || ShouldMakeTemp< TB , 2 >::value;

            };

            template < typename TA >
            struct QuatConjugate :
                Expression< QuatConjugate< TA > >
            {
                auto_creference< TA > a;

                QuatConjugate( auto_creference< TA > _a ) :
                    a( _a )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr );
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >* ptr ) const
                {
                    return a.has_same_root( ptr );
                }

                inline ret_type< TA > operator [] ( IndexType i ) const
                {
                    return GetBy( a , i ) * MAKE_SIGNED( i==3 );
                }
                inline ret_type< TA > Get( IndexType x , IndexType i ) const
                {
                    return GetBy( a , 0 , i ) * MAKE_SIGNED( i==3 );
                }
            };
            template < typename TA >
            struct Traits< QuatConjugate< TA > > : Traits< TA >
            {
                constexpr const static IndexType cols = 1;
                constexpr const static IndexType rows = 4;
                constexpr const static int operations = Traits< TA >::operations + 1;
                constexpr const static bool catch_reference = false;
            };

            template < typename TA , typename TB >
            struct QuatMult :
                Expression< QuatMult< TA , TB > >
            {
                //common_type< ret_type< TA > , ret_type< TB > > outs[4];
                typename ShouldMakeTemp< TA , 4 >::type a;
                typename ShouldMakeTemp< TB , 4 >::type b;

                QuatMult( auto_creference< TA > _a , auto_creference< TB > _b ) :
                    a( _a ) , b( _b )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr ) || HasSameRoot( b , ptr );
                }

                //
                // + a0b3 + a1b2 - a2b1 + a3b0
                // - a0b2 + a1b3 + a2b0 + a3b1
                // + a0b1 - a1b0 + a2b3 + a3b2
                // - a0b0 - a1b1 - a2b2 + a3b3

                inline typename std::common_type< ret_type< TA > , ret_type< TB > >::type operator [] ( IndexType i ) const
                {
                    return
                        -MAKE_SIGNED( i&1 )         *    GetBy( a , 0 , 0 )*GetBy( b , 0 , i^0b11 ) -
                        (int)( ( i & 0b10 ) - 1 )   *    GetBy( a , 0 , 1 )*GetBy( b , 0 , i^0b10 ) +
                        (int)( ((i+1)&0b10) - 1 )   *    GetBy( a , 0 , 2 )*GetBy( b , 0 , i^0b01 ) +
                                                         GetBy( a , 0 , 3 )*GetBy( b , 0 , i );
                }
                inline typename std::common_type< ret_type< TA > , ret_type< TB > >::type Get( IndexType x , IndexType i ) const
                {
                    return
                        -MAKE_SIGNED( i&1 )         *    GetBy( a , 0 , 0 )*GetBy( b , 0 , i^0b11 ) -
                        (int)( ( i & 0b10 ) - 1 )   *    GetBy( a , 0 , 1 )*GetBy( b , 0 , i^0b10 ) +
                        (int)( ((i+1)&0b10) - 1 )   *    GetBy( a , 0 , 2 )*GetBy( b , 0 , i^0b01 ) +
                                                         GetBy( a , 0 , 3 )*GetBy( b , 0 , i );
                }
            };
            template < typename TA , typename TB >
            struct Traits< QuatMult< TA , TB > > : TraitsCombine< TA , TB >
            {
                constexpr const static bool is_gettable =
                    ShouldMakeTemp< TA , 4 >::is_gettable || ShouldMakeTemp< TB , 4 >::is_gettable;
                constexpr const static IndexType cols = 1;
                constexpr const static IndexType rows = 4;
                typedef typename std::common_type< ret_type< TA > , ret_type< TB > >::type return_type;
                constexpr const static bool is_restrict =
                    ShouldMakeTemp< TA , 4 >::value && ShouldMakeTemp< TB , 4 >::value;
                constexpr const static int operations =
                    ( ShouldMakeTemp< TA , 4 >::operations + ShouldMakeTemp< TB , 4 >::operations )*4 + 4;
                constexpr const static bool catch_reference =
                    ShouldMakeTemp< TA , 4 >::value || ShouldMakeTemp< TB , 4 >::value;
            };

        };  // namespace Expression
    };  // namespace Matrix
};  // namespace EH
