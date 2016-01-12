#pragma once

#include "Global.h"
#include "expressions.h"
#include "expression_traits.h"

namespace EH
{
    namespace Matrix
    {
        namespace Complex
        {
            template < typename T >
            struct ComplexConjugate;
            template < typename TA , typename TB >
            struct ComplexMultiply;

            template < typename T = float , typename ARG >
            auto constexpr Complex( const ARG theta )
            {
                return Matrix< T , 2 , 1 >
                (
                    std::cos( theta ) , std::sin( theta )
                );
            }
            template < typename T >
            auto Conjugate( const expression_size_type< T , 2 , 1 >& c )
            {
                return ComplexConjugate< T >( c );
            }

            template < typename TA , typename TB >
            auto Multiply( const expression_size_type< TA , 2 , 1 >& _c1 ,
                           const expression_size_type< TB , 2 , 1 >& _c2 )
            {
                return ComplexMultiply< TA , TB >( _c1 , _c2 );
            }

            template < typename T >
            auto Arg( const expression_size_type< T , 2 , 1 >& c )
            {
                return std::atan2( GetBy( c , 0 , 1 ) , GetBy( c , 0 , 0 ) );
            }




            template < typename TA >
            struct ComplexConjugate :
                Expression< ComplexConjugate< TA > >
            {
                using typename Expression< ComplexConjugate< TA > >::result_type;

                auto_reference< TA > a;

                ComplexConjugate( auto_reference< TA > _a ) :
                    a( _a )
                {
                }

                _ehm_inline result_type Get( IndexType i ) const
                {
                    return -GetBy( a , i ) * MAKE_SIGNED( i );
                }
                _ehm_inline result_type Get( IndexType x , IndexType i ) const
                {
                    return -GetBy( a , 0 , i ) * MAKE_SIGNED( i );
                }
            };





            template < typename TA , typename TB >
            struct ComplexMultiply :
                Expression< ComplexMultiply< TA , TB > >
            {
                using typename Expression< ComplexMultiply< TA , TB > >::result_type;
                typename Expressions::ShouldMakeTemp< TA , 2 >::type a;
                typename Expressions::ShouldMakeTemp< TB , 2 >::type b;

                ComplexMultiply( auto_reference< TA > _a , auto_reference< TB > _b ) :
                    a( _a ) , b( _b )
                {
                }

                // a0b0 - a1b1
                // a0b1 + a1b0
                _ehm_inline result_type Get( IndexType i ) const
                {
                    return GetBy( a , 0 , 0 )*GetBy( b , 0 , i ) + MAKE_SIGNED( i )*GetBy( a , 0 , 1 )*GetBy( b , 0 , i==0 );
                }
                _ehm_inline result_type Get( IndexType x , IndexType i ) const
                {
                    return GetBy( a , 0 , 0 )*GetBy( b , 0 , i ) + MAKE_SIGNED( i )*GetBy( a , 0 , 1 )*GetBy( b , 0 , i==0 );
                }

            };


        };  // namespace Complex


        namespace Quaternion
        {
            template < typename TA >
            struct QuatConjugate;

            template < typename T >
            auto Quaternion( const expression_size_type< T , 3 , 1 >& v , typename expression_traits< T >::result_type w )
            {
                return Vector< typename expression_traits< T >::result_type , 4 >
                {
                    GetBy( v , 0 , 0 ) , GetBy( v , 0 , 1 ) , GetBy( v , 0 , 2 ) , w
                };
            }
            template < typename T >
            auto Quaternion( const Vector< T , 3 >& axis )
            {
                using result_type = typename expression_traits< T >::result_type;
                const result_type l = axis.Length();
                const result_type sinc =
                    l==result_type( 0 ) ?
                        result_type( 0 ) :
                        std::sin( l/2  )/l;
                return Quaternion( axis * sinc , std::cos( l/2 ) );
            }


            template < typename T >
            auto Conjugate( const expression_size_type< T , 4 , 1 >& q )
            {
                return QuatConjugate< T >( q );
            }

            template < typename TA , typename TB >
            auto Multiply( const expression_size_type< TA , 4 , 1 >& _q1 ,
                           const expression_size_type< TB , 4 , 1 >& _q2 )
            {
                const typename Expressions::ShouldMakeTemp< TA , 4 >::type q1( _q1 );
                const typename Expressions::ShouldMakeTemp< TB , 4 >::type q2( _q2 );

                return Matrix< typename expression_traits< TA , TB >::result_type , 4 , 1 >
                {
                    GetBy(q1,0,0)*GetBy(q2,0,3) + GetBy(q2,0,0)*GetBy(q1,0,3) + GetBy(q1,0,1)*GetBy(q2,0,2) - GetBy(q1,0,2)*GetBy(q2,0,1) ,
                    GetBy(q1,0,1)*GetBy(q2,0,3) + GetBy(q2,0,1)*GetBy(q1,0,3) + GetBy(q1,0,2)*GetBy(q2,0,0) - GetBy(q1,0,0)*GetBy(q2,0,2) ,
                    GetBy(q1,0,2)*GetBy(q2,0,3) + GetBy(q2,0,2)*GetBy(q1,0,3) + GetBy(q1,0,0)*GetBy(q2,0,1) - GetBy(q1,0,1)*GetBy(q2,0,0) ,
                   -GetBy(q1,0,0)*GetBy(q2,0,0) - GetBy(q1,0,1)*GetBy(q2,0,1) - GetBy(q1,0,2)*GetBy(q2,0,2) + GetBy(q1,0,3)*GetBy(q2,0,3)
                };
            }

            template < typename TA , typename TB , typename = void >
            auto Multiply( const expression_size_type< TA , 4 , 1 >& _q ,
                           const expression_size_type< TB , 3 , 1 >& _v )
            {
#define Q0 GetBy(q,0,0)
#define Q1 GetBy(q,0,1)
#define Q2 GetBy(q,0,2)
#define Q3 GetBy(q,0,3)
#define V0 GetBy(v,0,0)
#define V1 GetBy(v,0,1)
#define V2 GetBy(v,0,2)
                const typename Expressions::ShouldMakeTemp< TA , 8 >::type q( _q );
                const typename Expressions::ShouldMakeTemp< TB , 3 >::type v( _v );

                return Vector< typename expression_traits< TA , TB >::result_type , 3 >
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
            auto InverseMultiply( const expression_size_type< TA , 4 , 1 >& _q ,
                                  const expression_size_type< TB , 3 , 1 >& _v )
            {
#define Q0 GetBy(q,0,0)
#define Q1 GetBy(q,0,1)
#define Q2 GetBy(q,0,2)
#define Q3 GetBy(q,0,3)
#define V0 GetBy(v,0,0)
#define V1 GetBy(v,0,1)
#define V2 GetBy(v,0,2)
                const typename Expressions::ShouldMakeTemp< TA , 8 >::type q( _q );
                const typename Expressions::ShouldMakeTemp< TB , 3 >::type v( _v );

                return Vector< typename expression_traits< TA , TB >::result_type , 3 >
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


            template < typename CLS >
            auto Matrix( const expression_size_type< CLS , 4 , 1 >& _q )
            {
#define Q0 GetBy(q,0,0)
#define Q1 GetBy(q,0,1)
#define Q2 GetBy(q,0,2)
#define Q3 GetBy(q,0,3)
                const typename Expressions::ShouldMakeTemp< CLS , 8 >::type q( _q );

                return matrix_type< CLS >
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


            template < typename TA >
            struct QuatConjugate :
                Expression< QuatConjugate< TA > >
            {
                using typename Expression< QuatConjugate< TA > >::result_type;
                auto_reference< TA > a;

                QuatConjugate( auto_reference< TA > _a ) :
                    a( _a )
                {
                }

                _ehm_inline result_type Get( IndexType i ) const
                {
                    return GetBy( a , i ) * MAKE_SIGNED( i==3 );
                }
                _ehm_inline result_type Get( IndexType x , IndexType i ) const
                {
                    return GetBy( a , 0 , i ) * MAKE_SIGNED( i==3 );
                }
            };
        };  // namespace Quaternion


        template < typename TA >
        struct expression_traits< Complex::ComplexConjugate< TA > > : expression_traits< TA >
        {
            _ehm_const int operations = expression_traits< TA >::operations + 1;
            _ehm_const bool catch_reference = false;
        };
        template < typename TA , typename TB >
        struct expression_traits< Complex::ComplexMultiply< TA , TB > > : expression_traits< TA , TB >
        {
            using tempA = Expressions::ShouldMakeTemp< TA , 2 >;
            using tempB = Expressions::ShouldMakeTemp< TB , 2 >;
            _ehm_const bool is_single_index = true;
            _ehm_const bool is_restrict     = tempA::value && tempB::value;
            _ehm_const bool catch_reference = tempA::value || tempB::value;

            _ehm_const int operations = expression_traits< TA >::operations + expression_traits< TB >::operations + 2;
        };
        template < typename TA >
        struct expression_traits< Quaternion::QuatConjugate< TA > > : expression_traits< TA >
        {
            _ehm_const int operations = expression_traits< TA >::operations + 1;
            _ehm_const bool catch_reference = false;
        };

/*
            */
    };  // namespace Matrix
};  // namespace EH
