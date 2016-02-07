#pragma once

#include "Global.h"
#include "expressions.h"
#include "expression_traits.h"

#include "head.h"

namespace EH
{
    namespace Matrix
    {
        template < typename T = float >
        struct Complex : Vector< T , 2 >
        {
            using Vector< T , 2 >::operator=;
            constexpr explicit Complex( T arg ) :
                Vector< T , 2 >( std::cos( arg ) , std::sin( arg ) )
            {
            }
            template < typename TA >
            Complex( const expression_vector_size_type< TA , 2 >& v ) :
                Vector< T , 2 >( GetByVector( v , 0 ) , GetByVector( v , 1 ) )
            {
            }
            Complex( T r , T i ) :
                Vector< T , 2 >( r , i )
            {
            }

            using Vector< T , 2 >::Get;
            inline auto Conjugate() const
            {
                return Complex< T >( Get( 0 ) , -Get( 1 ) );
            }
            auto Arg() const
            {
                return std::atan2( Get( 1 ) , Get( 0 ) );
            }
        };
        template < typename TA >
        auto operator * ( const Complex< TA >& c1 , const Vector< TA , 2 >& c2 )
        {
            return Complex< TA >( c1[0]*c2[0] - c1[1]*c2[1] , c1[0]*c2[1] + c1[1]*c2[0] );
        }

        template < typename T = float >
        struct Quaternion : Vector< T , 4 >
        {
            using Vector< T , 4 >::operator =;
            explicit Quaternion( const Vector< T , 3 >& axis )
            {
                using result_type = T;
                const result_type l = axis.Length();
                const result_type sinc =
                    l==result_type( 0 ) ?
                        result_type( 0 ) :
                        std::sin( l/2  )/l;
                this->FillAggressive( sinc * axis[0] , sinc * axis[1] , sinc * axis[2] , std::cos( l/2 ) );
            }
            template < typename TA , typename = void >
            Quaternion( const expression_vector_size_type< TA , 4 >& v ) :
                Vector< T , 4 >( GetByVector( v , 0 ) , GetByVector( v , 1 ) , GetByVector( v , 2 ) , GetByVector( v , 3 ) )
            {
            }
            template < typename TA >
            Quaternion( const expression_vector_size_type< TA , 3 >& v3 , T angle )
            {
                T c = std::cos( angle / 2 );
                T s = std::sin( angle / 2 );

                this->FillAggressive( s*GetByVector( v3 , 0 ) , s*GetByVector( v3 , 1 ) , s*GetByVector( v3 , 2 ) , c );
            }
            Quaternion( T v0 , T v1 , T v2 , T v3 ) :
                Vector< T , 4 >( v0 , v1 , v2 , v3 )
            {
            }

            using Vector< T , 4 >::Get;
            inline auto Conjugate() const
            {
                return Quaternion< T >( -Get( 0 ) , -Get( 1 ) , -Get( 2 ) , Get( 3 ) );
            }

            using Vector< T , 4 >::x;
            using Vector< T , 4 >::y;
            using Vector< T , 4 >::z;
            using Vector< T , 4 >::w;
            auto matrix() const
            {
                return Matrix< T , 3 >
                    (
                     x*x - y*y - z*z + w*w , 2 * ( w*z + x*y ) , 2 * ( x*z - y*w ) ,
                     2 * ( x*y - z*w ) , -x*x + y*y - z*z + w*w , 2 * ( x*w + y*z ) ,
                     2 * ( x*z + y*w ) , 2 * ( -x*w + y*z ) , -x*x - y*y + z*z + w*w
                    );
            }
            auto operator () ( const Vector< T , 3 >& v ) const
            {
                return Vector< T , 3 >
                    (
                        (x*x - y*y - z*z + w*w)*v.x + 2 * ( x*y - z*w )*v.y + 2 * ( x*z + y*w )*v.z ,
                        2 * ( w*z + x*y )*v.x + (-x*x + y*y - z*z + w*w)*v.y +  2 * ( -x*w + y*z )*v.z ,
                        2 * ( x*z - y*w )*v.x + 2 * ( x*w + y*z )*v.y + (-x*x - y*y + z*z + w*w)*v.z
                    );
            }
        };

        template < typename TA >
        auto operator * ( const Quaternion< TA >& q1 , const Vector< TA , 4 >& q2 )
        {
            return Quaternion< TA >
                (
                 q1[0]*q2[3] + q1[3]*q2[0] + q1[1]*q2[2] - q1[2]*q2[1] ,
                 q1[1]*q2[3] + q1[3]*q2[1] + q1[2]*q2[0] - q1[0]*q2[2] ,
                 q1[2]*q2[3] + q1[3]*q2[2] + q1[0]*q2[1] - q1[1]*q2[0] ,
                 -q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] + q1[3]*q2[3]
                );
        }

        /*
        namespace Complex
        {

            template < typename TA , typename TB >
            auto Multiply( const expression_size_type< TA , 2 , 1 >& _c1 ,
                           const expression_size_type< TB , 2 , 1 >& _c2 )
            {
                return ComplexMultiply< TA , TB >( _c1 , _c2 );
            }




        namespace Quaternion
        {
            template < typename TA >
            struct QuatConjugate;

            template < typename TA , typename TB >
            auto Multiply( const expression_size_type< TA , 4 , 1 >& _q1 ,
                           const expression_size_type< TB , 4 , 1 >& _q2 )
            {
                typename Expressions::ShouldMakeTemp< const TA , 4 >::type q1( _q1 );
                typename Expressions::ShouldMakeTemp< const TB , 4 >::type q2( _q2 );

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
                const typename Expressions::ShouldMakeTemp< const TA , 8 >::type q( _q );
                typename Expressions::ShouldMakeTemp< const TB , 3 >::type v( _v );

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
            */
    };  // namespace Matrix
};  // namespace EH
