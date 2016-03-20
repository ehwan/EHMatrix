#pragma once

#include "Global.h"
#include "expression_traits.h"
#include "complex.h"

namespace EH
{
    namespace Matrix
    {
        namespace Util
        {
            // returns matrix version of y = ax + b
            template < IndexType OutN = 0 ,
                       typename TA , typename TB ,
                       IndexType N = expression_traits< TA >::rows >
            auto FirstOrderMatrix( const expression_size_type< TA , N , 1 >& a ,
                                   const expression_size_type< TB , N , 1 >& b )
            {
                // if OutN is 0 , automatically fit the least number
                constexpr const IndexType Out = ( OutN==0 ? N+1 : OutN );
                static_assert( (Out > N) , "dimension of out matrix is not fit" );

                using retT = typename expression_traits< TA , TB >::result_type;

                Matrix< retT , Out , Out > ret( 1 );

                for( IndexType i=0; i<N; ++i )
                {
                    RefBy( ret , i , i ) = GetBy( a , 0 , i );
                    RefBy( ret , Out-1 , i ) = GetBy( b , 0 , i );
                }

                return ret;
            }
            template < IndexType OutN = 0 ,
                       typename TA , typename TB , typename TC , typename TD ,
                       IndexType N >
            auto OrthologMatrix( const Vector< TA , N >& frommin ,
                                 const Vector< TB , N >& frommax ,
                                 const Vector< TC , N >& tomin ,
                                 const Vector< TD , N >& tomax )
            {
                // y = ( x - from_a ) * ( to_b - to_a )/( from_b - from_a ) + to_a;
                auto coeff = ( tomax - tomin ) / ( frommax - frommin );
                auto offset = tomin - coeff * frommin;

                // now, y = x * coeff + offset;
                return FirstOrderMatrix< OutN >( coeff , offset );
            }

            template < typename TA = float >
            auto PerspectiveMatrix( const Vector< float , 3 >& min ,
                                    const Vector< float , 3 >& max )
            {
                const Vector< float , 3 > sizeI = 1.0f / ( max - min );
                const Vector< float , 3 > mmid = ( min + max ) * sizeI;

                Matrix< float , 4 > ret =
                {
                    2 * min.z*sizeI.x , 0 , 0 , 0 ,
                    0 , 2 * min.z*sizeI.y , 0 , 0 ,
                    mmid.x , mmid.y , -mmid.z , -1 ,
                    0 , 0 , -2 * min.z*max.z*sizeI.z , 0
                };
                //mmid.z , 2 * sizei.Z;

                return ret;
            }
            template < typename TA = float >
            constexpr auto PerspectiveMatrix( float up , float down , float left , float right , float near , float far )
            {
                const float upt = std::tan( up );
                const float downt = std::tan( down );
                const float leftt = std::tan( left );
                const float rightt = std::tan( right );

                return PerspectiveMatrix( { -leftt * near , -downt * near , near } , { rightt * near , upt * near , far } );
            }

            // x pitch
            // y yaw
            // z roll
            template < typename TA = float >
            auto EyeMatrix( const Vector< float , 3 >& angle )
            {
                // pitch quaternion
                auto pq = Quaternion< float >( Vector< float , 3 >{ 1 , 0 , 0 } , angle[0] );
                // yaw quaternion
                auto yq = Quaternion< float >( Vector< float , 3 >{ 0 , 1 , 0 } , angle[1] );

                auto zq = yq * pq;

                auto zaxis = zq( { 0 , 0 , 1 } );
                // roll quaternion
                auto rq = Quaternion< float >( zaxis , angle[2] );

                auto xq = rq * yq;
                auto xaxis = xq( { 1 , 0 , 0 } );

                return Matrix< float , 3 >( xaxis , Cross( zaxis , xaxis ) , zaxis );
            }
            template < typename TA = float >
            auto EyeMatrix( const Matrix< float , 3 >& axis , const Vector< float , 3 >& pos )
            {
                const Matrix< float , 3 > tr = axis.Transpose();
                return Matrix< float , 4 >( tr , 0 , 0 , 0 , -tr * pos , 1 );
            }
        };  // namespace Util
    };  // namespace Matrix
}; //  namespace EH
