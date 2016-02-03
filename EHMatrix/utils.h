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
                       IndexType N = expression_traits< TA >::rows >
            auto OrthologMatrix( const expression_size_type< TA , N , 1 >& frommin ,
                                 const expression_size_type< TB , N , 1 >& frommax ,
                                 const expression_size_type< TC , N , 1 >& tomin ,
                                 const expression_size_type< TD , N , 1 >& tomax )
            {
                // y = ( x - from_a ) * ( to_b - to_a )/( from_b - from_a ) + to_a;
                const auto coeff = ( tomax - tomin ) / ( frommax - frommin );
                const auto offset = tomin - coeff * frommin;

                // now, y = x * coeff + offset;
                return FirstOrderMatrix< OutN >( coeff , offset );
            }

            template < typename TA , typename TB >
            auto PerspectiveMatrix( const expression_size_type< TA , 3 , 1 >& min ,
                                    const expression_size_type< TB , 3 , 1 >& max ,
                                    float nearl )
            {
                const Vector< float , 3 > sizeI = 1.0f / ( max - min );
                const Vector< float , 3 > mmid = ( min + max ) * sizeI;

                Matrix< float , 4 > ret =
                {
                    2 * nearl * sizeI.x , 0 , 0 , 0 ,
                    0 , 2 * nearl * sizeI.y , 0 , 0 ,
                    mmid.x , mmid.y , -mmid.z , -1 ,
                    0 , 0 , 2 * sizeI.z , -2 * min[2] * max[2] * sizeI.z , 0
                };
                //mmid.z , 2 * sizei.Z;

                return ret;
            }

            // x pitch
            // y yaw
            // z roll
            template < typename TA >
            auto EyeMatrix( const expression_size_type< TA , 3 , 1 >& angle )
            {
                const Vector< float , 4 > pq = Quaternion::Quaternion( Vector< float , 3 >( 1 , 0 , 0 ) , angle[0] );
                const Vector< float , 4 > yq = Quaternion::Quaternion( Vector< float , 3 >( 0 , 1 , 0 ) , angle[1] );
                const Vector< float , 4 > zq = Quaternion::Multiply( pq , yq );

                const Vector< float , 3 > zaxis = Quaternion::Multiply( zq , Vector< float , 3 >( 0 , 0 , 1 ) );
                const Vector< float , 4 > rq = Quaternion::Quaternion( zaxis , angle[2] );

                const Vector< float , 4 > xq = Quaternion::Multiply( rq , yq );
                //const Vector< float , 3 > xaxis = Quaternion::Multiply( xq , Vector< float , 3 >( 1 , 0 , 0 ) );
            }
        };  // namespace Util
    };  // namespace Matrix
}; //  namespace EH
