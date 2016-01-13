#pragma once

#include "Global.h"
#include "expression_traits.h"

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
        };  // namespace Util
    };  // namespace Matrix
};  // namespace EH
