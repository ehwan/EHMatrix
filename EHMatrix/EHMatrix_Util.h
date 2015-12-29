#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_head.h"

namespace EH
{
    namespace Matrix
    {
        namespace Util
        {
            // returns matrix version of y = ax + b
            template < IndexType OutN = 0 ,
                       typename TA , typename TB ,
                       IndexType N = Expression::Traits< TA >::rows >
            auto FirstOrderMatrix( const Expression::expression_size_type< TA , N , 1 >& a ,
                                   const Expression::expression_size_type< TB , N , 1 >& b )
            {
                // if OutN is 0 , automatically fit the least number
                constexpr const IndexType Out = ( OutN==0 ? N+1 : OutN );
                static_assert( (Out > N) , "dimension of out matrix is not fit" );

                typedef typename std::common_type< Expression::ret_type< TA > , Expression::ret_type< TB > >::type retT;

                Matrix< retT , Out , Out > ret;
                ret.template Fill< retT >( 0 );

                auto& diag = ret.Diagonal();
                auto& col  = ret.Column( Out-1 );
                for( IndexType i=0; i<N; ++i )
                {
                    Expression::GetByRef( diag , 0 , i ) = Expression::GetBy( a , 0 , i );
                    Expression::GetByRef( col , 0 , i )  = Expression::GetBy( b , 0 , i );
                }

                return ret;
            }
            template < IndexType OutN = 0 ,
                       typename TA , typename TB , typename TC , typename TD ,
                       IndexType N = Expression::Traits< TA >::rows >
            auto OrthologMatrix( const Expression::expression_size_type< TA , N , 1 >& frommin ,
                                 const Expression::expression_size_type< TB , N , 1 >& frommax ,
                                 const Expression::expression_size_type< TC , N , 1 >& tomin ,
                                 const Expression::expression_size_type< TD , N , 1 >& tomax )
            {
                // y = ( x - from_a ) * ( to_b - to_a )/( from_b - from_a ) + to_a;
                const auto& coeff = ( tomax - tomin ) / ( frommax - frommin );
                const auto& offset = tomin - coeff * frommin;

                // now, y = x * coeff + offset;
                return FirstOrderMatrix< OutN >( coeff , offset );
            }
        };  // namespace Util
    };  // namespace Matrix
};  // namespace EH
