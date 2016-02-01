#pragma once

#include "Global.h"
#include "head.h"
#include <limits>
#include <numeric>

namespace EH
{
    template < typename T , Matrix::IndexType N >
    struct AABB
    {
        using vec_type = Matrix::Vector< T , N >;
        using aabb_type = AABB< T , N >;

        vec_type min;
        vec_type max;

        AABB()
        {
            //min = 0;
            //max = 0;
            min = std::numeric_limits< T >::max();
            max = std::numeric_limits< T >::lowest();
        }
        AABB( const vec_type& a , const vec_type& b ) :
            min( a ) , max( b )
        {
        }

        void operator |= ( const vec_type& v )
        {
            min = Matrix::min( v , min );
            max = Matrix::max( v , max );
        }
        void operator |= ( const aabb_type& aabb )
        {
            min = Matrix::min( min , aabb.min );
            max = Matrix::max( max , aabb.max );
        }
        void operator &= ( const aabb_type& aabb )
        {
            min = Matrix::max( min , aabb.min );
            max = Matrix::min( max , aabb.max );
        }

        aabb_type operator | ( const vec_type& v ) const
        {
            return aabb_type
            {
                Matrix::min( v , min ) ,
                Matrix::max( v , max )
            };
        }
        aabb_type operator | ( const aabb_type& aabb ) const
        {
            return aabb_type
            {
                Matrix::min( aabb.min , min ) ,
                Matrix::max( aabb.max , max )
            };
        }
        aabb_type operator & ( const aabb_type& aabb ) const
        {
            return aabb_type
            {
                Matrix::max( aabb.min , min ) ,
                Matrix::min( aabb.max , max )
            };
        }
    };
};
