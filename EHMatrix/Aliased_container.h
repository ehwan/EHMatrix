#pragma once

namespace EH
{
    namespace Matrix
    {
        template < typename T , IndexType M , IndexType N >
        struct Matrix_aliased_container
        {
            T s[ N ];
        };
        template < typename T >
        struct Matrix_aliased_container< T , 2 , 1 >
        {
            union
            {
                T s[ 2 ];
                struct{ T x; T y; };
                struct{ T r; T g; };
            };
        };
        template < typename T >
        struct Matrix_aliased_container< T , 3 , 1 >
        {
            union
            {
                T s[ 3 ];
                struct{ T x; T y; T z; };
                struct{ T r; T g; T b; };
            };
        };
        template < typename T , IndexType M >
        struct Matrix_aliased_container< T , M , 1 >
        {
            union
            {
                T s[ 4 ];
                struct{ T x; T y; T z; T w; };
                struct{ T r; T g; T b; T a; };
            };
        };
    };  // namespace Matrix
};  // namespace EH
