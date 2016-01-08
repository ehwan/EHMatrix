#pragma once

namespace EH
{
    namespace Matrix
    {
        template < typename T , IndexType M , IndexType N >
        struct Matrix_aliased_container
        {
            T s[ M*N ];
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
                T s[ M ];
                struct{ T x; T y; T z; T w; };
                struct{ T r; T g; T b; T a; };
            };
        };


        template < typename T >
        struct Matrix_aliased_container< T , 1 , 2 >
        {
            union
            {
                T s[ 2 ];
                struct{ T x; T y; };
                struct{ T r; T g; };
            };
        };
        template < typename T >
        struct Matrix_aliased_container< T , 1 , 3 >
        {
            union
            {
                T s[ 3 ];
                struct{ T x; T y; T z; };
                struct{ T r; T g; T b; };
            };
        };
        template < typename T , IndexType M >
        struct Matrix_aliased_container< T , 1 , M >
        {
            union
            {
                T s[ M ];
                struct{ T x; T y; T z; T w; };
                struct{ T r; T g; T b; T a; };
            };
        };

        template < typename T >
        struct Matrix_aliased_container< T , 1 , 1 >
        {
            union
            {
                T s[ 1 ];
                T x;
                T r;
            };
        };
    };  // namespace Matrix
};  // namespace EH
