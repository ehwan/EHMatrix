#pragma once

#include <cassert>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <utility>

// the external libraries

#include "../../EHLog.h"
#include "../../static_sequence/static_sequence.h"

#define MAKE_SIGNED(x) (int)( 1 | ( (x)-1 ) )

#ifndef NDEBUG
#define _ehm_inline
#else
#define _ehm_inline inline
#endif

#define _ehm_const constexpr const static

namespace EH
{
    namespace Matrix
    {
        using EH::LOG;
        using EH::ERROR;
        using EH::LOGR;

        using IndexType = unsigned int;

        //column-based M x N matrix
        // n columns
        template < typename T , IndexType M , IndexType N = M >
        struct Matrix;
        template < typename T , IndexType N >
        using Vector = Matrix< T , N , 1 >;

        template < typename T >
        using vec2 = Vector< T , 2 >;
        template < typename T >
        using vec3 = Vector< T , 3 >;
        template < typename T >
        using vec4 = Vector< T , 4 >;

        template < typename T >
        using mat2 = Matrix< T , 2 >;
        template < typename T >
        using mat3 = Matrix< T , 3 >;
        template < typename T >
        using mat4 = Matrix< T , 4 >;

        using vec2f = vec2< float >;
        using vec3f = vec3< float >;
        using vec4f = vec4< float >;
        using mat2f = mat2< float >;
        using mat3f = mat3< float >;
        using mat4f = mat4< float >;

        template < typename CRTP >
        struct Expression;
        template < typename CRTP >
        struct WritableExpression;

        template < typename ... Ts >
        struct expression_traits;
    };  // namespace Matrix;
};   // namespaec EH

