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
        using EH::LOG::LOG;
        using EH::LOG::ERROR;
        using EH::LOG::LOGR;

        using IndexType = unsigned int;

        //column-based M x N matrix
        // n columns
        template < typename T , IndexType M , IndexType N = M >
        struct Matrix;
        template < typename T , IndexType N >
        using Vector = Matrix< T , N , 1 >;

        template < typename CRTP >
        struct Expression;
        template < typename CRTP >
        struct WritableExpression;

        template < typename ... Ts >
        struct expression_traits;
    };  // namespace Matrix;
};   // namespaec EH

