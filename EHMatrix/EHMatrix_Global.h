#pragma once

#include <cassert>
#include <cmath>
#include <type_traits>
#include <iostream>
//#include "../EHLog.h"

#define MAKE_SIGNED(x) (int)( 1 | ( (x)-1 ) )

#ifndef NDEBUG
#define _ehm_inline
#else
#define _ehm_inline inline
#endif

#define _ehm_const constexpr const static

#define EXPRESSION_ASSIGN_OPERATOR(parent_type) \
    template < typename ASSIGN_TYPE >\
    _ehm_inline void operator = ( ASSIGN_TYPE&& exp )\
    {\
        parent_type::operator = ( std::template forward< ASSIGN_TYPE >( exp ) );\
    }

namespace EH
{
    namespace Matrix
    {
        using IndexType = unsigned int;

        //column-based M x N matrix
        // n columns
        template < typename T , IndexType M , IndexType N = M >
        struct Matrix;
        template < typename T , IndexType N >
        using Vector = Matrix< T , N , 1 >;

        template < typename CLS >
        struct Expression;

        template < typename CLS ,
                   IndexType M  ,
                   IndexType N >
        struct Matrix_Interface;
    };  // namespace Matrix;
};   // namespaec EH

