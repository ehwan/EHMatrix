#pragma once

#include <cassert>
#include <cmath>
#include <type_traits>
#include <iostream>
//#include "../EHLog.h"

//#define FMA_OPTIMIZATION 1

#define MAKE_SIGNED(x) (int)( 1 | ( (x)-1 ) )

#define ENABLE_ERROR_EXPRESSION true

#ifndef NDEBUG
#define _ehm_inline
#else
#define _ehm_inline inline
#endif

#define EXPRESSION_ASSIGN_EXPRESSION(parent_type) \
    template < typename EXP_CLS >\
    _ehm_inline void operator = ( const EH::Matrix::Expression::Expression< EXP_CLS >& m )\
    {\
        parent_type::operator = ( m );\
    }
#define EXPRESSION_ASSIGN_LST(parent_type) \
    template < typename LST_TYPE >\
    _ehm_inline void operator = ( std::initializer_list< LST_TYPE > lst )\
    {\
        parent_type::operator = ( lst );\
    }
#define EXPRESSION_ASSIGN_SCALAR(parent_type) \
    template < typename SCALAR_TYPE >\
    _ehm_inline \
    typename std::enable_if< std::is_arithmetic< SCALAR_TYPE >::value >::type\
    operator = ( const SCALAR_TYPE _scalar_ )\
    {\
        parent_type::operator = ( _scalar_ );\
    }
#define EXPRESSION_ASSIGN_OPERATOR(parent_type) \
    EXPRESSION_ASSIGN_EXPRESSION(parent_type) \
    EXPRESSION_ASSIGN_LST(parent_type) \
    EXPRESSION_ASSIGN_SCALAR(parent_type)

namespace EH
{
    namespace Matrix
    {
        typedef unsigned int IndexType;
        //column-based M x N matrix
        // n columns
        template < typename T , IndexType M , IndexType N = M >
        struct Matrix;
        template < typename T , IndexType N >
        using Vector = Matrix< T , N , 1 >;

        namespace Expression
        {
            template< typename EXP >
            struct Traits;
            template < typename CLS >
            struct Expression;

        };

        template < typename CLS ,
                   IndexType M = Expression::Traits< CLS >::rows ,
                   IndexType N = Expression::Traits< CLS >::cols >
        struct Matrix_Interface;

        template < typename T >
        using remove_cr = typename std::remove_reference< typename std::remove_const< T >::type >::type;

        template< typename T >
        using auto_reference = typename std::conditional<
              std::is_arithmetic< T >::value || Expression::Traits< T >::catch_reference==false ,
              T ,
              typename std::add_lvalue_reference< T >::type
          >::type;
        template < typename T >
        using auto_creference = auto_reference< typename std::add_const< T >::type >;

    };  // namespace Matrix;
};   // namespaec EH

