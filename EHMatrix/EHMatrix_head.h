#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_Expression.h"
#include "Aliased_container.h"

#include <type_traits>
#include <initializer_list>

namespace EH
{
    namespace Matrix
    {
        template < typename THIS , IndexType M , IndexType N >
        struct Matrix_Interface
        {
            using typename THIS::result_type;

            _ehm_const IndexType rows = THIS::rows;
            _ehm_const IndexType cols = THIS::cols;

            // some operator implements here;
        };
    };  // namespace Matrix
};  //namespace EH
