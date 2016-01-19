#pragma once

#include "Global.h"

namespace EH
{
    namespace Matrix
    {
        template < typename T >
        struct is_scalar :
            std::is_arithmetic< typename std::remove_reference< T >::type >
        {
        };
        template < typename T >
        struct is_column_vector :
            std::integral_constant< bool , expression_traits< T >::cols == 1 >
        {
        };
        template < typename T >
        struct is_row_vector :
            std::integral_constant< bool , expression_traits< T >::rows == 1 >
        {
        };
        template < typename T >
        struct is_vector :
            std::integral_constant< bool , is_column_vector< T >::value || is_row_vector< T >::value >
        {
        };
        template < typename T >
        struct is_square :
            std::integral_constant< bool , expression_traits< T >::rows == expression_traits< T >::cols >
        {
        };
        template < typename T >
        struct is_expression :
            std::integral_constant< bool , expression_traits< T >::is_expression >
        {
        };
        template < typename TA , typename TB >
        struct is_same_cols :
            std::integral_constant<
                bool ,
                expression_traits< TA >::cols == expression_traits< TB >::cols
            >
        {
        };
        template < typename TA , typename TB >
        struct is_same_rows :
            std::integral_constant<
                bool ,
                expression_traits< TA >::rows == expression_traits< TB >::rows
            >
        {
        };
        template < typename TA , typename TB >
        struct is_same_size :
            std::integral_constant<
                bool ,
                is_same_cols< TA , TB >::value && is_same_rows< TA , TB >::value
            >
        {
        };
        template < typename TA , typename TB >
        struct is_same_vector :
            std::integral_constant<
                bool ,
                is_vector< TA >::value && is_vector< TB >::value &&
                is_column_vector< TA >::value == is_column_vector< TB >::value &&
                is_row_vector< TA >::value == is_row_vector< TB >::value
            >
        {
        };
        template < typename TA >
        struct vector_size :
            std::integral_constant<
                IndexType ,
                std::max( expression_traits< TA >::rows , expression_traits< TA >::cols )
            >
        {
        };

        template < typename DST , typename SRC >
        struct is_assign_restrict :
            std::integral_constant<
                bool ,
                ( expression_traits< DST >::is_restrict && expression_traits< SRC >::is_restrict )
                //|| ( std::is_same< expression_traits< DST >::root_type , expression_traits< SRC >::root_type >::value == false )
            >
        {
        };

        template < typename T >
        using matrix_type = Matrix< typename expression_traits< T >::result_type ,
                                    expression_traits< T >::rows ,
                                    expression_traits< T >::cols >;


        template < typename T >
        struct expression_traits< T >
        {
            using result_type = typename std::remove_const< typename std::remove_reference< T >::type >::type;
            using root_type   = result_type;

            _ehm_const IndexType cols            = 1;
            _ehm_const IndexType rows            = 1;

            _ehm_const bool      is_single_index = true;
            _ehm_const bool      is_restrict     = true;
            _ehm_const bool      catch_reference = false;

            _ehm_const bool      is_expression   = false;

            _ehm_const int       operations      = 0;


        };
        template < typename TA , typename TB >
        struct expression_traits< TA , TB >
        {
            using result_type = typename std::common_type< typename expression_traits< TA >::result_type ,
                                                           typename expression_traits< TB >::result_type >::type;

            _ehm_const IndexType cols = std::max( expression_traits< TA >::cols , expression_traits< TB >::cols );
            _ehm_const IndexType rows = std::max( expression_traits< TA >::rows , expression_traits< TB >::rows );

            _ehm_const bool is_single_index = expression_traits< TA >::is_single_index && expression_traits< TB >::is_single_index;
            _ehm_const bool is_restrict     = expression_traits< TA >::is_restrict && expression_traits< TB >::is_restrict;

            _ehm_const bool is_expression   = true;
        };
        template < typename T >
        struct expression_traits< Expression< T > > : expression_traits< T >
        {
        };
        template < typename T >
        struct expression_traits< WritableExpression< T > > : expression_traits< T >
        {
        };

        template < typename T >
        struct expression_traits< const T > : expression_traits< T >
        {
        };
        template < typename T >
        struct expression_traits< T& > : expression_traits< T >
        {
        };
        template < typename T >
        struct expression_traits< T const& > : expression_traits< T >
        {
        };



        template < typename TA >
        using auto_reference = typename std::conditional<
                                    expression_traits< TA >::catch_reference ,
                                    typename std::add_lvalue_reference< TA >::type ,
                                    typename std::remove_reference< TA >::type
                                >::type;


        template < typename TA >
        typename std::enable_if< is_scalar< TA >::value , typename expression_traits< TA >::result_type >::type
        constexpr _ehm_inline
        GetBy( TA&& a , IndexType x , IndexType y )
        {
            return a;
        }
        template < typename TA >
        typename std::enable_if< is_scalar< TA >::value , typename expression_traits< TA >::result_type >::type
        constexpr _ehm_inline
        GetBy( TA&& a , IndexType i )
        {
            return a;
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            expression_traits< TA >::is_single_index ,
            typename expression_traits< TA >::result_type
        >::type
        constexpr _ehm_inline
        GetBy( TA&& a , IndexType i )
        {
            return a.Get( i );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            expression_traits< TA >::is_single_index == false ,
            typename expression_traits< TA >::result_type
        >::type
        constexpr _ehm_inline
        GetBy( TA&& a , IndexType i )
        {
            ERROR( "Single-index acces is not valid for this expression" );
            return a.Get( i );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            expression_traits< TA >::is_single_index == false ,
            typename expression_traits< TA >::result_type
        >::type
        constexpr _ehm_inline
        GetBy( TA&& a , IndexType x , IndexType y )
        {
            return a.Get( x , y );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            expression_traits< TA >::is_single_index ,
            typename expression_traits< TA >::result_type
        >::type
        constexpr _ehm_inline
        GetBy( TA&& a , IndexType x , IndexType y )
        {
            return a.Get( y + x * expression_traits< TA >::rows );
        }

        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            is_column_vector< TA >::value ,
            typename expression_traits< TA >::result_type
        >::type
        constexpr _ehm_inline
        GetByVector( TA&& a , IndexType i )
        {
            return GetBy( a , 0 , i );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            is_row_vector< TA >::value ,
            typename expression_traits< TA >::result_type
        >::type
        constexpr _ehm_inline
        GetByVector( TA&& a , IndexType i )
        {
            return GetBy( a , i , 0 );
        }





        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            expression_traits< TA >::is_single_index ,
            typename expression_traits< TA >::result_type&
        >::type
        constexpr _ehm_inline
        RefBy( TA&& a , IndexType i )
        {
            return a.Ref( i );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            expression_traits< TA >::is_single_index == false ,
            typename expression_traits< TA >::result_type&
        >::type
        constexpr _ehm_inline
        RefBy( TA&& a , IndexType i )
        {
            ERROR( "single-index access is not valid for this expression" );
            return a.Ref( i );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            expression_traits< TA >::is_single_index ,
            typename expression_traits< TA >::result_type&
        >::type
        constexpr _ehm_inline
        RefBy( TA&& a , IndexType x , IndexType y )
        {
            return a.Ref( y + x * expression_traits< TA >::rows );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            expression_traits< TA >::is_single_index == false ,
            typename expression_traits< TA >::result_type&
        >::type
        constexpr _ehm_inline
        RefBy( TA&& a , IndexType x , IndexType y )
        {
            return a.Ref( x , y );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            is_column_vector< TA >::value ,
            typename expression_traits< TA >::result_type
        >::type
        constexpr _ehm_inline
        RefByVector( TA&& a , IndexType i )
        {
            return RefBy( a , 0 , i );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            is_row_vector< TA >::value ,
            typename expression_traits< TA >::result_type
        >::type
        constexpr _ehm_inline
        RefByVector( TA&& a , IndexType i )
        {
            return RefBy( a , i , 0 );
        }

    };
};  // namespace EH
