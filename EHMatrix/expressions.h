#pragma once

#include "Global.h"
#include "expression_traits.h"

namespace EH
{
    namespace Matrix
    {
        namespace Expressions
        {

            template < typename T , int FACTOR >
            struct ShouldMakeTemp
            {
                /*
                 * for Expression that has many access to same index; so make it to calculate many times
                 * decide whether to make temp variable to cache the result
                 * expression's size & operations are counting with it; ( how many times they access? )
                 * Access_Factor would be how many times you will access same index;
                 */
                _ehm_const bool value =
                    expression_traits< T >::operations * ( FACTOR  - 1 ) > expression_traits< T >::rows*expression_traits< T >::cols;

                // if should make temp object, returns matrix type
                // else returns expression-reference type
                typedef typename std::conditional<
                            value ,
                            const Matrix< typename expression_traits< T >::result_type ,
                                          expression_traits< T >::rows ,
                                          expression_traits< T >::cols
                                        > ,
                            auto_reference< T >
                        >::type type;

                _ehm_const bool is_single_index = value ? true : expression_traits< T >::is_single_index;
                _ehm_const bool is_restrict     = value ? true : expression_traits< T >::is_restrict;

                _ehm_const int  operations      = value ? 0    : expression_traits< T >::operations;
            };

            template < typename TA >
            struct Diagonal :
                WritableExpression< Diagonal< TA > >
            {
                using parent = WritableExpression< Diagonal< TA > >;
                using parent::operator=;

                auto_reference< TA > a;

                Diagonal( auto_reference< TA > _a ) :
                    a( _a )
                {
                }


                _ehm_inline auto Get( IndexType i ) const
                {
                    return GetBy( a , i , i );
                }
                _ehm_inline auto Get( IndexType x , IndexType y ) const
                {
                    return GetBy( a , y , y );
                }
                _ehm_inline auto& Ref( IndexType i )
                {
                    return RefBy( a , i , i );
                }
                _ehm_inline auto& Ref( IndexType x , IndexType y )
                {
                    return RefBy( a , y , y );
                }
            };
            template < typename TA >
            struct Transpose :
                WritableExpression< Transpose< TA > >
            {
                using parent = WritableExpression< Transpose< TA > >;
                using parent::operator=;
                using typename WritableExpression< Transpose< TA > >::result_type;

                auto_reference< TA > a;

                Transpose( auto_reference< TA > _a ) :
                    a( _a )
                {
                }

                template < typename SFINE = TA >
                typename std::enable_if< is_row_vector< SFINE >::value , result_type >::type
                _ehm_inline
                Get( IndexType i ) const
                {
                    return GetBy( a , i , 0 );
                }
                template < typename SFINE = TA >
                typename std::enable_if< is_column_vector< SFINE >::value , result_type >::type
                _ehm_inline
                Get( IndexType i ) const
                {
                    return GetBy( a , 0 , i );
                }
                _ehm_inline result_type Get( IndexType x , IndexType y ) const
                {
                    return GetBy( a , y , x );
                }
                template < typename SFINE = TA >
                typename std::enable_if< is_row_vector< SFINE >::value , result_type& >::type
                _ehm_inline
                Ref( IndexType i )
                {
                    return RefBy( a , i , 0 );
                }
                template < typename SFINE = TA >
                typename std::enable_if< is_column_vector< SFINE >::value , result_type& >::type
                _ehm_inline
                Ref( IndexType i )
                {
                    return RefBy( a , 0 , i );
                }
                _ehm_inline result_type& Ref( IndexType x , IndexType y )
                {
                    return RefBy( a , y , x );
                }
            };

            template < typename TA , IndexType M , IndexType N >
            struct SubMatrix :
                WritableExpression< SubMatrix< TA , M , N > >
            {
                using parent = WritableExpression< SubMatrix< TA , M , N > >;
                using parent::operator=;

                auto_reference< TA > a;

                const IndexType X , Y;

                SubMatrix( auto_reference< TA >  _a , const IndexType _x , const IndexType _y ) :
                    a( _a ) , X( _x ) , Y( _y )
                {
                }

                constexpr _ehm_inline auto Get( IndexType i ) const
                {
                    return GetBy( a , i + X*expression_traits< TA >::rows );
                }
                constexpr _ehm_inline auto Get( IndexType x , IndexType y ) const
                {
                    return GetBy( a , x + X , y + Y );
                }
                constexpr _ehm_inline auto& Ref( IndexType i )
                {
                    return RefBy( a , i + X*expression_traits< TA >::rows );
                }
                constexpr _ehm_inline auto& Ref( IndexType x , IndexType y )
                {
                    return RefBy( a , x + X , y + Y );
                }
            };

            template < typename TA >
            struct Vector2Square :
                Expression< Vector2Square< TA > >
            {
                using Expression< Vector2Square< TA > >::rows;
                using Expression< Vector2Square< TA > >::cols;
                using typename Expression< Vector2Square< TA > >::result_type;

                auto_reference< TA > a;

                Vector2Square( auto_reference< TA > _a ) :
                    a( _a )
                {
                }

                //template < typename TTA >
                //Vector2Square( TTA&& _a ) :
                    //a( _a )
                //{
                //}

                constexpr _ehm_inline
                result_type
                Get( IndexType i ) const
                {
                    constexpr const IndexType pow = std::log2( rows );
                    return ( i % rows ) == 0 ? GetBy( a , i >> pow ) : 0;
                }
                constexpr _ehm_inline
                result_type
                Get( IndexType x , IndexType y ) const
                {
                    return x == y ? GetBy( a , y ) : 0;
                }
            };


            template < typename TA , typename TB >
            struct MatMatMult :
                Expression< MatMatMult< TA , TB > >
            {
                using typename Expression< MatMatMult< TA , TB > >::result_type;

                const typename ShouldMakeTemp< TA , expression_traits< TB >::cols >::type a;
                const typename ShouldMakeTemp< TB , expression_traits< TA >::rows >::type b;

                MatMatMult( auto_reference< TA > _a , auto_reference< TB > _b ) :
                    a( _a ) , b( _b )
                {
                }

                // vec-mat mult
                template < typename _TA = TA , typename _TB = TB >
                typename std::enable_if< is_row_vector< _TA >::value , result_type >::type
                Get( IndexType i ) const
                {
                    result_type sum = result_type( 0 );

                    for( IndexType j=0; j<expression_traits< TA >::cols; ++j )
                    {
                        sum += GetBy( a , j , 0 ) * GetBy( b , i , j );
                    }

                    return sum;
                }
                // mat-vec mult
                template < typename _TA = TA , typename _TB = TB >
                typename std::enable_if< is_column_vector< _TB >::value , result_type >::type
                Get( IndexType i ) const
                {
                    result_type sum = result_type( 0 );

                    for( IndexType j=0; j<expression_traits< TA >::cols; ++j )
                    {
                        sum += GetBy( a , j , i ) * GetBy( b , 0 , j );
                    }

                    return sum;
                }
                // and the others
                auto Get( IndexType x , IndexType y ) const
                {
                    result_type sum = result_type( 0 );

                    for( IndexType j=0; j<expression_traits< TA >::cols; ++j )
                    {
                        sum += GetBy( a , j , y ) * GetBy( b , x , j );
                    }

                    return sum;
                }

            };


        };  // namespace Expression


        template < typename TA >
        struct expression_traits< Expressions::Diagonal< TA > > : expression_traits< TA >
        {
            _ehm_const IndexType rows            = std::min( expression_traits< TA >::rows , expression_traits< TA >::cols );
            _ehm_const IndexType cols            = 1;
            _ehm_const bool      catch_reference = false;
            _ehm_const bool      is_single_index = true;
            _ehm_const bool      is_restrict     = false;
        };
        template < typename TA >
        struct expression_traits< Expressions::Transpose< TA > > : expression_traits< TA >
        {
            _ehm_const IndexType rows            = expression_traits< TA >::cols;
            _ehm_const IndexType cols            = expression_traits< TA >::rows;
            _ehm_const bool      catch_reference = false;
            _ehm_const bool      is_single_index = is_vector< TA >::value;
            _ehm_const bool      is_restrict     = false;
        };

        template < typename TA , IndexType M , IndexType N >
        struct expression_traits< Expressions::SubMatrix< TA , M , N > > : expression_traits< TA >
        {
            // case it's index accessable;
            //   -original is index accessable AND
            //      -original is vector type
            //      -same M( rows count ) as original
            _ehm_const bool      is_single_index = expression_traits< TA >::is_single_index &&
                                                   ( M == expression_traits< TA >::rows || is_vector< TA >::value );
            _ehm_const bool      catch_reference = false;
            _ehm_const bool      is_restrict     = false;

            _ehm_const IndexType rows            = M;
            _ehm_const IndexType cols            = N;

            //_ehm_const int       operations      = expression_traits< TA >::operations + 1;
        };

        template < typename TA >
        struct expression_traits< Expressions::Vector2Square< TA > > : expression_traits< TA >
        {
            _ehm_const bool catch_reference = false;
            _ehm_const bool is_restrict     = false;
            _ehm_const bool is_single_index = expression_traits< TA >::is_single_index &&
                                              ( expression_traits< TA >::rows & ( expression_traits< TA >::rows + 1 ) ) == 0;

            _ehm_const IndexType rows = vector_size< TA >::value;
            _ehm_const IndexType cols = rows;

            _ehm_const int operations = expression_traits< TA >::operations + 1;
        };



        template < typename TA , typename TB >
        struct expression_traits< Expressions::MatMatMult< TA , TB > > : expression_traits< TA , TB >
        {
            using tempA = Expressions::ShouldMakeTemp< TA , expression_traits< TB >::cols >;
            using tempB = Expressions::ShouldMakeTemp< TB , expression_traits< TA >::rows >;

            _ehm_const IndexType cols = expression_traits< TB >::cols;
            _ehm_const IndexType rows = expression_traits< TA >::rows;

            // if return is vector ( either row or col ) can be index-accessed
            _ehm_const bool is_single_index = ( cols == 1 || rows == 1 );
            _ehm_const bool is_restrict     = tempA::value && tempB::value;
            _ehm_const bool catch_reference = tempA::value || tempB::value;


            _ehm_const int  operations      = ( tempA::operations + tempB::operations ) * expression_traits< TA >::cols;
        };


    };  // namespace Matrix
};  // namespace EH
