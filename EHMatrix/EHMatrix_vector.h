#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_head.h"

#include <initializer_list>
#include <cmath>

namespace EH
{
    namespace Matrix
    {
        template < typename THIS , IndexType N >
        struct Matrix_Interface< THIS , N , 1 >
        {
            constexpr const static IndexType M = Expression::Traits< THIS >::rows;
            typedef Expression::ret_type< THIS > T;

            _ehm_inline void operator *= ( const T a )
            {
                static_cast< THIS& >( *this ).Multiply( a );
            }
            _ehm_inline void operator /= ( const T a )
            {
                static_cast< THIS& >( *this ).Divide( a );
            }

            template < typename CLS >
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , M , 1 >& exp )
            {
                static_cast< THIS& >( *this ).PlusAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( exp ) );
            }
            template < typename CLS >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , M , 1 >& exp )
            {
                static_cast< THIS& >( *this ).MinusAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( exp ) );
            }
            _ehm_inline
            void
            operator += ( std::initializer_list< T > lst )
            {
                static_cast< THIS& >( *this ).Plus( lst.begin() , lst.end() );
            }
            _ehm_inline
            void
            operator -= ( std::initializer_list< T > lst )
            {
                static_cast< THIS& >( *this ).Minus( lst.begin() , lst.end() );
            }

            _ehm_inline void operator += ( const T a )
            {
                static_cast< THIS& >( *this ).Plus( a );
            }
            _ehm_inline void operator -= ( const T a )
            {
                static_cast< THIS& >( *this ).Minus( a );
            }

            template < typename CLS >
            _ehm_inline void operator *= ( const Expression::expression_size_type< CLS , M , 1 >& v )
            {
                static_cast< THIS& >( *this ).MultiplyAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( v ) );
            }
            template < typename CLS >
            _ehm_inline void operator /= ( const Expression::expression_size_type< CLS , M , 1 >& v )
            {
                static_cast< THIS& >( *this ).DivideAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::tpye( v ) );
            }
            _ehm_inline
            void
            operator *= ( std::initializer_list< T > lst )
            {
                static_cast< THIS& >( *this ).Multiply( lst.begin() , lst.end() );
            }
            _ehm_inline
            void
            operator /= ( std::initializer_list< T > lst )
            {
                static_cast< THIS& >( *this ).Divide( lst.begin() , lst.end() );
            }

            template < typename EXP_CLS >
            _ehm_inline void operator = ( const Expression::Expression< EXP_CLS >& m )
            {
                static_cast< THIS& >( *this ).FillAggressive( typename Expression::AssignShouldMakeTemp< THIS , EXP_CLS >::type( m ) );
            }
            _ehm_inline
            void
            operator = ( std::initializer_list< T > lst )
            {
                static_cast< THIS& >( *this ).Fill( lst.begin() , lst.end() );
            }
            _ehm_inline
            void
            operator = ( const T a )
            {
                static_cast< THIS& >( *this ).Fill( a );
            }

            template < typename SFINE = THIS >
            typename std::enable_if< Expression::Traits< SFINE >::operations <= 1 , T >::type
            LengthSquared() const
            {
                using Expression::GetBy;
                T sum = T( 0 );
                for( IndexType i=0; i<M; ++i )
                {
                    sum += GetBy( static_cast< const THIS& >( *this ) , 0 , i ) * GetBy( static_cast< const THIS& >( *this ) , 0 , i );
                }
                return sum;
            }

            template < typename SFINE = THIS >
            typename std::enable_if< (Expression::Traits< SFINE >::operations > 1) , T >::type
            LengthSquared() const
            {
                using Expression::GetBy;
                T sum = T( 0 );
                for( IndexType i=0; i<M; ++i )
                {
                    const T t = GetBy( static_cast< const THIS& >( *this ) , 0 , i );
                    sum += t * t;
                }
                return sum;
            }
            auto Length() const
            {
                return std::sqrt( LengthSquared() );
            }
            auto Normalize()
            {
                using Expression::GetByRef;
                const auto l = Length();
                const auto invl = decltype( l )( 1 ) / l;
                for( IndexType i=0; i<M; ++i )
                {
                    GetByRef( static_cast< THIS& >( *this ) , 0 , i ) *= invl;
                }
                return l;
            }

            void Log() const
            {
                using Expression::GetBy;
                std::cout << "Vector Log\n";
                for( IndexType m=0; m<M; ++m )
                {
                    std::cout << "(\t";
                    std::cout << GetBy( static_cast< const THIS& >( *this ) , 0 , m );
                    std::cout << "\t)\n";
                }
            }
        };

        namespace Expression
        {
            // zero extra typename=void is for square-scalar operations
            template < typename TA , typename = void >
            auto
            operator + ( const expression_size_type< TA , 0 , 1 >& v ,
                         const ret_type< TA > s )
            {
                return UnaryExp< TA >(
                        v ,
                        [ s ]( const auto x )
                        {
                            return x + s;
                        }
                    );
            }
            template < typename TA , typename = void >
            auto
            operator - ( const expression_size_type< TA , 0 , 1 >& v ,
                         const ret_type< TA > s )
            {
                return UnaryExp< TA >(
                        v ,
                        [ s ]( const auto x )
                        {
                            return x - s;
                        }
                    );
            }
            template < typename TA , typename = void >
            auto
            operator + ( const ret_type< TA > s ,
                         const expression_size_type< TA , 0 , 1 >& v )
            {
                return UnaryExp< TA >(
                        v ,
                        [ s ]( const auto x )
                        {
                            return x + s;
                        }
                    );
            }
            template < typename TA , typename = void >
            auto
            operator - ( const ret_type< TA > s ,
                         const expression_size_type< TA , 0 , 1 >& v )
            {
                return UnaryExp< TA >(
                        v ,
                        [ s ]( const auto x )
                        {
                            return x - s;
                        }
                    );
            }

            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator * ( const expression_size_type< TA , M , 1 >& v1 ,
                         const expression_size_type< TB , M , 1 >& v2 )
            {
                return BinaryExp< TA , TB >(
                        v1 , v2 ,
                        []( const auto x , const auto y )
                        {
                            return x * y;
                        }
                        );
            }
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows >
            auto
            operator / ( const expression_size_type< TA , M , 1 >& v1 ,
                         const expression_size_type< TB , M , 1 >& v2 )
            {
                return BinaryExp< TA , TB >(
                        v1 , v2 ,
                        []( const auto x , const auto y )
                        {
                            return x / y;
                        }
                    );
            }
        };  // namespace Expression


        template < typename TA , typename TB , typename = void >
        typename std::common_type< Expression::ret_type< TA > , Expression::ret_type< TB > >::type
        Cross( const Expression::expression_size_type< TA , 2 , 1 >& v1 ,
               const Expression::expression_size_type< TB , 2 , 1 >& v2 )
        {
            using Expression::GetBy;
            return GetBy( v1 , 0 , 0 )*GetBy( v2 , 0 , 1 ) - GetBy( v1 , 0 , 1 )*GetBy( v2 , 0 , 0 );
        }
        template < typename TA >
        auto
        Cross( const Expression::expression_size_type< TA , 2 , 1 >& v ,
               const Expression::ret_type< TA > a )
        {
            return make_unary_expression(
                    Expression::OuterProduct< 2 , TA >( v ) ,
                    [ a ]( const auto x )
                    {
                        return x * a;
                    }
                );
        }
        template < typename TA >
        auto
        Cross( const Expression::ret_type< TA > a ,
               const Expression::expression_size_type< TA , 2 , 1 >& v )
        {
            return Expression::UnaryExp< Expression::OuterProduct< 2 , TA > >(
                    Expression::OuterProduct< 2 , TA >( v ) ,
                    [ a ]( const auto x )
                    {
                        return - x * a;
                    }
                );
        }
        template < typename TA >
        auto
        Cross( const Expression::expression_size_type< TA , 2 , 1 >& v )
        {
            return Expression::OuterProduct< 2 , TA >( v );
        }
        template < typename TA , typename TB >
        auto
        Cross( const Expression::expression_size_type< TA , 3 , 1 >& _a ,
               const Expression::expression_size_type< TB , 3 , 1 >& _b )
        {
            using Expression::GetBy;
            typename Expression::ShouldMakeTemp< TA , 2 >::type a( _a );
            typename Expression::ShouldMakeTemp< TB , 2 >::type b( _b );
            return Vector< typename std::common_type< Expression::ret_type< TA > , Expression::ret_type< TB > >::type , 3 >
            {
                {
                    GetBy( a , 0 , 1 )*GetBy( b , 0 , 2 ) - GetBy( a , 0 , 2 )*GetBy( b , 0 , 1 ) ,
                    GetBy( a , 0 , 2 )*GetBy( b , 0 , 0 ) - GetBy( a , 0 , 0 )*GetBy( b , 0 , 2 ) ,
                    GetBy( a , 0 , 0 )*GetBy( b , 0 , 1 ) - GetBy( a , 0 , 1 )*GetBy( b , 0 , 0 )
                }
            };
        }
    };  // namespace Matrix
};  // namespace EH
