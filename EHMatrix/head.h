# pragma once

#include "Global.h"
#include "expression_interface.h"
#include "expression_traits.h"
#include "aliased_container.h"

#include <type_traits>
#include <initializer_list>

namespace EH
{
    namespace Matrix
    {
        template < typename T , IndexType M , IndexType N >
        struct Matrix :
            aliased_container< T , M , N > ,
            //Expression< Matrix< T , M , N > >
            WritableExpression< Matrix< T , M , N > >
        {
            using parent = WritableExpression< Matrix< T , M , N > >;
            using typename parent::result_type;

            using parent::operator=;

            constexpr _ehm_inline T Get( IndexType i ) const
            {
                return aliased_container< T , M , N >::s[ i ];
            }
            constexpr _ehm_inline T Get( IndexType x , IndexType y ) const
            {
                return aliased_container< T , M , N >::s[ y + x*M ];
            }
            constexpr _ehm_inline T& Ref( IndexType i )
            {
                return aliased_container< T , M , N >::s[ i ];
            }
            constexpr _ehm_inline T& Ref( IndexType x , IndexType y )
            {
                return aliased_container< T , M , N >::s[ y + x*M ];
            }

            constexpr Matrix(){}

            template < typename ... Ts ,
                       typename = typename std::enable_if<
                           EH::static_sequence< IndexType , matrix_size< Ts >::value ... >::sum() == M*N &&
                           sizeof...( Ts ) >= 1
                        >::type >
            constexpr Matrix( Ts&& ... args )
            {
                parent::FillAggressive( std::forward< Ts >( args )... );
            }

            template < typename ... Ts , typename = void ,
                       typename = typename std::enable_if<
                           M == N &&
                           EH::static_sequence< IndexType , matrix_size< Ts >::value ... >::sum() == M &&
                           sizeof...( Ts ) >= 1
                        >::type >
            constexpr Matrix( Ts&& ... args )
            {
                parent::Fill( 0 );
                parent::Diagonal().FillAggressive( std::forward< Ts >( args )... );
            }

            template < IndexType M0 = M , IndexType N0 = N ,
                       typename = typename std::enable_if< N0 == M0 >::type >
            constexpr Matrix( const T a )
            {
                parent::Fill( 0 );
                parent::Diagonal().Fill( a );
            }
            template < IndexType M0 = M , IndexType N0 = N , typename = void ,
                       typename = typename std::enable_if< M0 == 1 || N0 == 1 >::type >
            constexpr Matrix( const T a )
            {
                parent::Fill( a );
            }

            constexpr _ehm_inline
            T* begin()
            {
                return aliased_container< T , M , N >::s;
            }
            constexpr _ehm_inline
            T* end()
            {
                return aliased_container< T , M , N >::s + M*N;
            }
            constexpr _ehm_inline
            const T* begin() const
            {
                return aliased_container< T , M , N >::s;
            }
            constexpr _ehm_inline
            const T* end() const
            {
                return aliased_container< T , M , N >::s + M*N;
            }

            constexpr _ehm_inline
            T* data()
            {
                return aliased_container< T , M , N >::s;
            }
            constexpr _ehm_inline
            const T* data() const
            {
                return aliased_container< T , M , N >::s;
            }
        };

        template < typename T , IndexType M , IndexType N >
        struct expression_traits< Matrix< T , M , N > >
        {
            using result_type = T;
            using root_type   = Matrix< T , M , N >;

            _ehm_const IndexType cols            = N;
            _ehm_const IndexType rows            = M;

            _ehm_const bool      is_single_index = true;
            _ehm_const bool      is_restrict     = true;
            _ehm_const bool      catch_reference = true;

            _ehm_const bool      is_expression   = true;

            _ehm_const int       operations      = 0;
        };


        template < typename TA , typename TB , typename = void >
        auto Cross( Vector< TA , 2 >& v1 ,
                    Vector< TB , 2 >& v2 )
        {
            return GetBy( v1 , 0 , 0 )*GetBy( v2 , 0 , 1 ) - GetBy( v1 , 0 , 1 )*GetBy( v2 , 0 , 0 );
        }
        template < typename TA , typename = typename std::enable_if< is_expression< TA >::value &&
                                                                     is_vector< TA >::value &&
                                                                     vector_size< TA >::value == 2
                                                                   >::type
                 >
        auto
        Cross( TA&& v ,
               typename expression_traits< TA >::result_type a )
        {
            return Expressions::make_unary( Expressions::Vector2Skew< TA >( std::forward< TA >( v ) ) ,
                    [ a ]( const auto x )
                    {
                        return x * a;
                    }
                );
        }
        template < typename TA , typename = typename std::enable_if< is_expression< TA >::value &&
                                                                     is_vector< TA >::value &&
                                                                     vector_size< TA >::value == 2
                                                                   >::type
                 >
        auto
        Cross( typename expression_traits< TA >::result_type a ,
               TA&& v )
        {
            return Expressions::make_unary( Expressions::Vector2Skew< TA >( std::forward< TA >( v ) ) ,
                    [ a ]( const auto x )
                    {
                        return - x * a;
                    }
                );
        }
        template < typename TA >
        typename std::enable_if<
            is_expression< TA >::value &&
            is_vector< TA >::value &&
            vector_size< TA >::value == 2 ,
            Expressions::Vector2Skew< TA >
        >::type
        Cross( TA&& v )
        {
            return Expressions::Vector2Skew< TA >( std::forward< TA >( v ) );
        }
        template < typename TA , typename TB >
        typename std::enable_if<
            is_expression< TA >::value && is_expression< TB >::value &&
            is_column_vector< TA >::value && is_column_vector< TB >::value &&
            vector_size< TA >::value == 3 && vector_size< TB >::value == 3 ,
            Matrix< typename expression_traits< TA , TB >::result_type ,
                    3 , 1 >
        >::type
        Cross( const expression_size_type< TA , 3 , 1 >& _a ,
               const expression_size_type< TB , 3 , 1 >& _b )
        {
            typename Expressions::ShouldMakeTemp< const TA , 2 >::type a( _a );
            typename Expressions::ShouldMakeTemp< const TB , 2 >::type b( _b );
            return Vector< typename expression_traits< TA , TB >::result_type , 3 >
            {
                GetBy( a , 0 , 1 )*GetBy( b , 0 , 2 ) - GetBy( a , 0 , 2 )*GetBy( b , 0 , 1 ) ,
                GetBy( a , 0 , 2 )*GetBy( b , 0 , 0 ) - GetBy( a , 0 , 0 )*GetBy( b , 0 , 2 ) ,
                GetBy( a , 0 , 0 )*GetBy( b , 0 , 1 ) - GetBy( a , 0 , 1 )*GetBy( b , 0 , 0 )
            };
        }

        template < typename TA >
        typename std::enable_if<
            is_square< TA >::value ,
            typename expression_traits< TA >::value
        >::type
        _ehm_inline
        Det( const Expression< TA >& exp )
        {
        }

        template < typename TA >
        auto Inverse( const expression_size_type< TA , 2 , 2 >& e1 )
        {
            typename Expressions::ShouldMakeTemp< const TA , 2 >::type temp( e1 );
            using result_type = typename expression_traits< TA >::result_type;
            const result_type d = result_type( 1 ) /
                ( GetBy( temp , 0 , 0 )*GetBy( temp , 1 , 1 ) - GetBy( temp , 1 , 0 )*GetBy( temp , 0 , 1 ) );

            return Matrix< result_type , 2 , 2 >
            {
                d*GetBy( temp , 1 , 1 ) , -d*GetBy( temp , 0 , 1 ) , -d*GetBy( temp , 1 , 0 ) , d*GetBy( temp , 0 , 0 )
            };
        }
        template < typename TA >
        auto Inverse( const Matrix< TA , 3 , 3 >& e1 )
        {
            const Vector< TA , 3 > c12 = Cross( e1.Column( 1 ) , e1.Column( 2 ) );
            const Vector< TA , 3 > c20 = Cross( e1.Column( 2 ) , e1.Column( 0 ) );
            const Vector< TA , 3 > c01 = Cross( e1.Column( 0 ) , e1.Column( 1 ) );
            TA det = e1.Column( 0 ).Transpose() * c12;

            TA detI = TA( 1 ) / det;
            return Matrix< TA , 3 , 3 >
            {
                detI * c12.Transpose() ,
                detI * c20.Transpose() ,
                detI * c01.Transpose()
            };
        }

        template < typename TA ,
                   typename = typename std::enable_if<
                       is_expression< TA >::value
                   >::type
                 >
        _ehm_inline auto
        floor( TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    []( auto x )
                    {
                        return std::floor( x );
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if<
                       is_expression< TA >::value
                   >::type
                 >
        _ehm_inline auto
        ceil( TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    []( auto x )
                    {
                        return std::ceil( x );
                    }
                );
        }
        template < typename TA ,
                   typename = typename std::enable_if<
                       is_expression< TA >::value
                   >::type
                 >
        _ehm_inline auto
        round( TA&& exp )
        {
            return Expressions::make_unary( std::forward< TA >( exp ) ,
                    []( auto x )
                    {
                        return std::round( x );
                    }
                );
        }
        template < typename TA >
        _ehm_inline auto
        abs( const Expression< TA >& exp )
        {
            return Expressions::make_unary( exp ,
                    []( auto x )
                    {
                        return std::abs( x );
                    }
                );
        }
        template < typename TA , typename TB >
        _ehm_inline auto
        max( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )
                    {
                        return std::max( a , b );
                    }
                );
        }
        template < typename TA , typename TB >
        _ehm_inline auto
        min( const Expression< TA >& e1 , const Expression< TB >& e2 )
        {
            return Expressions::make_binary( e1 , e2 ,
                    []( auto a , auto b )
                    {
                        return std::min( a , b );
                    }
                );
        }

        template < typename TA >
        _ehm_inline auto
        sqrt( const Expression< TA >& e1 )
        {
            return Expressions::make_unary( e1 ,
                    []( auto a )
                    {
                        return std::sqrt( a );
                    }
                );
        }

        template < typename T >
        inline auto
        clamp( T val , T min , T max )
        {
            return std::max( std::min( val , max ) );
        }
        template < typename T >
        inline auto
        constexpr cycle( T val , T min , T max )
        {
            if( val > max ){ return val - ( max - min ); }
            else if( val < min ){ return val + ( max - min ); }
            return val;
        }
        template < typename T >
        inline auto
        clamp( const Expression< T >& e1 , typename expression_traits< T >::result_type min , typename expression_traits< T >::result_type max )
        {
            return Expressions::make_unary( e1 ,
                    [ min , max ]( auto x )
                    {
                        return clamp( x , min , max );
                    }
                );
        }
        template < typename T >
        inline auto
        cycle( const Expression< T >& e1 , typename expression_traits< T >::result_type min , typename expression_traits< T >::result_type max )
        {
            return Expressions::make_unary( e1 ,
                    [ min , max ]( auto x )
                    {
                        return cycle( x , min , max );
                    }
                );
        }
    };  // namespace Matrix
};  //namespace EH
