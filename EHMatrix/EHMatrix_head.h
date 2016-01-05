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
            void operator *= ( const Expression::expression_size_type< CLS , N , N >& m )
            {
                static_cast< THIS& >( *this ).FillAggressive( Expression::mat_type< THIS >( static_cast< const THIS& >( *this ) * m ) );
            }

            template < typename CLS >
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , M , N >& exp )
            {
                static_cast< THIS& >( *this ).PlusAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( exp ) );
            }
            template < typename CLS >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , M , N >& exp )
            {
                static_cast< THIS& >( *this ).MinusAggressive( typename Expression::AssignShouldMakeTemp< THIS , CLS >::type( exp ) );
            }
            template < typename LST_TYPE >
            _ehm_inline
            void
            operator += ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Plus( lst.begin() , lst.end() );
            }
            template < typename LST_TYPE >
            _ehm_inline
            void
            operator -= ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Minus( lst.begin() , lst.end() );
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

            void Log() const
            {
                using Expression::GetBy;
                std::cout << "Matrix Log" << '\n';
                for( IndexType m=0; m<M; ++m )
                {
                    std::cout << "(\t" ;
                    for( IndexType n=0; n<N; ++n )
                    {
                        std::cout << GetBy( static_cast< const THIS& >( *this ) , n , m );
                        if( n!=N-1 ){ std::cout << " , "; }
                    }
                    std::cout << "\t))\n";
                }
            }
        };
        template < typename T , IndexType M , IndexType N >
        struct Matrix :
            Matrix_aliased_container< T , M , N > ,
            Expression::Expression< Matrix< T , M , N > >
        {
            typedef Expression::Expression< Matrix< T , M , N > > parent;

            _ehm_inline T& operator [] ( IndexType i )
            {
                return Matrix_aliased_container< T , M , N >::s[ i ];
            }
            _ehm_inline T& Get( IndexType x , IndexType y )
            {
                return Matrix_aliased_container< T , M , N >::s[ x*M + y ];
            }
            _ehm_inline T operator [] ( IndexType i ) const
            {
                return Matrix_aliased_container< T , M , N >::s[ i ];
            }
            _ehm_inline T Get( IndexType x , IndexType y ) const
            {
                return Matrix_aliased_container< T , M , N >::s[ x*M + y ];
            }

            // default constructor; does nothing
            Matrix(){}

            template < typename ... Ts ,
                       typename = typename std::enable_if< Expression::AggressiveSize< Ts... >::value == M*N >::type >
            Matrix( Ts&& ... args )
            {
                parent::FillAggressive(
                    std::template forward< Ts >( args )...
                );
            }

            // copy from iterator
            template < typename IterType , typename = typename std::iterator_traits< IterType >::value_type >
            explicit Matrix( IterType begin , IterType end )
            {
                parent::Fill( begin , end );
            }
            // copy from initializer-list
            //Matrix( std::initializer_list< T > lst )
            //{
                //parent::Fill( lst.begin() , lst.end() );
            //}

            // for square-matrix
            // vector assign
            // fill diagonal vector with parameter 'v'
            // else fill 0
            template < typename CLS ,
                       IndexType M2 = M ,
                       IndexType N2 = N ,
                       typename = typename std::enable_if< M2==N2 >::type >
            Matrix( const Expression::expression_size_type< CLS , M , 1 >& v )
            {
                parent::Fill( 0 );
                parent::Diagonal().FillAggressive( v );
            }

            // square matrix scalar assign
            // same as function above
            template < IndexType M2 = M ,
                       IndexType N2 = N ,
                       typename = typename std::enable_if< M2==N2 >::type >
            Matrix( const T sc )
            {
                parent::Fill( 0 );
                parent::Diagonal().Fill( sc );
            }

            // vector scalar assign
            // fill with given scalar
            template < IndexType N2 = N ,
                       typename = typename std::enable_if< N2==1 >::type >
            Matrix( const T sc )
            {
                parent::Fill( sc );
            }

            Matrix( std::initializer_list< T > lst )
            {
                parent::Fill( lst.begin() , lst.end() );
            }

            EXPRESSION_ASSIGN_OPERATOR( parent )

            constexpr
            _ehm_inline
            bool has_same_root( const Matrix< T , M , N >& ptr ) const
            {
                return this == &ptr;
            }
            template < typename T2 , IndexType M2 , IndexType N2 >
            constexpr
            _ehm_inline
            bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
            {
                return false;
            }

            // the memory functions
            //
            _ehm_inline T* data()
            {
                return Matrix_aliased_container< T , M , N >::s;
            }
            _ehm_inline const T* data() const
            {
                return Matrix_aliased_container< T , M , N >::s;
            }
            _ehm_inline T* begin()
            {
                return Matrix_aliased_container< T , M , N >::s;
            }
            _ehm_inline T* end()
            {
                return Matrix_aliased_container< T , M , N >::s + M*N;
            }
            _ehm_inline T* begin() const
            {
                return Matrix_aliased_container< T , M , N >::s;
            }
            _ehm_inline T* end() const
            {
                return Matrix_aliased_container< T , M , N >::s + M*N;
            }
            _ehm_inline T* cbegin() const
            {
                return Matrix_aliased_container< T , M , N >::s;
            }
            _ehm_inline T* cend() const
            {
                return Matrix_aliased_container< T , M , N >::s + M*N;
            }
        };


        namespace Expression
        {
            template < typename T , IndexType M , IndexType N >
            struct Traits< Matrix< T , M , N > >
            {
                constexpr const static bool is_gettable = false;
                constexpr const static IndexType cols = N;
                constexpr const static IndexType rows = M;
                typedef T return_type;
                constexpr const static bool is_restrict = true;
                constexpr const static int operations = 0;
                constexpr const static bool catch_reference = true;

                typedef Matrix< T , M , N > root_type;

                template < typename CLS >
                struct has_same_root
                {
                    constexpr const static bool value = std::is_same< CLS , root_type >::value;
                };
            };

            //normal mat-mat mult;
            //except vector dot-product
            template < typename TA , typename TB ,
                       IndexType M  = Traits< TA >::rows ,
                       IndexType N  = Traits< TA >::cols ,
                       IndexType N2 = Traits< TB >::cols ,
                       typename = typename std::enable_if<
                           ( M == 1 && N2 == 1 )==false
                           >::type >
            auto
            operator * ( const expression_size_type< TA , M , N  >& m1 ,
                         const expression_size_type< TB , N , N2 >& m2 )
            {
                return MatMatMult< TA , TB >( m1 , m2 );
            }
            //vec-vec mult; dot product , inner product
            template < typename TA , typename TB , typename = void ,
                       IndexType M  = Traits< TA >::rows ,
                       IndexType N  = Traits< TA >::cols ,
                       IndexType N2 = Traits< TB >::cols ,
                       typename = typename std::enable_if<
                           M == 1 && N != 1 && N2 == 1
                           >::type >
            auto
            operator * ( const expression_size_type< TA , M , N  >& v1 ,
                         const expression_size_type< TB , N , N2 >& v2 )
            {
                typedef typename std::common_type< ret_type< TA > , ret_type< TB > >::type ret;

                ret sum = ret( 0 );
                for( IndexType j=0; j<N; ++j )
                {
                    sum += GetBy( v1 , 0 , j ) * GetBy( v2 , 0 , j );
                }
                return sum;
            }

            // general plus operator
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TB >::cols >
            auto
            operator + ( const expression_size_type< TA , M , N >& a ,
                         const expression_size_type< TB , M , N >& b )
            {
                return BinaryExp< TA , TB >(
                        a , b ,
                        []( const auto x , const auto y )
                        {
                            return x + y;
                        }
                    );
            }
            // general minus operator
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TB >::cols >
            auto
            operator - ( const expression_size_type< TA , M , N >& a ,
                         const expression_size_type< TB , M , N >& b )
            {
                return BinaryExp< TA , TB >(
                        a , b ,
                        []( const auto x , const auto y )
                        {
                            return x - y;
                        }
                    );
            }

            // general scalar scale
            template < typename TA , typename TB ,
                       typename = typename std::enable_if<
                           is_scalar< TA >::value
                           >::type >
            auto operator * ( const TA s ,
                              const Expression< TB >& m )
            {
                return UnaryExp< TB >(
                        m ,
                        [ s ]( const auto x )
                        {
                            return x * s;
                        }
                    );
            }
            // general scalar scale
            template < typename TA , typename TB ,
                       typename = typename std::enable_if<
                           std::is_arithmetic< TA >::value
                           >::type >
            auto operator * ( const Expression< TB >& m ,
                              const TA s )
            {
                return UnaryExp< TB >(
                        m ,
                        [ s ]( const auto x )
                        {
                            return x * s;
                        }
                    );
            }

            // general scalar scale
            template < typename TA , typename TB ,
                       typename = typename std::enable_if<
                           std::is_arithmetic< TA >::value
                           >::type >
            auto operator / ( const Expression< TB >& m ,
                              const TA s )
            {
                return UnaryExp< TB >(
                        m ,
                        [ s ]( const auto x )
                        {
                            return x / s;
                        }
                    );
            }

            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TB >::rows >
            auto operator * ( const expression_size_type< TA , M , N+1 >& m ,
                              const expression_size_type< TB , N , 1 >& v )
            {
                auto submat = m.template SubMatrix< M , N >( 0 , 0 );
                auto col = m.Column( N );

                return BinaryExp< decltype( submat ) , decltype( col ) >(
                        submat , col ,
                        []( const auto x , const auto y )
                        {
                            return x + y;
                        }
                    );
            }
        };  // namespace Expression


        template < typename CLS >
        _ehm_inline auto floor( const Expression::Expression< CLS >& exp )
        {
            return Expression::UnaryExp< CLS >(
                    exp ,
                    []( const auto x )
                    {
                        return std::floor( x );
                    }
                );
        }
        template < typename CLS >
        _ehm_inline auto ceil( const Expression::Expression< CLS >& exp )
        {
            return Expression::UnaryExp< CLS >(
                    exp ,
                    []( const auto x )
                    {
                        return std::ceil( x );
                    }
                );
        }
        template < typename CLS >
        _ehm_inline auto round( const Expression::Expression< CLS >& exp )
        {
            return Expression::UnaryExp< CLS >(
                    exp ,
                    []( const auto x )
                    {
                        return std::round( x );
                    }
                );
        }
    };  // namespace Matrix
};  //namespace EH
