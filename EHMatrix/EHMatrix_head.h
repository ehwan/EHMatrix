#pragma once

#include "EHMatrix_Global.h"
//#include "EHMatrix_iterator.h"
//#include "EHMatrix_container.h"
#include "EHMatrix_Expression.h"
//#include "EHMatrix_bitset.h"
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
                static_cast< THIS& >( *this ).template Multiply< T >( a );
            }
            _ehm_inline void operator /= ( const T a )
            {
                static_cast< THIS& >( *this ).template Divide< T >( a );
            }

            template < typename CLS >
            void operator *= ( const Expression::expression_size_type< CLS , N , N >& m )
            {
                static_cast< THIS& >( *this ).Fill_Safe( static_cast< const THIS& >(*this) * m );
            }

            template < typename CLS >
            _ehm_inline void operator += ( const Expression::expression_size_type< CLS , M , N >& exp )
            {
                static_cast< THIS& >( *this ).Plus_Safe( exp );
            }
            template < typename CLS >
            _ehm_inline void operator -= ( const Expression::expression_size_type< CLS , M , N >& exp )
            {
                static_cast< THIS& >( *this ).Minus_Safe( exp );
            }
            template < typename LST_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< LST_TYPE , T >::value >::type
            operator += ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Plus( lst.begin() , lst.end() );
            }
            template < typename LST_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< LST_TYPE , T >::value >::type
            operator -= ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Minus( lst.begin() , lst.end() );
            }

            template < typename EXP_CLS >
            _ehm_inline void operator = ( const Expression::Expression< EXP_CLS >& m )
            {
                static_cast< THIS& >( *this ).Fill_Safe( m );
            }
            template < typename LST_TYPE >
            _ehm_inline
            typename std::enable_if< std::is_convertible< LST_TYPE , T >::value >::type
            operator = ( std::initializer_list< LST_TYPE > lst )
            {
                static_cast< THIS& >( *this ).Fill( lst.begin() , lst.end() );
            }

            void Log() const
            {
                using Expression::GetBy;
                std::cout << "Matrix Log" << '\n';
                for( IndexType m=0; m<N; ++m )
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
                return Matrix_aliased_container< T , M , N >::s[ y ];
            }
            _ehm_inline T operator [] ( IndexType i ) const
            {
                return Matrix_aliased_container< T , M , N >::s[ i ];
            }
            _ehm_inline T Get( IndexType x , IndexType y ) const
            {
                return Matrix_aliased_container< T , M , N >::s[ y ];
            }

            // default constructor; does nothing
            Matrix(){}

            // copy constructor
            template < typename CLS >
            Matrix( const Expression::expression_size_type< CLS , M , N >& exp )
            {
                parent::template Fill<CLS>( exp );
            }

            // copy from iterator
            template < typename IterType >
            explicit Matrix( IterType begin , IterType end )
            {
                parent::Fill( begin , end );
            }
            // copy from initializer-list
            Matrix( std::initializer_list< T > lst )
            {
                parent::Fill( lst.begin() , lst.end() );
            }

            // for square-matrix
            // vector assign
            // fill diagonal vector with parameter 'v'
            // else fill 0
            template < typename SFINE = Matrix< T , M , N > , typename CLS ,
                       IndexType M2 = Expression::Traits< SFINE >::rows ,
                       IndexType N2 = Expression::Traits< SFINE >::cols ,
                       typename = typename std::enable_if< M2==N2 >::type >
            Matrix( const Expression::expression_size_type< CLS , M , 1 >& v )
            {
                parent::template Fill< T >( 0 );
                parent::Diagonal().template Fill< CLS >( v );
            }

            // square matrix scalar assign
            // same as function above
            template < typename SFINE = Matrix< T , M , N > ,
                       IndexType M2 = Expression::Traits< SFINE >::rows ,
                       IndexType N2 = Expression::Traits< SFINE >::cols ,
                       typename = typename std::enable_if< M2==N2 >::type >
            Matrix( const T sc )
            {
                parent::template Fill< T >( 0 );
                parent::Diagonal().template Fill< T >( sc );
            }

            // vector scalar assign
            // fill with given scalar
            template < typename SFINE = Matrix< T , M , N > , typename = void ,
                       IndexType M2 = Expression::Traits< SFINE >::rows ,
                       IndexType N2 = Expression::Traits< SFINE >::cols ,
                       typename = typename std::enable_if< N2==1 >::type >
            Matrix( const T sc )
            {
                parent::template Fill< T >( sc );
            }

            EXPRESSION_ASSIGN_OPERATOR( parent )

            _ehm_inline bool has_same_root( const Matrix< T , M , N >* ptr ) const
            {
                return this == ptr;
            }
            template < typename T2 , IndexType M2 , IndexType N2 >
            _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >* ptr ) const
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
                return Plus< TA , TB >( a , b );
            }
            // general minus operator
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TB >::cols >
            auto
            operator - ( const expression_size_type< TA , M , N >& a ,
                         const expression_size_type< TB , M , N >& b )
            {
                return Plus< TA , NegativeTo< TB > >( a , b.Negative() );
            }

            // general scalar scale
            template < typename TA , typename TB ,
                       typename = typename std::enable_if<
                           std::is_arithmetic< TA >::value
                           >::type >
            auto operator * ( const TA s ,
                              const Expression< TB >& m )
            {
                return Multiply< TA , TB >( s , m );
            }
            // general scalar scale
            template < typename TA , typename TB ,
                       typename = typename std::enable_if<
                           std::is_arithmetic< TA >::value
                           >::type >
            auto operator * ( const Expression< TB >& m ,
                              const TA s )
            {
                return Multiply< TB , TA >( m , s );
            }

            // general scalar scale
            template < typename TA , typename TB ,
                       typename = typename std::enable_if<
                           std::is_arithmetic< TA >::value
                           >::type >
            auto operator / ( const Expression< TB >& m ,
                              const TA s )
            {
                return Divide< TB , TA >( m , s );
            }

            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TB >::rows >
            auto operator * ( const expression_size_type< TA , M , N+1 >& m ,
                              const expression_size_type< TB , N , 1 >& v )
            {
                auto submat = m.template SubMatrix< 0 , 0 , M , N >();
                auto col = m.template Column< N >();
                return Plus< MatMatMult< remove_cr< decltype( submat ) > , TB > , remove_cr< decltype( col ) > >
                    (
                        MatMatMult< remove_cr< decltype( submat ) > , TB >( submat , v ) ,
                        col
                    );
            }
        };  // namespace Expression

        ////specialize last component is 1
        ////offset matrix multiplying
        //template < typename T1 , typename T2 , IndexType M , IndexType N , IndexType N2 ,
            //typename = typename std::enable_if< CHECK_TYPE_VALUE(T1)&&CHECK_TYPE_VALUE(T2) > >
        //BaseMatrix< typename std::common_type< T1 , T2 >::type , M , N2 > operator * ( const BaseMatrix< T1 , M , N+1 >& m1 , const BaseMatrix< T2 , N , N2 >& m2 );


        template < typename CLS >
        _ehm_inline auto floor( const Expression::Expression< CLS >& exp )
        {
            return Expression::FloorExp< CLS >( exp );
        }
        template < typename CLS >
        _ehm_inline auto ceil( const Expression::Expression< CLS >& exp )
        {
            return Expression::CeilExp< CLS >( exp );
        }
        template < typename CLS >
        _ehm_inline auto round( const Expression::Expression< CLS >& exp )
        {
            return Expression::RoundExp< CLS >( exp );
        }



        //template < typename T >
        //T Det( const BaseMatrix< T , 2 , 2 >& m );
        //template < typename T >
        //T Det( const BaseMatrix< T , 3 , 3 >& m );
        //template < typename T , IndexType N >
        //T Det( const BaseMatrix< T , N , N >& m );

        template < typename TA , typename TB , typename TC ,
                 typename = typename std::enable_if<
                     Expression::Traits< TA >::cols == 1 &&
                     Expression::Traits< TA >::cols == Expression::Traits< TB >::cols &&
                     Expression::Traits< TC >::cols == Expression::Traits< TB >::cols &&
                     Expression::Traits< TA >::rows == Expression::Traits< TB >::rows &&
                     Expression::Traits< TC >::rows == Expression::Traits< TB >::rows
                     >::type >
        auto fma( const Expression::Expression< TA >& v1 ,
                  const Expression::Expression< TB >& v2 ,
                  const Expression::Expression< TC >& v3 )
        {
            return Expression::FMA< TA , TB , TC >( v1 , v2 , v3 );
        }
        template < typename TA , typename TB ,
                 typename = typename std::enable_if<
                     Expression::Traits< TA >::cols == 1 &&
                     Expression::Traits< TA >::cols == Expression::Traits< TB >::cols &&
                     Expression::Traits< TA >::rows == Expression::Traits< TB >::rows
                     >::type >
        auto fma( const Expression::Expression< TA >& v1 ,
                  const Expression::Expression< TB >& v2 ,
                  const Expression::ret_type< TA > c )
        {
            return Expression::FMA< TA , TB , Expression::ret_type< TA > >( v1 , v2 , c );
        }
        template < typename TA , typename TC ,
                 typename = typename std::enable_if<
                     Expression::Traits< TA >::cols == 1 &&
                     Expression::Traits< TA >::cols == Expression::Traits< TC >::cols &&
                     Expression::Traits< TA >::rows == Expression::Traits< TC >::rows
                     >::type >
        auto fma( const Expression::Expression< TA >& v1 ,
                  const Expression::ret_type< TA > b ,
                  const Expression::Expression< TC >& v3 )
        {
            return Expression::FMA< TA , Expression::ret_type< TA > , TC >( v1 , b , v3 );
        }


    };  // namespace Matrix
};  //namespace EH
