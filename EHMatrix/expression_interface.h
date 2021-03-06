#pragma once

#include "Global.h"

#include "expression_traits.h"
#include "expressions.h"
#include "expressions_lambda.h"

#include <utility>
#include <iterator>
#include <type_traits>
#include <bitset>

namespace EH
{
    namespace Matrix
    {
        template < typename CLS , IndexType M , IndexType N ,
                   typename Enabled = typename std::enable_if<
                       ( M == 0 || M == expression_traits< CLS >::rows ) &&
                       ( N == 0 || N == expression_traits< CLS >::cols )
                   >::type
                 >
        using expression_size_type = Expression< CLS >;
        template < typename CLS , IndexType N ,
                   typename Enabled = typename std::enable_if<
                       is_vector< CLS >::value && ( N == 0 || vector_size< CLS >::value == N )
                   >::type
                 >
        using expression_vector_size_type = Expression< CLS >;

        template < typename CRTP >
        struct Expression
        {
            using result_type = typename expression_traits< CRTP >::result_type;

            _ehm_const bool      is_single_index = expression_traits< CRTP >::is_single_index;
            _ehm_const bool      is_restrict     = expression_traits< CRTP >::is_restrict;
            _ehm_const IndexType rows            = expression_traits< CRTP >::rows;
            _ehm_const IndexType cols            = expression_traits< CRTP >::cols;

            constexpr inline operator       CRTP& ()       { return static_cast<       CRTP& >( *this ); }
            constexpr inline operator const CRTP& () const { return static_cast< const CRTP& >( *this ); }

            // reading

            template < bool SFINE = is_single_index >
            typename std::enable_if< SFINE , result_type >::type
            constexpr inline
            Get( IndexType i ) const
            {
                return GetBy( static_cast< const CRTP& >( *this ) , i );
            }
            result_type
            constexpr inline
            Get( IndexType x , IndexType y ) const
            {
                return GetBy( static_cast< const CRTP& >( *this ) , x , y );
            }


            constexpr inline
            result_type operator [] ( IndexType i ) const
            {
                return Get( i );
            }


            inline auto
            operator - () const
            {
                return Expressions::make_unary( *this ,
                        []( auto x )
                        {
                            return -x;
                        }
                    );
            }
            auto
            inline
            Diagonal() const
            {
                return Expressions::Diagonal< const CRTP >( *this );
            }
            auto
            inline
            Transpose() const
            {
                return Expressions::Transpose< const CRTP >( *this );
            }

            template < typename TO ,
                       typename = typename std::enable_if< std::is_same< TO , result_type >::value == false >::type >
            auto
            inline
            Convert() const
            {
                return Expressions::make_unary( *this ,
                        []( auto x )
                        {
                            return static_cast< TO >( x );
                        }
                    );
            }
            template < typename TO , typename = void ,
                       typename = typename std::enable_if< std::is_same< TO , result_type >::value >::type >
            auto
            inline
            Convert() const
            {
                return *this;
            }

            template < IndexType M , IndexType N >
            auto
            constexpr inline
            SubMatrix( const IndexType X , IndexType Y ) const
            {
                return Expressions::SubMatrix< const CRTP , M , N >( *this , X , Y );
            }
            auto
            constexpr inline
            Column( const IndexType X ) const
            {
                return Expressions::SubMatrix< const CRTP , rows , 1 >( *this , X , 0 );
            }
            auto
            constexpr inline
            Row( const IndexType Y ) const
            {
                return Expressions::SubMatrix< const CRTP , 1 , cols >( *this , 0 , Y );
            }

            template < typename FUNC , typename SFINE = CRTP >
            typename std::enable_if< is_single_access< SFINE >::value >::type
            constexpr inline
            ForeachConst( FUNC&& func ) const
            {
                for( IndexType i=0; i<cols*rows; ++i )
                {
                    func( GetBy( *this , i ) );
                }
            }
            template < typename FUNC , typename SFINE = CRTP >
            typename std::enable_if< is_single_access< SFINE >::value == false >::type
            constexpr inline
            ForeachConst( FUNC&& func ) const
            {
                for( IndexType x=0; x<cols; ++x )
                {
                    for( IndexType y=0; y<rows; ++y )
                    {
                        func( GetBy( *this , x , y ) );
                    }
                }
            }

            inline
            result_type sum( result_type zero = result_type(0) ) const
            {
                ForeachConst(
                        [ &zero ]( auto x )
                        {
                            zero += x;
                        }
                    );
                return zero;
            }

            constexpr inline
            matrix_type< CRTP >  matrix() const
            {
                return Matrix< result_type , rows , cols >( *this );
            }
            template < typename TO >
            inline
            typename std::conditional< is_assign_restrict< TO , CRTP >::value ,
                                       const CRTP& ,
                                       Matrix< result_type , rows , cols >
                     >::type
            assign_temp() const
            {
                return *this;
            }


            template < typename SFINE = CRTP >
            typename std::enable_if< is_vector< SFINE >::value , result_type >::type
            LengthSquared() const
            {
                result_type sum( 0 );
                ForeachConst(
                        [ &sum ]( auto x )
                        {
                            sum += x*x;
                        }
                    );
                return sum;
            }
            template < typename SFINE = CRTP >
            typename std::enable_if< is_vector< SFINE >::value , decltype( std::sqrt( result_type(0) ) ) >::type
            inline
            Length() const
            {
                return std::sqrt( LengthSquared() );
            }

            template < typename SFINE = result_type >
            typename std::enable_if< std::is_same< bool , SFINE >::value , bool >::type
            all() const
            {
                bool ret = true;
                ForeachConst(
                        [ &ret ]( auto x )
                        {
                            ret &= x;
                        }
                    );
                return ret;
            }
            template < typename SFINE = result_type >
            typename std::enable_if< std::is_same< bool , SFINE >::value , bool >::type
            none() const
            {
                bool ret = false;
                ForeachConst(
                        [ &ret ]( auto x )
                        {
                            ret |= x;
                        }
                    );
                return ret == false;
            }
            template < typename SFINE = result_type >
            typename std::enable_if< std::is_same< bool , SFINE >::value , bool >::type
            any() const
            {
                bool ret = false;
                ForeachConst(
                        [ &ret ]( auto x )
                        {
                            ret |= x;
                        }
                    );
                return ret;
            }

            template < typename SFINE = result_type , typename = typename std::enable_if< std::is_same< bool , SFINE >::value >::type >
            inline operator bool () const { return all(); }

            template < typename SFINE = result_type >
            typename std::enable_if< std::is_same< bool , SFINE >::value &&
                                     is_single_index ,
                                   std::bitset< rows*cols > >::type
            bitset() const
            {
                std::bitset< rows*cols > ret;
                for( IndexType i=0; i<rows*cols; ++i )
                {
                    ret.set( i , GetBy( *this , i ) );
                }
                return ret;
            }
            template < typename SFINE = result_type >
            typename std::enable_if< std::is_same< bool , SFINE >::value &&
                                     is_single_index == false ,
                                   std::bitset< rows*cols > >::type
            bitset() const
            {
                std::bitset< rows*cols > ret;
                for( IndexType x=0; x<cols; ++x )
                {
                    for( IndexType y=0; y<rows; ++y )
                    {
                        ret.set( y + x*rows , GetBy( *this , x , y ) );
                    }
                }
                return ret;
            }

            template < typename SFINE = result_type , typename = typename std::enable_if< std::is_same< bool , SFINE >::value >::type >
            inline operator std::bitset< rows*cols >() const { return bitset(); }


            const CRTP& Log() const
            {
                LOG( "Matrix Log :" );
                for( IndexType y=0; y<rows; ++y )
                {
                    LOGR( "(\t" );
                    for( IndexType x=0; x<cols-1; ++x )
                    {
                        LOGR( GetBy( *this , x , y ) , " , " );
                    }
                    LOG( GetBy( *this , cols-1 , y ) , "\t)" );
                }
                return *this;
            }
        }; // struct Expression

        template < typename CRTP >
        struct WritableExpression :
            Expression< CRTP >
        {
            using parent = Expression< CRTP >;
            using typename parent::result_type;

            using parent::rows;
            using parent::cols;
            using parent::is_single_index;
            using parent::is_restrict;


            constexpr inline operator       CRTP& ()       { return static_cast<       CRTP& >( *this ); }
            constexpr inline operator const CRTP& () const { return static_cast< const CRTP& >( *this ); }


            using parent::operator[];
            constexpr inline
            result_type& operator [] ( IndexType i )
            {
                return Ref( i );
            }


            // get reference
            constexpr inline
            result_type&
            Ref( IndexType i )
            {
                return RefBy( *this , i );
            }

            constexpr inline
            result_type&
            Ref( IndexType x , IndexType y )
            {
                return RefBy( *this , x , y );
            }



            using Expression< CRTP >::Diagonal;
            auto
            inline
            Diagonal()
            {
                return Expressions::Diagonal< CRTP >( *this );
            }
            using Expression< CRTP >::Transpose;
            auto
            inline
            Transpose()
            {
                return Expressions::Transpose< CRTP >( *this );
            }


            using Expression< CRTP >::SubMatrix;
            template < IndexType M , IndexType N >
            auto
            constexpr inline
            SubMatrix( IndexType X , IndexType Y )
            {
                return Expressions::SubMatrix< CRTP , M , N >( *this , X , Y );
            }
            using Expression< CRTP >::Column;
            auto
            constexpr inline
            Column( const IndexType X )
            {
                return Expressions::SubMatrix< CRTP , rows , 1 >( *this , X , 0 );
            }
            using Expression< CRTP >::Row;
            auto
            constexpr inline
            Row( const IndexType Y )
            {
                return Expressions::SubMatrix< CRTP , 1 , cols >( *this , 0 , Y );
            }



            template < typename FUNC , typename SFINE = CRTP >
            typename std::enable_if< is_single_access< SFINE >::value >::type
            constexpr inline
            Foreach( FUNC&& func )
            {
                for( IndexType i=0; i<cols*rows; ++i )
                {
                    func( RefBy( *this , i ) );
                }
            }
            template < typename FUNC , typename SFINE = CRTP >
            typename std::enable_if< is_single_access< SFINE >::value == false >::type
            constexpr inline
            Foreach( FUNC&& func )
            {
                for( IndexType x=0; x<cols; ++x )
                {
                    for( IndexType y=0; y<rows; ++y )
                    {
                        func( RefBy( *this , x , y ) );
                    }
                }
            }
            template < typename FUNC , typename SFINE = CRTP , typename IterType >
            typename std::enable_if< is_single_access< SFINE >::value >::type
            constexpr inline
            Foreach( IterType&& begin , IterType&& end , FUNC&& func )
            {
                assert( std::distance( begin , end ) == rows * cols );
                for( IndexType i=0; i<cols*rows; ++i )
                {
                    func( RefBy( *this , i ) , *(begin++) );
                }
            }
            template < typename FUNC , typename SFINE = CRTP , typename IterType >
            typename std::enable_if< is_single_access< SFINE >::value == false >::type
            constexpr inline
            Foreach( IterType&& begin , IterType&& end , FUNC&& func )
            {
                assert( std::distance( begin , end ) == rows * cols );
                for( IndexType x=0; x<cols; ++x )
                {
                    for( IndexType y=0; y<rows; ++y )
                    {
                        func( RefBy( *this , x , y ) , *(begin++) );
                    }
                }
            }
            template < IndexType OX , IndexType OY ,
                       typename T0 , typename FUNC >
            typename std::enable_if<
                OY == 0 && expression_traits< T0 >::rows == rows &&
                is_single_access< T0 >::value
            >::type
            constexpr inline
            SubForeach( FUNC&& func , T0&& arg0 )
            {
                constexpr const IndexType off = OX * expression_traits< T0 >::rows;
                for( IndexType i=0; i<matrix_size< T0 >::value; ++i )
                {
                    func( RefBy( *this , off + i ) , GetBy( arg0 , i ) );
                }
            }
            template < IndexType OX , IndexType OY ,
                       typename T0 , typename FUNC >
            typename std::enable_if<
                ( OY != 0 || expression_traits< T0 >::rows != rows ) &&
                expression_traits< T0 >::cols == 1 &&
                is_single_access< T0 >::value
            >::type
            constexpr inline
            SubForeach( FUNC&& func , T0&& arg0 )
            {
                for( IndexType i=0; i<matrix_size< T0 >::value; ++i )
                {
                    func( RefBy( *this , OX , OY + i ) , GetBy( arg0 , i ) );
                }
            }
            template < IndexType OX , IndexType OY ,
                       typename T0 , typename FUNC >
            typename std::enable_if<
                ( OY != 0 || expression_traits< T0 >::rows != rows ) &&
                expression_traits< T0 >::cols != 1 &&
                expression_traits< T0 >::rows == 1 &&
                is_single_access< T0 >::value
            >::type
            constexpr inline
            SubForeach( FUNC&& func , T0&& arg0 )
            {
                for( IndexType i=0; i<matrix_size< T0 >::value; ++i )
                {
                    func( RefBy( *this , OX + i , OY ) , GetBy( arg0 , i ) );
                }
            }
            template < IndexType OX , IndexType OY ,
                       typename T0 , typename FUNC >
            typename std::enable_if<
                ( OY != 0 || expression_traits< T0 >::rows != rows ) &&
                expression_traits< T0 >::cols != 1 &&
                expression_traits< T0 >::rows != 1 &&
                is_single_access< T0 >::value
            >::type
            constexpr inline
            SubForeach( FUNC&& func , T0&& arg0 )
            {
                for( IndexType x=0; x<expression_traits< T0 >::cols; ++x )
                {
                    for( IndexType y=0; y<expression_traits< T0 >::rows; ++y )
                    {
                        func( RefBy( *this , OX + x , OY + y ) , GetBy( arg0 , x , y ) );
                    }
                }
            }


            template < IndexType OX , IndexType OY ,
                       typename T0 , typename FUNC >
            typename std::enable_if<
                is_single_access< T0 >::value == false
            >::type
            constexpr inline
            SubForeach( FUNC&& func , T0&& arg0 )
            {
                for( IndexType x=0; x<expression_traits< T0 >::cols; ++x )
                {
                    for( IndexType y=0; y<expression_traits< T0 >::rows; ++y )
                    {
                        func( RefBy( *this , OX + x , OY + y ) , GetBy( arg0 , x , y ) );
                    }
                }
            }


            /*
             * case single-acessable submatrix
             *  - vector
             *  - same rows
             */

            template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                       typename T0 , typename ... Ts , typename FUNC >
            typename std::enable_if<
                _LEFT == 0 &&
                expression_traits< T0 >::rows == M0 &&
                expression_traits< T0 >::cols == N0
            >::type
            constexpr inline
            AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
            {
                SubForeach< OX , OY >( std::forward< FUNC >( func ) , std::forward< T0 >( arg0 ) );
            }

            template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT ,
                       typename T0 , typename ... Ts , typename FUNC >
            typename std::enable_if< _LEFT!=0 >::type
            constexpr inline
            AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
            {
                constexpr const IndexType _SIZE = expression_traits< T0 >::rows * expression_traits< T0 >::cols;
                static_assert( _SIZE <= _LEFT , "Invalid expression size for aggressive assign" );
                AggressiveForeach< OX , OY , M0 , N0 , _LEFT - _SIZE >(
                    std::forward< FUNC >( func ) ,
                    std::forward< Ts >( args )...
                );
            }
            template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                       typename T0 , typename ... Ts , typename FUNC >
            typename std::enable_if<
                _LEFT == 0 &&
                expression_traits< T0 >::rows == M0 &&
                expression_traits< T0 >::cols != N0
            >::type
            constexpr inline
            AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
            {
                constexpr const IndexType _COLS = expression_traits< T0 >::cols;
                AggressiveForeach< OX , OY , M0 , _COLS >(
                        std::forward< FUNC >( func ) ,
                        std::forward< T0 >( arg0 )
                );
                AggressiveForeach< OX + _COLS , OY , M0 , N0 - _COLS >(
                        std::forward< FUNC >( func ) ,
                        std::forward< Ts >( args )...
                );
            }
            template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                       typename T0 , typename ... Ts , typename FUNC >
            typename std::enable_if<
                _LEFT == 0 &&
                expression_traits< T0 >::rows != M0 &&
                expression_traits< T0 >::cols == N0
            >::type
            constexpr inline
            AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
            {
                constexpr const IndexType _ROWS = expression_traits< T0 >::rows;
                AggressiveForeach< OX , OY , _ROWS , N0 >(
                        std::forward< FUNC >( func ) ,
                        std::forward< T0 >( arg0 )
                );
                AggressiveForeach< OX , OY + _ROWS , M0 - _ROWS , N0 >(
                        std::forward< FUNC >( func ) ,
                        std::forward< Ts >( args )...
                );
            }
            template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                       typename T0 , typename ... Ts , typename FUNC >
            typename std::enable_if<
                _LEFT == 0 &&
                expression_traits< T0 >::rows != M0 &&
                expression_traits< T0 >::cols != N0
            >::type
            constexpr inline
            AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
            {
                constexpr const IndexType _COLS = expression_traits< T0 >::cols;
                constexpr const IndexType _ROWS = expression_traits< T0 >::rows;
                AggressiveForeach< OX , OY , _ROWS , _COLS >(
                        std::forward< FUNC >( func ) ,
                        std::forward< T0 >( arg0 )
                );
                AggressiveForeach< OX , OY + _ROWS , M0 - _ROWS , _COLS >(
                        std::forward< FUNC >( func ) ,
                        std::forward< Ts >( args )...
                );
                AggressiveForeach< OX + _COLS , OY , M0 , N0 - _COLS , _COLS * ( M0 - _ROWS ) >(
                        std::forward< FUNC >( func ) ,
                        std::forward< Ts >( args )...
                );
            }

            constexpr inline
            void Fill( const result_type i )
            {
                Foreach(
                        [ i ]( auto& a )
                        {
                            a = i;
                        }
                    );
            }
            template < typename ... Ts >
            typename std::enable_if<
                EH::static_sequence< std::size_t , expression_traits< Ts >::rows*expression_traits< Ts >::cols ... >::sum() == rows*cols
            >::type
            inline
            FillAggressive( Ts&& ... args )
            {
                AggressiveForeach(
                        []( auto& a , auto b )
                        {
                            a = b;
                        } ,
                        std::forward< Ts >( args )...
                    );
            }

            template < typename SFINE = CRTP >
            typename std::enable_if< is_vector< SFINE >::value , result_type >::type
            Normalize()
            {
                const auto L = parent::Length();
                const auto invL = decltype( L )( 1 )/L;
                if( L != 0 )
                {
                    Foreach(
                            [ invL ]( auto& x )
                            {
                                x *= invL;
                            }
                        );
                }
                return L;
            }

            template < typename TA >
            typename std::enable_if<
                is_same_size< TA , CRTP >::value ,
                CRTP&
            >::type
            inline
            operator = ( const Expression< TA >& exp )
            {
                FillAggressive( exp.template assign_temp< CRTP >() );
                return *this;
            }

            template < typename TA >
            typename std::enable_if<
                is_square< CRTP >::value &&
                is_column_vector< TA >::value &&
                vector_size< TA >::value == rows ,
                CRTP&
            >::type
            inline
            operator = ( const Expression< TA >& exp )
            {
                Fill( 0 );
                Diagonal() = exp;
                return *this;
            }

            template < typename SFINE = CRTP >
            typename std::enable_if<
                is_square< SFINE >::value ,
                CRTP&
            >::type
            inline
            operator = ( result_type scalar )
            {
                Fill( 0 );
                Diagonal().Fill( scalar );
                return *this;
            }
            template < typename SFINE = CRTP >
            typename std::enable_if<
                is_vector< SFINE >::value ,
                CRTP&
            >::type
            inline
            operator = ( result_type scalar )
            {
                Fill( scalar );
                return *this;
            }


            template < typename TA >
            typename std::enable_if<
                is_same_size< TA , CRTP >::value ,
                CRTP&
            >::type
            inline
            operator += ( const Expression< TA >& exp )
            {
                AggressiveForeach(
                        []( auto& a , auto b )
                        {
                            a += b;
                        } ,
                        exp.template assign_temp< CRTP >()
                    );
                return *this;
            }
            template < typename TA >
            typename std::enable_if<
                is_square< CRTP >::value &&
                is_column_vector< TA >::value &&
                vector_size< TA >::value == rows ,
                CRTP&
            >::type
            inline
            operator += ( const Expression< TA >& exp )
            {
                Diagonal() += exp;
                return *this;
            }
            /*
            template < typename Container , typename = decltype( std::begin( std::declval< Container >() ) ) >
            inline CRTP&
            operator += ( Container&& container )
            {
                LOG( "container plus" );
                Foreach( std::begin( std::forward< Container >( container ) ) , std::end( std::forward< Container >( container ) ) ,
                        []( auto& a , auto b )
                        {
                            a += b;
                        }
                    );

                return *this;
            }
            */
            template < typename SFINE = CRTP >
            typename std::enable_if<
                is_square< SFINE >::value ,
                CRTP&
            >::type
            inline
            operator += ( result_type scalar )
            {
                Diagonal() += scalar;
                return *this;
            }
            template < typename SFINE = CRTP >
            typename std::enable_if<
                is_vector< SFINE >::value ,
                CRTP&
            >::type
            inline
            operator += ( result_type scalar )
            {
                Foreach(
                        [ scalar ]( auto& a )
                        {
                            a += scalar;
                        }
                    );
                return *this;
            }



            template < typename TA >
            typename std::enable_if<
                is_same_size< CRTP , TA >::value ,
                CRTP&
            >::type
            inline
            operator -= ( const Expression< TA >& exp )
            {
                AggressiveForeach(
                        []( auto& a , auto b )
                        {
                            a -= b;
                        } ,
                        exp.template assign_temp< CRTP >()
                    );
                return *this;
            }
            template < typename TA >
            typename std::enable_if<
                is_square< CRTP >::value &&
                is_column_vector< TA >::value &&
                vector_size< TA >::value == rows ,
                CRTP&
            >::type
            inline
            operator -= ( const Expression< TA >& exp )
            {
                Diagonal() -= exp;
                return *this;
            }
            template < typename SFINE = CRTP >
            typename std::enable_if<
                is_square< SFINE >::value ,
                CRTP&
            >::type
            inline
            operator -= ( result_type scalar )
            {
                Diagonal() -= scalar;
                return *this;
            }
            template < typename SFINE = CRTP >
            typename std::enable_if<
                is_vector< SFINE >::value ,
                CRTP&
            >::type
            inline
            operator -= ( result_type scalar )
            {
                Foreach(
                        [ scalar ]( auto& a )
                        {
                            a -= scalar;
                        }
                    );
                return *this;
            }

            inline CRTP& operator *= ( result_type scalar ) {
                Foreach(
                        [ scalar ]( auto& a )
                        {
                            a *= scalar;
                        }
                    );
                return *this;
            }
            template < typename TA >
            typename std::enable_if<
                is_same_vector< TA , CRTP >::value &&
                vector_size< TA >::value == vector_size< CRTP >::value ,
                CRTP&
            >::type
            inline operator *= ( const Expression< TA >& exp )
            {
                AggressiveForeach(
                        []( auto& a , auto b )
                        {
                            a *= b;
                        } ,
                        exp.template assign_temp< CRTP >()
                    );
                return *this;
            }
            template < typename TA >
            typename std::enable_if<
                is_square< TA >::value &&
                expression_traits< TA >::cols == cols &&
                is_column_vector< CRTP >::value == false ,
                CRTP&
            >::type
            inline operator *= ( const Expression< TA >& exp )
            {
                FillAggressive( Matrix< result_type , rows , cols >( (*this) * exp ) );
                return *this;
            }


            inline CRTP& operator /= ( result_type scalar )
            {
                Foreach(
                        [ scalar ]( auto& a )
                        {
                            a /= scalar;
                        }
                    );
                return *this;
            }
            template < typename TA >
            typename std::enable_if<
                is_same_vector< TA , CRTP >::value &&
                vector_size< TA >::value == vector_size< CRTP >::value ,
                CRTP&
            >::type
            inline operator /= ( const Expression< TA >& exp )
            {
                AggressiveForeach(
                        []( auto& a , auto b )
                        {
                            a /= b;
                        } ,
                        exp.template assign_temp< CRTP >()
                    );
                return *this;
            }

        }; // struct WritableExpression

    };  // namespace Matrix
};  // namespace EH
