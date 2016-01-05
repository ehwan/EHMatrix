#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_Expression.h"
#include <utility>
#include <bitset>

namespace EH
{
    namespace Matrix
    {
        namespace Expression
        {
            template < typename CLS >
            struct Expression :
                Matrix_Interface< Expression< CLS > >
            {
                static_assert( Traits< CLS >::rows > 0 , "Invalid rows for expression" );
                static_assert( Traits< CLS >::cols > 0 , "Invalid cols for expression" );

                typedef Matrix_Interface< Expression< CLS > > interface_type;

                typedef typename Traits< CLS >::return_type RET;
                constexpr const static IndexType rows = Traits< CLS >::rows;
                constexpr const static IndexType cols = Traits< CLS >::cols;

                template < typename TA >
                using temp_type = typename std::conditional<
                                        AssignShouldMakeTemp< CLS , TA >::value ,
                                        Matrix< RET , rows , cols > ,
                                        typename std::add_const< typename std::add_lvalue_reference< remove_cr<TA> >::type >::type
                                    >::type;


                _ehm_inline RET& operator [] ( IndexType i )
                {
                    return static_cast< CLS& >( *this )[ i ];
                }
                _ehm_inline RET& Get( IndexType x , IndexType y )
                {
                    return static_cast< CLS& >( *this ).Get( x , y );
                }
                _ehm_inline RET operator [] ( IndexType i ) const
                {
                    return static_cast< const CLS& >( *this )[ i ];
                }
                _ehm_inline RET Get( IndexType x , IndexType y ) const
                {
                    return static_cast< const CLS& >( *this ).Get( x , y );
                }

                _ehm_inline operator       CLS& ()       { return static_cast<       CLS& >( *this ); }
                _ehm_inline operator const CLS& () const { return static_cast< const CLS& >( *this ); }


                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& m ) const
                {
                    return static_cast< const CLS& >( *this ).has_same_root( m );
                }

                Expression(){}
                //template < typename THIS = CLS ,
                           //typename = typename std::enable_if< std::is_same< THIS , Matrix< RET , rows , cols > >::value >
                           //>
                //Expression( std::initializer_list< RET > lst )
                //{
                    //EH::LOG::LOG( " Expression Init list" );
                    //Fill( lst.begin() , lst.end() );
                //}

                // sub-matrix functions
                //{
                    constexpr _ehm_inline
                    auto
                    Column( IndexType i )
                    {
                        return SubMatrixRuntimeExp< CLS , rows , 1 >( *this , i , 0 );
                    }
                    constexpr _ehm_inline
                    auto
                    Column( IndexType i ) const
                    {
                        return SubMatrixRuntimeExp< typename std::add_const< CLS >::type , rows , 1 >( *this , i , 0 );
                    }
                    constexpr _ehm_inline
                    auto
                    Row( IndexType i )
                    {
                        return SubMatrixRuntimeExp< CLS , 1 , cols >( *this , 0 , i );
                    }
                    constexpr _ehm_inline
                    auto
                    Row( IndexType i ) const
                    {
                        return SubMatrixRuntimeExp< typename std::add_const< CLS >::type , 1 , cols >( *this , 0 , i );
                    }
                    template < IndexType M ,IndexType N >
                    constexpr _ehm_inline
                    auto
                    SubMatrix( IndexType x , IndexType y )
                    {
                        return SubMatrixRuntimeExp< CLS , M , N >( *this , x , y );
                    }
                    template < IndexType M ,IndexType N >
                    constexpr _ehm_inline
                    auto
                    SubMatrix( IndexType x , IndexType y ) const
                    {
                        return SubMatrixRuntimeExp< typename std::add_const< CLS >::type , M , N >( *this , x , y );
                    }
                //}

                template < typename TO >
                typename std::enable_if< std::is_same< RET , TO >::value , const CLS& >::type
                _ehm_inline
                Convert() const
                {
                    return *this;
                }
                template < typename TO ,
                           typename = typename std::enable_if< std::is_same< RET , TO >::value == false >::type
                         >
                auto
                _ehm_inline
                Convert() const
                {
                    return UnaryExp< CLS , TO , 0 >(
                            *this ,
                            []( const auto x )
                            {
                                return static_cast< TO >( x );
                            }
                        );
                }

                _ehm_inline auto Negative() const
                {
                    return UnaryExp< CLS , void , 0 >(
                            *this ,
                            []( const auto a )
                            {
                                return -a;
                            }
                        );
                }
                _ehm_inline auto operator - () const
                {
                    return Negative();
                }

                _ehm_inline auto Transpose()
                {
                    return TransposeTo< CLS >( *this );
                }
                _ehm_inline auto Transpose() const
                {
                    return TransposeTo< typename std::add_const< CLS >::type >( *this );
                }

                _ehm_inline auto Diagonal()
                {
                    return DiagonalExp< CLS >( *this );
                }
                _ehm_inline auto Diagonal() const
                {
                    return DiagonalExp< typename std::add_const< CLS >::type >( *this );
                }

                template < typename TA >
                auto
                _ehm_inline
                Shuffle( const expression_major_type< TA , IndexType , 0 , 0 >& idx )
                {
                    return ShuffleExp< CLS , TA >( *this , idx );
                }
                template < typename TA >
                auto
                _ehm_inline
                Shuffle( const expression_major_type< TA , IndexType , 0 , 0 >& idx ) const
                {
                    return ShuffleExp< typename std::add_const< CLS >::type , TA >( *this , idx );
                }


                // the template of expression;
                // use this functions at derived class; eg. Matrix.
                // T can be scalar or Expression
                //
                // fill with given expression;
                // DOES NOT MAKE TEMP OBJECT
                // -> be sure, it may accidently overwrite on the memory ( if later one is part of former one )
                //
                // has versions of gettable & non-gettable for performance;


                template < typename FUNC , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable >::type
                _ehm_inline
                Foreach( FUNC&& func )
                {
                    for( IndexType x=0; x<cols; ++x )
                    {
                        for( IndexType y=0; y<rows; ++y )
                        {
                            func( GetByRef( *this , x , y ) );
                        }
                    }
                }
                template < typename FUNC , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable == false >::type
                _ehm_inline
                Foreach( FUNC&& func )
                {
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        func( GetByRef( *this , i ) );
                    }
                }
                template < typename T0 , typename FUNC , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable >::type
                _ehm_inline
                Foreach( FUNC&& func , T0&& s )
                {
                    for( IndexType x=0; x<cols; ++x )
                    {
                        for( IndexType y=0; y<rows; ++y )
                        {
                            func( GetByRef( *this , x , y ) , GetBy( s , x , y ) );
                        }
                    }
                }
                template < typename T0 , typename FUNC , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable == false >::type
                _ehm_inline
                Foreach( FUNC&& func , T0&& s )
                {
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        func( GetByRef( *this , i ) , GetBy( s , i ) );
                    }
                }
                template < typename IterType , typename FUNC , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable >::type
                Foreach( FUNC&& func , IterType&& begin , IterType&& end )
                {
                    assert( std::distance( begin , end ) == rows * cols );
                    for( IndexType x=0; x<cols; ++x )
                    {
                        for( IndexType y=0; y<rows; ++y )
                        {
                            func( GetByRef( *this , x , y ) , *(begin++) );
                        }
                    }
                }
                template < typename IterType , typename FUNC , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable == false >::type
                Foreach( FUNC&& func , IterType&& begin , IterType&& end )
                {
                    assert( std::distance( begin , end ) == rows * cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        func( GetByRef( *this , i ) , *(begin++) );
                    }
                }


                template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT ,
                           typename T0 , typename ... Ts , typename FUNC >
                typename std::enable_if< _LEFT!=0 >::type
                _ehm_inline
                AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
                {
                    constexpr const IndexType _SIZE = Traits< T0 >::rows * Traits< T0 >::cols;
                    static_assert( _SIZE <= _LEFT , "Invalid expression size for aggressive assign" );
                    AggressiveForeach< OX , OY , M0 , N0 , _LEFT - _SIZE >(
                        std::template forward< FUNC >( func ) ,
                        std::template forward< Ts >( args )...
                    );
                }
                template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                           typename T0 , typename ... Ts , typename FUNC >
                typename std::enable_if<
                    _LEFT == 0 &&
                    Traits< T0 >::rows == M0 &&
                    Traits< T0 >::cols == N0 &&

                    ( OY != 0 || M0 != rows ||
                    Traits< CLS >::is_gettable || Traits< T0 >::is_gettable )
                >::type
                _ehm_inline
                AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
                {
                    for( IndexType x=0; x<N0; ++x )
                    {
                        for( IndexType y=0; y<M0; ++y )
                        {
                            func( GetByRef( *this , x + OX , y + OY ) , GetBy( arg0 , x , y ) );
                        }
                    }
                }
                template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                           typename T0 , typename ... Ts , typename FUNC >
                typename std::enable_if<
                    _LEFT == 0 &&
                    Traits< T0 >::rows == M0 &&
                    Traits< T0 >::cols == N0 &&

                    OY == 0 && M0 == rows &&
                    Traits< CLS >::is_gettable == false && Traits < T0 >::is_gettable == false
                >::type
                _ehm_inline
                AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
                {
                    for( IndexType i=0; i<M0*N0; ++i )
                    {
                        func( GetByRef( *this , i + OX*rows ) , GetBy( arg0 , i ) );
                    }
                }
                template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                           typename T0 , typename ... Ts , typename FUNC >
                typename std::enable_if<
                    _LEFT == 0 &&
                    Traits< T0 >::rows == M0 &&
                    Traits< T0 >::cols != N0
                >::type
                _ehm_inline
                AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
                {
                    constexpr const IndexType _COLS = Traits< T0 >::cols;
                    AggressiveForeach< OX , OY , M0 , _COLS      >(
                            std::template forward< FUNC >( func ) ,
                            std::template forward< T0 >( arg0 )
                    );
                    AggressiveForeach< OX + _COLS , OY , M0 , N0 - _COLS >(
                            std::template forward< FUNC >( func ) ,
                            std::template forward< Ts >( args )...
                    );
                }
                template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                           typename T0 , typename ... Ts , typename FUNC >
                typename std::enable_if<
                    _LEFT == 0 &&
                    Traits< T0 >::rows != M0 &&
                    Traits< T0 >::cols == N0
                >::type
                _ehm_inline
                AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
                {
                    constexpr const IndexType _ROWS = Traits< T0 >::rows;
                    AggressiveForeach< OX , OY , _ROWS , N0 >(
                            std::template forward< FUNC >( func ) ,
                            std::template forward< T0 >( arg0 )
                    );
                    AggressiveForeach< OX , OY + _ROWS , M0 - _ROWS , N0 >(
                            std::template forward< FUNC >( func ) ,
                            std::template forward< Ts >( args )...
                    );
                }
                template < IndexType OX = 0 , IndexType OY = 0 , IndexType M0 = rows , IndexType N0 = cols , IndexType _LEFT = 0 ,
                           typename T0 , typename ... Ts , typename FUNC >
                typename std::enable_if<
                    _LEFT == 0 &&
                    Traits< T0 >::rows != M0 &&
                    Traits< T0 >::cols != N0
                >::type
                _ehm_inline
                AggressiveForeach( FUNC&& func , T0&& arg0 , Ts&& ... args )
                {
                    constexpr const IndexType _COLS = Traits< T0 >::cols;
                    constexpr const IndexType _ROWS = Traits< T0 >::rows;
                    AggressiveForeach< OX , OY , _ROWS , _COLS >(
                            std::template forward< FUNC >( func ) ,
                            std::template forward< T0 >( arg0 )
                    );
                    AggressiveForeach< OX , OY + _ROWS , M0 - _ROWS , _COLS >(
                            std::template forward< FUNC >( func ) ,
                            std::template forward< Ts >( args )...
                    );
                    AggressiveForeach< OX + _COLS , OY , M0 , N0 - _COLS , _COLS * ( M0 - _ROWS ) >(
                            std::template forward< FUNC >( func ) ,
                            std::template forward< Ts >( args )...
                    );
                }


                _ehm_inline
                void
                Fill( const RET a )
                {
                    Foreach(
                        [ a ]( auto& x )
                        {
                             x = a;
                        }
                    );
                }
                template < typename IterType >
                _ehm_inline
                void
                Fill( IterType&& begin , IterType&& end )
                {
                    Foreach(
                        []( auto& a , const auto b )
                        {
                             a = b;
                        } ,
                        std::template forward< IterType >( begin ) , std::template forward< IterType >( end )
                    );
                }
                template < typename ... Ts >
                _ehm_inline
                void
                FillAggressive( Ts&& ... args )
                {
                    AggressiveForeach(
                            []( auto& a , const auto b )
                            {
                                a = b;
                            } ,
                            std::template forward< Ts >( args )...
                    );
                }
                template < typename ... Ts >
                _ehm_inline
                void
                MultiplyAggressive( Ts&& ... args )
                {
                    AggressiveForeach(
                            []( auto& a , const auto b )
                            {
                                a *= b;
                            } ,
                            std::template forward< Ts >( args )...
                    );
                }
                template < typename ... Ts >
                _ehm_inline
                void
                DivideAggressive( Ts&& ... args )
                {
                    AggressiveForeach(
                            []( auto& a , const auto b )
                            {
                                a /= b;
                            } ,
                            std::template forward< Ts >( args )...
                    );
                }
                template < typename ... Ts >
                _ehm_inline
                void
                PlusAggressive( Ts&& ... args )
                {
                    AggressiveForeach(
                            []( auto& a , const auto b )
                            {
                                a += b;
                            } ,
                            std::template forward< Ts >( args )...
                    );
                }
                template < typename ... Ts >
                _ehm_inline
                void
                MinusAggressive( Ts&& ... args )
                {
                    AggressiveForeach(
                            []( auto& a , const auto b )
                            {
                                a -= b;
                            } ,
                            std::template forward< Ts >( args )...
                    );
                }

                template < typename T >
                _ehm_inline
                void
                Multiply( T&& s )
                {
                    Foreach(
                        [ s ]( auto& a )
                        {
                            a *= s;
                        }
                    );
                }
                template < typename T >
                _ehm_inline
                void
                Divide( T&& s )
                {
                    Foreach(
                        [ s ]( auto& a )
                        {
                            a /= s;
                        }
                    );
                }
                template < typename T >
                _ehm_inline
                void
                Plus( T&& s )
                {
                    Foreach(
                        [ s ]( auto& a )
                        {
                            a += s;
                        }
                    );
                }
                template < typename T >
                _ehm_inline
                void
                Minus( T&& s )
                {
                    Foreach(
                        [ s ]( auto& a )
                        {
                            a -= s;
                        }
                    );
                }

                template < typename IterType >
                _ehm_inline
                void
                Multiply( IterType&& begin , IterType&& end )
                {
                    Foreach(
                        []( auto& a , const auto b )
                        {
                            a *= b;
                        } ,
                        std::template forward< IterType >( begin ) , std::template forward< IterType >( end )
                    );
                }
                template < typename IterType >
                _ehm_inline
                void
                Divide( IterType&& begin , IterType&& end )
                {
                    Foreach(
                        []( auto& a , const auto b )
                        {
                            a /= b;
                        } ,
                        std::template forward< IterType >( begin ) , std::template forward< IterType >( end )
                    );
                }
                template < typename IterType >
                _ehm_inline
                void
                Plus( IterType&& begin , IterType&& end )
                {
                    Foreach(
                        []( auto& a , const auto b )
                        {
                            a += b;
                        } ,
                        std::template forward< IterType >( begin ) , std::template forward< IterType >( end )
                    );
                }
                template < typename IterType >
                _ehm_inline
                void
                Minus( IterType&& begin , IterType&& end )
                {
                    Foreach(
                        []( auto& a , const auto b )
                        {
                            a -= b;
                        } ,
                        std::template forward< IterType >( begin ) , std::template forward< IterType >( end )
                    );
                }


                template < typename SFINE = RET >
                typename std::enable_if<
                    std::is_same< bool , SFINE >::value &&
                    Traits< CLS >::is_gettable ,
                    bool
                >::type
                all() const
                {
                    bool t = true;
                    for( IndexType x=0; x<cols; ++x )
                    {
                        for( IndexType y=0; y<rows; ++y )
                        {
                            t &= GetBy( *this , x , y );
                        }
                    }
                    return t;
                }
                template < typename SFINE = RET >
                typename std::enable_if<
                    std::is_same< bool , SFINE >::value &&
                    Traits< CLS >::is_gettable == false ,
                    bool
                >::type
                all() const
                {
                    bool t = true;
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        t &= GetBy( *this , i );
                    }
                    return t;
                }
                template < typename SFINE = RET >
                typename std::enable_if<
                    std::is_same< bool , SFINE >::value &&
                    Traits< CLS >::is_gettable ,
                    bool
                >::type
                none() const
                {
                    bool t = false;
                    for( IndexType x=0; x<cols; ++x )
                    {
                        for( IndexType y=0; y<rows; ++y )
                        {
                            t |= GetBy( *this , x , y );
                        }
                    }
                    return t == false;
                }
                template < typename SFINE = RET >
                typename std::enable_if<
                    std::is_same< bool , SFINE >::value &&
                    Traits< CLS >::is_gettable == false ,
                    bool
                >::type
                none() const
                {
                    bool t = false;
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        t |= GetBy( *this , i );
                    }
                    return t == false;
                }
                template < typename SFINE = RET >
                typename std::enable_if<
                    std::is_same< bool , SFINE >::value &&
                    Traits< CLS >::is_gettable ,
                    bool
                >::type
                any() const
                {
                    bool t = false;
                    for( IndexType x=0; x<cols; ++x )
                    {
                        for( IndexType y=0; y<rows; ++y )
                        {
                            t |= GetBy( *this , x , y );
                        }
                    }
                    return t;
                }
                template < typename SFINE = RET >
                typename std::enable_if<
                    std::is_same< bool , SFINE >::value &&
                    Traits< CLS >::is_gettable == false ,
                    bool
                >::type
                any() const
                {
                    bool t = false;
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        t |= GetBy( *this , i );
                    }
                    return t;
                }

                template < typename SFINE = RET >
                typename std::enable_if<
                    std::is_same< bool , SFINE >::value &&
                    Traits< CLS >::is_gettable ,
                    std::bitset< rows*cols >
                >::type
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
                template < typename SFINE = RET >
                typename std::enable_if<
                    std::is_same< bool , SFINE >::value &&
                    Traits< CLS >::is_gettable == false ,
                    std::bitset< rows*cols >
                >::type
                bitset() const
                {
                    std::bitset< rows*cols > ret;
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        ret.set( i , GetBy( *this , i ) );
                    }
                    return ret;
                }

                template < typename SFINE = RET >
                _ehm_inline
                operator typename std::enable_if<
                    std::is_same< SFINE , bool >::value ,
                    std::bitset< rows*cols >
                >::type
                () const
                {
                    return bitset();
                }
                template < typename SFINE = RET >
                _ehm_inline
                operator typename std::enable_if<
                    std::is_same< SFINE , bool >::value ,
                    bool
                >::type
                () const
                {
                    return all();
                }

                // simply forwarding assign operator from interface
                EXPRESSION_ASSIGN_OPERATOR( interface_type )
            };  // struct Expression
        };  // namespace Expression
    };  // namespace Matrix
};  // namespace EH
