#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_Expression.h"
#include <utility>

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

                template < typename TO >
                _ehm_inline
                typename std::enable_if< std::is_same< TO , RET >::value == false , ConvertTo< CLS , TO > >::type
                Convert() const
                {
                    return ConvertTo< CLS , TO >( *this );
                }
                // if TO is same type; then return THIS
                template < typename TO >
                _ehm_inline
                typename std::enable_if< std::is_same< TO , RET >::value , const CLS& >::type
                Convert() const
                {
                    return *this;
                }

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


                _ehm_inline auto Negative() const
                {
                    return NegativeTo< CLS >( *this );
                }
                _ehm_inline auto operator - () const
                {
                    return NegativeTo< CLS >( *this );
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

                // fill with given iterator
                // gettable
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable >::type
                Fill( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator assign" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) = *( begin++ );
                        }
                    }
                }
                // fill with given iterator
                // non-gettable
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable==false >::type
                Fill( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator assign" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) = *( begin++ );
                    }
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

                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Fill( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "assign invalid size of expression" );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) = GetBy( a , i , j );
                        }
                    }
                }
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable==false && Traits< T >::is_gettable==false >::type
                Fill( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "assign invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) = GetBy( a , i );
                    }
                }

                template < IndexType OFFSET = 0 , IndexType MAX = rows , typename SFINE = CLS , typename T0 , typename ... Ts >
                typename std::enable_if< is_column_vector< SFINE >::value >::type
                FillAggressive( T0&& a , Ts&& ... args )
                {
                    static_assert( is_column_vector< T0 >::value , "only vector can perform aggressive-assignment" );

                    constexpr const IndexType _ROWS = Traits< T0 >::rows;
                    for( IndexType i=0; i<_ROWS; ++i )
                    {
                        GetByRef( *this , 0 , OFFSET + i ) = GetBy( a , 0 , i );
                    }
                    FillAggressive< OFFSET + _ROWS , MAX >( std::template forward< Ts >( args )... );
                }
                template < IndexType OFFSET = 0 , IndexType MAX = cols , typename SFINE = CLS , typename T0 , typename ... Ts >
                typename std::enable_if< is_column_vector< SFINE >::value == false >::type
                FillAggressive( T0&& a , Ts&& ... args )
                {
                    static_assert( Traits< T0 >::rows == rows , "invalid row size for non-vector aggresive assign" );

                    constexpr const IndexType _COLS = Traits< T0 >::cols;

                    for( IndexType x=0; x<_COLS; ++x )
                    {
                        for( IndexType i=0; i<rows; ++i )
                        {
                            GetByRef( *this , OFFSET+x , i ) = GetBy( a , x , i );
                        }
                    }
                    FillAggressive< OFFSET + _COLS , MAX >( std::template forward< Ts >( args )... );
                }
                template < IndexType OFFSET = 0 , IndexType MAX >
                _ehm_inline void FillAggressive()
                {
                    static_assert( OFFSET == MAX , "Invalid size for aggresive assignment" );
                }


                // multiply , divide
                // same as Fill()
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Multiply( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Multiply invalid size of expression" );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) *= GetBy( a , i , j );
                        }
                    }
                }
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable==false && Traits< T >::is_gettable==false >::type
                Multiply( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Multiply invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) *= GetBy( a , i );
                    }
                }
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable >::type
                Multiply( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator multiply" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) *= *(begin++);
                        }
                    }
                }
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable==false >::type
                Multiply( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator multiply" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) *= *( begin++ );
                    }
                }

                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Divide( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Divide invalid size of expression" );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) /= GetBy( a , i , j );
                        }
                    }
                }
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable==false && Traits< T >::is_gettable==false >::type
                Divide( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Divide invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) /= GetBy( a , i );
                    }
                }
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable >::type
                Divide( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator divide" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) /= *(begin++);
                        }
                    }
                }
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable == false >::type
                Divide( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator divide" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) /= *( begin++ );
                    }
                }

                // and plus-assign , minus-assign
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Plus( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Plus invalid size of expression" );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) += GetBy( a , i , j );
                        }
                    }
                }
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable==false && Traits< T >::is_gettable==false >::type
                Plus( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Plus invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) += GetBy( a , i );
                    }
                }
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable >::type
                Plus( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator plus" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) += *(begin++);
                        }
                    }
                }
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable == false >::type
                Plus( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator plus" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) += *( begin++ );
                    }
                }

                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Minus( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Minus invalid size of expression" );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) -= GetBy( a , i , j );
                        }
                    }
                }
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable==false && Traits< T >::is_gettable==false >::type
                Minus( T&& a )
                {
                    static_assert( is_scalar< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Minus invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) -= GetBy( a , i );
                    }
                }
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable >::type
                Minus( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator minus" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) -= *(begin++);
                        }
                    }
                }
                template < typename IterType , typename SFINE = CLS >
                typename std::enable_if< Traits< SFINE >::is_gettable == false >::type
                Minus( IterType begin , IterType end )
                {
                    static_assert( std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value ,
                            "invalid type of iterator minus" );
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) -= *( begin++ );
                    }
                }

                // simply forwarding assign operator from interface
                EXPRESSION_ASSIGN_OPERATOR( interface_type )
            };  // struct Expression
        };  // namespace Expression
    };  // namespace Matrix
};  // namespace EH
