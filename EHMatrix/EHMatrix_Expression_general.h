#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_Expression.h"

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
                                        (Traits< TA >::is_restrict==false || Traits< CLS >::is_restrict==false) &&
                                        Traits< TA >::template has_same_root< typename Traits< CLS >::root_type >::value ,
                                        Matrix< ret_type< TA > , Traits< TA >::rows , Traits< TA >::cols > ,
                                        TA
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
                    _ehm_inline
                    auto
                    Column( IndexType i )
                    {
                        return SubMatrixRuntimeExp< CLS , rows , 1 >( *this , i , 0 );
                    }
                    _ehm_inline
                    auto
                    Column( IndexType i ) const
                    {
                        return SubMatrixRuntimeExp< typename std::add_const< CLS >::type , rows , 1 >( *this , i , 0 );
                    }
                    template < IndexType I >
                    _ehm_inline
                    auto
                    Column()
                    {
                        return SubMatrixExp< CLS , I , 0 , rows , 1 >( *this );
                    }
                    template < IndexType I >
                    _ehm_inline
                    auto
                    Column() const
                    {
                        return SubMatrixExp< typename std::add_const< CLS >::type , I , 0 , rows , 1 >( *this );
                    }
                    _ehm_inline
                    auto
                    Row( IndexType i )
                    {
                        return SubMatrixRuntimeExp< CLS , 1 , cols >( *this , 0 , i );
                    }
                    _ehm_inline
                    auto
                    Row( IndexType i ) const
                    {
                        return SubMatrixRuntimeExp< typename std::add_const< CLS >::type , 1 , cols >( *this , 0 , i );
                    }
                    template < IndexType I >
                    _ehm_inline
                    auto
                    Row()
                    {
                        return SubMatrixExp< CLS , 0 , I , 1 , cols >( *this );
                    }
                    template < IndexType I >
                    _ehm_inline
                    auto
                    Row() const
                    {
                        return SubMatrixExp< typename std::add_const< CLS >::type , 0 , I , 1 , cols >( *this );
                    }

                    template < IndexType X , IndexType Y , IndexType M , IndexType N >
                    _ehm_inline
                    auto
                    SubMatrix()
                    {
                        return SubMatrixExp< CLS , X , Y , M , N >( *this );
                    }
                    template < IndexType X , IndexType Y , IndexType M , IndexType N >
                    _ehm_inline
                    auto
                    SubMatrix() const
                    {
                        return SubMatrixExp< typename std::add_const< CLS >::type , X , Y , M , N >( *this );
                    }
                    template < IndexType M ,IndexType N >
                    _ehm_inline
                    auto
                    SubMatrix( IndexType x , IndexType y )
                    {
                        return SubMatrixRuntimeExp< CLS , M , N >( *this , x , y );
                    }
                    template < IndexType M ,IndexType N >
                    _ehm_inline
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
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable
                >::type
                Fill( IterType begin , IterType end )
                {
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
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable==false
                >::type
                Fill( IterType begin , IterType end )
                {
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
                // -> be sure, it may accidently overwrite on the memory
                //
                // or use *_Safe functions to make( or not by it's attributes ) temp object
                //
                // has versions of gettable & non-gettable for performance;

                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Fill( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
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
                Fill( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "assign invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) = GetBy( a , i );
                    }
                }
                template < typename T >
                _ehm_inline
                void
                Fill_Safe( const Expression< T >& exp )
                {
                    Fill< temp_type< T > >( exp );
                    //Fill< typename Expression< T >::temp_type >( exp );
                }

                template < IndexType OFFSET = 0 , typename T0 , typename ... Ts >
                typename std::enable_if<
                    std::is_arithmetic< T0 >::value && is_column_vector< CLS >::value
                >::type
                FillAggressive( const T0 a , Ts&& ... args )
                {
                    GetByRef( *this , 0 , OFFSET ) = a;
                    FillAggressive< OFFSET + 1 >( std::template forward< Ts >( args )... );
                }
                template < IndexType OFFSET = 0 , typename T0 , typename ... Ts >
                typename std::enable_if<
                    std::is_arithmetic< T0 >::value==false && is_column_vector< T0 >::value && is_column_vector< CLS >::value
                >::type
                FillAggressive( const Expression< T0 >& a , Ts&& ... args )
                {
                    for( IndexType i=0; i<(Traits< T0 >::rows); ++i )
                    {
                        GetByRef( *this , 0 , OFFSET + i ) = GetBy( a , 0 , i );
                    }
                    FillAggressive< OFFSET + (Traits< T0 >::rows) >( std::template forward< Ts >( args )... );
                }
                template < IndexType OFFSET = 0 >
                _ehm_inline void FillAggressive()
                {
                    static_assert( OFFSET == rows , "Invalid size for aggresive assignment" );
                }


                // multiply , divide
                // same as Fill()
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Multiply( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
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
                Multiply( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Multiply invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) *= GetBy( a , i );
                    }
                }
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable
                >::type
                Multiply( IterType begin , IterType end )
                {
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) *= *(begin++);
                        }
                    }
                }
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable==false
                >::type
                Multiply( IterType begin , IterType end )
                {
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) *= *( begin++ );
                    }
                }
                template < typename T >
                _ehm_inline
                void
                Multiply_Safe( const Expression< T >& exp )
                {
                    Multiply< temp_type< T > >( exp );
                }

                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Divide( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
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
                Divide( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Divide invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) /= GetBy( a , i );
                    }
                }
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable
                >::type
                Divide( IterType begin , IterType end )
                {
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) /= *(begin++);
                        }
                    }
                }
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable==false
                >::type
                Divide( IterType begin , IterType end )
                {
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) /= *( begin++ );
                    }
                }
                template < typename T >
                _ehm_inline
                void
                Divide_Safe( const Expression< T >& exp )
                {
                    Divide< temp_type< T > >( exp );
                }

                // and plus-assign , minus-assign
                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Plus( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
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
                Plus( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Plus invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) += GetBy( a , i );
                    }
                }
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable
                >::type
                Plus( IterType begin , IterType end )
                {
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) += *(begin++);
                        }
                    }
                }
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable==false
                >::type
                Plus( IterType begin , IterType end )
                {
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) += *( begin++ );
                    }
                }
                template < typename T >
                _ehm_inline
                void
                Plus_Safe( const Expression< T >& exp )
                {
                    Plus< temp_type< T > >( exp );
                }

                template < typename T >
                typename std::enable_if< Traits< CLS >::is_gettable || Traits< T >::is_gettable >::type
                Minus( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
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
                Minus( auto_creference< T > a )
                {
                    static_assert( std::is_arithmetic< T >::value || ( Traits< T >::rows == rows && Traits< T >::cols == cols ) ,
                            "Minus invalid size of expression" );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) -= GetBy( a , i );
                    }
                }
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable
                >::type
                Minus( IterType begin , IterType end )
                {
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols; ++i )
                    {
                        for( IndexType j=0; j<rows; ++j )
                        {
                            GetByRef( *this , i , j ) -= *(begin++);
                        }
                    }
                }
                template < typename IterType >
                typename std::enable_if<
                    std::is_convertible< typename std::iterator_traits< IterType >::value_type , RET >::value &&
                    Traits< CLS >::is_gettable==false
                >::type
                Minus( IterType begin , IterType end )
                {
                    assert( std::distance( begin , end ) == rows*cols );
                    for( IndexType i=0; i<cols*rows; ++i )
                    {
                        GetByRef( *this , i ) -= *( begin++ );
                    }
                }
                template < typename T >
                _ehm_inline
                void
                Minus_Safe( const Expression< T >& exp )
                {
                    Minus< temp_type< T > >( exp );
                }

                // simply forwarding assign operator from interface
                EXPRESSION_ASSIGN_OPERATOR( interface_type )
            };  // struct Expression
        };  // namespace Expression
    };  // namespace Matrix
};  // namespace EH
