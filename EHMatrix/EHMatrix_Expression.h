#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_Expression_functions.h"

#include <initializer_list>
#include <functional>



namespace EH
{
    namespace Matrix
    {
        namespace Expression
        {
            // default structure for Traits
            template < typename T >
            struct Traits
            {
                typedef typename std::enable_if< is_scalar< T >::value , int >::type Enabled;
                constexpr const static Enabled enabled_v = 0;

                constexpr const static bool is_gettable = false;
                constexpr const static IndexType cols = 1;
                constexpr const static IndexType rows = 1;
                typedef T return_type;


                // if restrict is false;
                // at assign make temporary matrix to copy from;
                // otherwise just assigning wit index
                constexpr const static bool is_restrict = true;

                // number of operations performed by this expression
                constexpr const static int operations = 0;

                constexpr const static bool catch_reference = false;

                typedef T root_type;

                template < typename CLS >
                struct has_same_root
                {
                    constexpr const static bool value = false;
                };
            };

            template < typename TA , typename RET = void , int OPADD = 1 >
            struct UnaryExp :
                Expression< UnaryExp< TA , RET , OPADD > >
            {
                typedef typename Traits< UnaryExp< TA , RET , OPADD > >::return_type _RET;
                auto_creference< TA > a;

                std::function< _RET( typename Traits< TA >::return_type ) > func;

                template < typename FUNC >
                UnaryExp( auto_creference< TA > _a , FUNC&& _func )
                    : a( _a ) , func( _func )
                {
                }

                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr );
                }

                _ehm_inline _RET operator [] ( IndexType i ) const
                {
                    //return GetBy( a , i );
                    return func( GetBy( a , i ) );
                }
                _ehm_inline _RET Get( IndexType x , IndexType y ) const
                {
                    //return GetBy( a  ,x , y );
                    return func( GetBy( a , x , y ) );
                }
            };
            template < typename TA , typename RET , int OPADD >
            struct Traits< UnaryExp< TA , RET , OPADD > > : Traits< TA >
            {
                typedef typename std::conditional<
                            std::is_same< RET , void >::value ,
                            typename Traits< TA >::return_type ,
                            RET
                        >::type return_type;
                constexpr const static bool catch_reference = false;
                constexpr const static int operations = Traits< TA >::operations + OPADD;
            };


            template < typename TA , typename TB , typename RET = void , int OPADD = 1 >
            struct BinaryExp :
                Expression< BinaryExp< TA , TB , RET , OPADD > >
            {
                typedef typename Traits< BinaryExp< TA , TB , RET , OPADD > >::return_type _RET;

                auto_creference< TA > a;
                auto_creference< TB > b;

                std::function< _RET( typename Traits< TA >::return_type , typename Traits< TB >::return_type ) > func;

                template < typename FUNC >
                BinaryExp( auto_creference< TA > _a , auto_creference< TB > _b , FUNC&& _func )
                    : a( _a ) , b( _b ) , func( _func )
                {
                }

                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr ) || HasSameRoot( b , ptr );
                }

                _ehm_inline auto operator [] ( IndexType i ) const
                {
                    //return GetBy( a , i );
                    return func( GetBy( a , i ) , GetBy( b , i ) );
                }
                _ehm_inline auto Get( IndexType x , IndexType y ) const
                {
                    //return GetBy( a , x , y );
                    return func( GetBy( a , x , y ) , GetBy( b , x , y ) );
                }
            };
            template < typename TA , typename TB , typename RET , int OPADD >
            struct Traits< BinaryExp< TA , TB , RET , OPADD > > : TraitsCombine< TA , TB >
            {
                typedef typename std::conditional<
                            std::is_same< RET , void >::value ,
                            typename std::common_type<
                                typename Traits< TA >::return_type , typename Traits< TB >::return_type
                            >::type ,
                            RET
                        >::type return_type;
                constexpr const static bool catch_reference = false;
                constexpr const static int operations = Traits< TA >::operations + Traits< TB >::operations + OPADD;
            };

            template < typename TA >
            struct TransposeTo;

            template < typename TA , typename TB >
            struct SquarePlus;

            // fills 0 if it's out of bound
            // fills 1 if it's out of bound && x==y
            // else copy from given matrix
            template < typename TA , IndexType M , IndexType N >
            struct GeneralResize;
            template < typename TA , IndexType M , IndexType N >
            struct SubMatrixRuntimeExp;
            template < typename TA >
            struct DiagonalExp;
            template < typename TA , typename TB >
            struct ShuffleExp;

            template< IndexType N , typename ... Ts >
            struct OuterProduct;

            template < typename TA , typename TB >
            struct ComplexMultiply;
            template < typename TA >
            struct ComplexConjugate;
            template < typename TA , typename TB >
            struct QuatMult;
            template < typename TA >
            struct QuatConjugate;

            // CLS : class where operator[] is defined;
            template < typename TA >
            struct TransposeTo :
                Expression< TransposeTo< TA > >
            {
                auto_reference< TA > a;

                TransposeTo( auto_reference< TA > _a )
                    : a( _a )
                {
                }

                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr );
                }
                _ehm_inline auto& operator [] ( IndexType i )
                {
                    return GetByRef( a , i );
                }
                _ehm_inline auto& Get( IndexType x , IndexType y )
                {
                    return GetByRef( a , y , x );
                }
                _ehm_inline auto operator [] ( IndexType i ) const
                {
                    return GetBy( a , i );
                }
                _ehm_inline auto Get( IndexType x , IndexType y ) const
                {
                    return GetBy( a , y , x );
                }

                typedef Expression< TransposeTo< TA > > parent;
                EXPRESSION_ASSIGN_OPERATOR( parent )
            };

            template < typename TA >
            struct Traits< TransposeTo< TA > > : Traits< TA >
            {
                constexpr const static IndexType cols = Traits< TA >::rows;
                constexpr const static IndexType rows = Traits< TA >::cols;
                constexpr const static bool is_gettable = Traits< TA >::is_gettable || is_vector< TA >::value == false;
                constexpr const static bool is_restrict = false;
                constexpr const static bool catch_reference = false;
            };

            template < typename TA , typename TB >
            struct MatMatMult :
                Expression< MatMatMult< TA , TB > >
            {
                typedef typename std::common_type< ret_type< TA > , ret_type< TB > >::type return_type;

                typename ShouldMakeTemp< TA , Traits< TB >::cols >::type a;
                typename ShouldMakeTemp< TB , Traits< TA >::rows >::type b;

                MatMatMult( auto_creference< TA > _a , auto_creference< TB > _b )
                    : a( _a ) , b( _b )
                {
                }

                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr ) || HasSameRoot( b , ptr );
                }
                // operator [] can be only used when return is vector;
                //
                // temporary template typename for SFINE


                // vec-mat mult
                template < typename _TA = TA , typename _TB = TB >
                typename std::enable_if< Traits< _TA >::rows == 1 && Traits< _TB >::cols != 1 , return_type >::type
                operator [] ( IndexType i ) const
                {
                    return_type sum = return_type( 0 );

                    for( IndexType j=0; j<Traits< TA >::cols; ++j )
                    {
                        sum += GetBy( a , j , 0 ) * GetBy( b , i , j );
                    }

                    return sum;
                }
                // mat-vec mult
                template < typename _TA = TA , typename _TB = TB >
                typename std::enable_if< Traits< _TA >::rows != 1 && Traits< _TB >::cols == 1 , return_type >::type
                operator [] ( IndexType i ) const
                {
                    return_type sum = return_type( 0 );

                    for( IndexType j=0; j<Traits< TA >::cols; ++j )
                    {
                        sum += GetBy( a , j , i ) * GetBy( b , 0 , j );
                    }

                    return sum;
                }
                // and the others
                return_type Get( IndexType x , IndexType y ) const
                {
                    return_type sum = return_type( 0 );

                    for( IndexType j=0; j<Traits< TA >::cols; ++j )
                    {
                        sum += GetBy( a , j , y ) * GetBy( b , x , j );
                    }

                    return sum;
                }
            };  // struct MatMatMult

            template < typename TA , typename TB >
            struct Traits< MatMatMult< TA , TB > > : TraitsCombine< TA , TB >
            {
                constexpr const static IndexType cols = Traits< TB >::cols;
                constexpr const static IndexType rows = Traits< TA >::rows;
                // if return is vector ( either row or col ) can be index-accessed
                constexpr const static bool is_gettable =
                    cols != 1 && rows != 1;
                // restrict is true only if both made temp object
                constexpr const static bool is_restrict =
                    ( ShouldMakeTemp< TA , Traits< TB >::cols >::value &&
                      ShouldMakeTemp< TB , Traits< TA >::rows >::value );
                constexpr const static int operations =
                    ( ShouldMakeTemp< TA , Traits< TB >::cols >::operations + ShouldMakeTemp< TB , Traits< TA >::rows >::operations + 1 )
                    * Traits< TA >::cols - 1;
                constexpr const static bool catch_reference =
                    ShouldMakeTemp< TA , Traits< TB >::cols >::value || ShouldMakeTemp< TB , Traits< TA >::rows >::value;
            };

            template < typename TA , typename TB >
            MatMatMult< remove_cr< TA > , remove_cr< TB > >
            make_matmatmult( TA&& a , TB&& b )
            {
                return MatMatMult< remove_cr< TA > , remove_cr< TB > >(
                        std::template forward< TA >( a ) ,
                        std::template forward< TB >( b )
                    );
            }



            template < typename TA , IndexType M , IndexType N >
            struct GeneralResize :
                Expression< DiagonalExp< TA > >
            {
                auto_creference< TA > a;

                GeneralResize( auto_creference< TA > _a ) :
                    a( _a )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr );
                }

                template < typename SFINE = TA >
                _ehm_inline
                typename std::enable_if< (N <= Traits< SFINE >::cols) , ret_type< TA > >::type
                operator [] ( IndexType i ) const
                {
                    return GetBy( a , i );
                }
                template < typename SFINE = TA >
                _ehm_inline
                typename std::enable_if< (N > Traits< SFINE >::cols) , ret_type< TA > >::type
                operator [] ( IndexType i ) const
                {
                    if( i < Traits< TA >::rows*Traits< TA >::cols )
                    {
                        return GetBy( a , i );
                    }
                    return (i/Traits< TA >::rows)==(i%(Traits< TA >::rows)) ? 1 : 0;
                }
                template < typename SFINE = TA >
                _ehm_inline
                typename std::enable_if< ( M<=Traits< SFINE >::rows && N<=Traits< SFINE >::cols ) , ret_type< TA > >::type
                Get ( IndexType x , IndexType y ) const
                {
                    return GetBy( a , x , y );
                }
                template < typename SFINE = TA >
                _ehm_inline
                typename std::enable_if< ( M>Traits< SFINE >::rows || N>Traits< SFINE >::cols ) , ret_type< TA > >::type
                Get ( IndexType x , IndexType y ) const
                {
                    if( (x < Traits< TA >::cols) && (y < Traits< TA >::rows) )
                    {
                        return GetBy( a , x , y );
                    }
                    return x==y ? 1 : 0;
                }
            };
            template < typename TA , IndexType M , IndexType N >
            struct Traits< GeneralResize< TA , M , N > > : Traits< TA >
            {
                // case it's index-accessable :
                //  - is not gettable AND
                //  -same height( rows count ) AND
                //      -output's width is less-same than input OR
                //      -rows+1 is power of two( thus can parse x and y coordinated easily )
                constexpr const static bool is_gettable = Traits< TA >::is_gettable ||
                    M != Traits< TA >::rows || ( N > Traits< TA >::cols && M&(M+1) != 0 );
                constexpr const static bool is_restrict = false;
                constexpr const static IndexType cols = N;
                constexpr const static IndexType rows = M;
                constexpr const static bool catch_reference = false;
            };

            template < typename TA >
            struct DiagonalExp :
                Expression< DiagonalExp< TA > >
            {
                auto_reference< TA > a;

                DiagonalExp( auto_reference< TA > _a ) :
                    a( _a )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr );
                }

                _ehm_inline ret_type< TA >& operator [] ( IndexType i )
                {
                    return GetByRef( a , i , i );
                }
                _ehm_inline ret_type< TA >& Get( IndexType x , IndexType y )
                {
                    return GetByRef( a , y , y );
                }
                _ehm_inline ret_type< TA > operator [] ( IndexType i ) const
                {
                    return GetBy( a , i , i );
                }
                _ehm_inline ret_type< TA > Get( IndexType x , IndexType y ) const
                {
                    return GetBy( a , y , y );
                }

                typedef Expression< DiagonalExp< TA > > parent;
                EXPRESSION_ASSIGN_OPERATOR( parent );
            };
            template < typename TA >
            struct Traits< DiagonalExp< TA > > : Traits< TA >
            {
                constexpr const static bool is_gettable = false;
                constexpr const static IndexType cols = 1;
                constexpr const static IndexType rows = std::min( Traits< TA >::cols , Traits< TA >::rows );
                constexpr const static bool is_restrict = false;
                constexpr const static bool catch_reference = false;
            };

            template < typename TA , IndexType M , IndexType N >
            struct SubMatrixRuntimeExp :
                Expression< SubMatrixRuntimeExp< TA , M , N > >
            {
                auto_reference< TA > a;
                const IndexType X;
                const IndexType Y;

                constexpr SubMatrixRuntimeExp( auto_reference< TA > _a , const IndexType _x , const IndexType _y ) :
                    a( _a ) , X( _x ) , Y( _y )
                {
                }

                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr );
                }
                constexpr _ehm_inline ret_type< TA >& operator [] ( IndexType i )
                {
                    return GetByRef( a , i + X*Traits< TA >::rows );
                }
                constexpr _ehm_inline ret_type< TA >& Get( IndexType x , IndexType y )
                {
                    return GetByRef( a , x + X , y + Y );
                }
                constexpr _ehm_inline ret_type< TA > operator [] ( IndexType i ) const
                {
                    return GetBy( a , i  + X*Traits< TA >::rows );
                }
                constexpr _ehm_inline ret_type< TA > Get( IndexType x , IndexType y ) const
                {
                    return GetBy( a , x + X , y + Y );
                }

                typedef Expression< SubMatrixRuntimeExp< TA , M , N > > parent;
                EXPRESSION_ASSIGN_OPERATOR( parent );
            };

            template < typename TA , IndexType M , IndexType N >
            struct Traits< SubMatrixRuntimeExp< TA , M , N > > : Traits< TA >
            {
                // case it's index accessable;
                //   -original is index accessable AND
                //      -original is vector type
                //      -same M( rows count ) as original

                constexpr const static bool is_gettable = Traits< TA >::is_gettable ||
                    ( M != Traits< TA >::rows && is_vector< TA >::value == false );
                constexpr const static bool is_restrict = false;
                constexpr const static IndexType cols = N;
                constexpr const static IndexType rows = M;
                constexpr const static bool catch_reference = false;
            };


            template < typename TA , typename TB >
            struct ShuffleExp :
                Expression< ShuffleExp< TA , TB > >
            {
                // make sure operand is index-accessable;
                typedef typename std::enable_if< Traits< TA >::is_gettable==false > Enabled;

                auto_reference< TA > a;
                auto_creference< TB > b;

                ShuffleExp( auto_reference< TA > _a , auto_creference< TB > _b ) :
                    a( _a ) , b( _b )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                constexpr _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >& ptr ) const
                {
                    return HasSameRoot( a , ptr ) || HasSameRoot( b , ptr );
                }

                _ehm_inline ret_type< TA >& operator [] ( IndexType i )
                {
                    return GetByRef( a , GetBy( b , i ) );
                }
                _ehm_inline ret_type< TA >& Get( IndexType x , IndexType y )
                {
                    return GetByRef( a , GetBy( b , x , y ) );
                }
                _ehm_inline ret_type< TA > operator [] ( IndexType i ) const
                {
                    return GetBy( a , GetBy( b , i ) );
                }
                _ehm_inline ret_type< TA > Get( IndexType x , IndexType y ) const
                {
                    return GetBy( a , GetBy( b , x , y ) );
                }

                typedef Expression< ShuffleExp< TA , TB > > parent;
                EXPRESSION_ASSIGN_OPERATOR( parent );
            };
            template < typename TA , typename TB >
            struct Traits< ShuffleExp< TA , TB > > : Traits< TA >
            {
                constexpr const static bool is_gettable = Traits< TB >::is_gettable;
                constexpr const static bool is_restrict = false;
                constexpr const static IndexType cols = Traits< TB >::cols;
                constexpr const static IndexType rows = Traits< TA >::rows;
                constexpr const static IndexType operations = Traits< TA >::operations + Traits< TB >::operations;
                constexpr const static bool catch_reference = false;
            };

        };  // namespace Expression

    };  // namespace Matrix
};  // namespace EH
