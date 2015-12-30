#pragma once

#include "EHMatrix_Global.h"
#include "EHMatrix_Expression.h"
#include <bitset>

namespace EH
{
    namespace Matrix
    {
        namespace Expression
        {
            template < typename TA , typename TB , typename CRTP >
            struct bitsetExp
            {
                constexpr const static bool is_gettable =
                    Traits< TA >::is_gettable || Traits< TB >::is_gettable;
                constexpr const static IndexType rows =
                    std::max( Traits< TA >::rows , Traits< TB >::rows );
                constexpr const static IndexType cols =
                    std::max( Traits< TA >::cols , Traits< TB >::cols );
                constexpr const static IndexType N = rows * cols;

                typedef std::bitset< N > stdbitset;

                _ehm_inline
                bool
                test( IndexType i ) const
                {
                    return static_cast< const CRTP& >( *this ).test( i );
                }
                template < typename THIS = bitsetExp< TA , TB , CRTP > >
                typename std::enable_if< THIS::is_gettable == false , bool >::type
                all() const
                {
                    bool bo = true;
                    for( IndexType i=0; i<N; ++i )
                    {
                        //if( test( i ) == false ){ return false; }
                        bo &= test( i );
                    }
                    //return true;
                    return bo;
                }
                template < typename THIS = bitsetExp< TA , TB , CRTP > >
                typename std::enable_if< THIS::is_gettable == false , bool >::type
                any() const
                {
                    bool bo = false;
                    for( IndexType i=0; i<N; ++i )
                    {
                        //if( test( i ) ){ return true; }
                        bo |= test( i );
                    }
                    //return false;
                    return bo;
                }
                template < typename THIS = bitsetExp< TA , TB , CRTP > >
                typename std::enable_if< THIS::is_gettable == false , bool >::type
                none() const
                {
                    bool bo = false;
                    for( IndexType i=0; i<N; ++i )
                    {
                        bo |= test( i );
                    }
                    return !bo;
                }
                template < typename THIS = bitsetExp< TA , TB , CRTP > , typename = void ,
                           typename = typename std::enable_if< THIS::is_gettable == false >::type
                           >
                operator stdbitset () const
                {
                    stdbitset bit;
                    for( IndexType i=0; i<N; ++i )
                    {
                        bit.set( i , test( i ) );
                    }
                    return bit;
                }



                // for gettable structures;
                _ehm_inline
                bool
                test( IndexType x , IndexType y ) const
                {
                    return static_cast< const CRTP& >( *this ).test( x , y );
                }
                template< typename THIS = bitsetExp< TA , TB , CRTP > >
                typename std::enable_if< THIS::is_gettable , bool >::type
                all() const
                {
                    bool bo = true;
                    for( IndexType c=0; c<cols; ++c )
                    {
                        for( IndexType r=0; r<rows; ++r )
                        {
                            bo &= test( c , r );
                        }
                    }
                    return bo;
                }
                template< typename THIS = bitsetExp< TA , TB , CRTP > >
                typename std::enable_if< THIS::is_gettable , bool >::type
                any() const
                {
                    bool bo = false;
                    for( IndexType c=0; c<cols; ++c )
                    {
                        for( IndexType r=0; r<rows; ++r )
                        {
                            bo |= test( c , r );
                        }
                    }
                    return bo;
                }
                template< typename THIS = bitsetExp< TA , TB , CRTP > >
                typename std::enable_if< THIS::is_gettable , bool >::type
                none() const
                {
                    bool bo = false;
                    for( IndexType c=0; c<cols; ++c )
                    {
                        for( IndexType r=0; r<rows; ++r )
                        {
                            bo |= test( c , r );
                        }
                    }
                    return !bo;
                }

                template < typename THIS = bitsetExp< TA , TB , CRTP > ,
                           typename = typename std::enable_if< THIS::is_gettable >::type
                           >
                operator stdbitset () const
                {
                    stdbitset bit;
                    for( IndexType c=0; c<cols; ++c )
                    {
                        for( IndexType r=0; r<rows; ++r )
                        {
                            bit.set( r + c*rows , test( c , r ) );
                        }
                    }
                    return bit;
                }

                _ehm_inline operator bool () const
                {
                    return all();
                }
            };

            template < typename TA , typename TB , bool MASK >
            struct IsEquals :
                bitsetExp< TA , TB , IsEquals< TA , TB , MASK > > ,
                Expression< IsEquals< TA , TB , MASK > >
            {
                auto_creference< TA > a;
                auto_creference< TB > b;

                IsEquals( auto_creference< TA > _a , auto_creference< TB > _b ) :
                    a( _a ) , b( _b )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >* ptr ) const
                {
                    return a.has_same_root( ptr ) || b.has_same_root( ptr );
                }

                _ehm_inline
                bool
                operator [] ( IndexType i ) const
                {
                    return test( i );
                }
                _ehm_inline
                bool
                Get( IndexType x , IndexType y ) const
                {
                    return test( x , y );
                }

                _ehm_inline bool test( IndexType i ) const
                {
                    return ( GetBy( a , i ) == GetBy( b , i ) )==MASK;
                }
                _ehm_inline bool test( IndexType x , IndexType y ) const
                {
                    return ( GetBy( a , x , y ) == GetBy( b , x , y ) )==MASK;
                }
            };

            template < typename TA , typename TB , bool MASK >
            struct IsLess :
                bitsetExp< TA , TB , IsLess< TA , TB , MASK > > ,
                Expression< IsLess< TA , TB , MASK > >
            {
                auto_creference< TA > a;
                auto_creference< TB > b;

                IsLess( auto_creference< TA > _a , auto_creference< TB > _b ) :
                    a( _a ) , b( _b )
                {
                }
                template < typename T2 , IndexType M2 , IndexType N2 >
                _ehm_inline bool has_same_root( const Matrix< T2 , M2 , N2 >* ptr ) const
                {
                    return a.has_same_root( ptr ) || b.has_same_root( ptr );
                }

                _ehm_inline
                bool
                operator [] ( IndexType i ) const
                {
                    return test( i );
                }
                _ehm_inline
                bool
                Get( IndexType x , IndexType y ) const
                {
                    return test( x , y );
                }
                _ehm_inline bool test( IndexType i ) const
                {
                    return ( GetBy( a , i ) < GetBy( b , i ) )==MASK;
                }
                _ehm_inline bool test( IndexType x , IndexType y ) const
                {
                    return ( GetBy( a , x , y ) < GetBy( b , x , y ) )==MASK;
                }
            };
        };  // namespace Expression;
        namespace Expression
        {
            template < typename TA , typename TB , bool MASK >
            struct Traits< IsEquals< TA , TB , MASK > > : TraitsCombine< TA , TB >
            {
                typedef bool return_type;
                constexpr const static bool catch_reference = false;
            };
            template < typename TA , typename TB , bool MASK >
            struct Traits< IsLess< TA , TB , MASK > > : TraitsCombine< TA , TB >
            {
                typedef bool return_type;
                constexpr const static bool catch_reference = false;
            };

            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator == ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return IsEquals< TA , TB , true >( m1 , m2 );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator == ( const expression_size_type< TA , M , N >& m1 ,
                               const ret_type< TA > s )
            {
                return IsEquals< TA , ret_type< TA > , true >( m1 , s );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator == ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return IsEquals< ret_type< TA > , TA , true >( s , m1 );
            }

            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator != ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return IsEquals< TA , TB , false >( m1 , m2 );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator != ( const expression_size_type< TA , M , N >& m1 ,
                               const ret_type< TA > s )
            {
                return IsEquals< TA , ret_type< TA > , false >( m1 , s );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator != ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return IsEquals< ret_type< TA > , TA , false >( s , m1 );
            }


            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator < ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return IsLess< TA , TB , true >( m1 , m2 );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator < ( const expression_size_type< TA , M , N >& m1 ,
                               const ret_type< TA > s )
            {
                return IsLess< TA , ret_type< TA > , true >( m1 , s );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator < ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return IsLess< ret_type< TA > , TA , true >( s , m1 );
            }


            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator <= ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return IsLess< TB , TA , false >( m2 , m1 );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator <= ( const expression_size_type< TA , M , N >& m1 ,
                               const ret_type< TA > s )
            {
                return IsLess< ret_type< TA > , TA , false >( s , m1 );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator <= ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return IsLess< TA , ret_type< TA > , false >( m1 , s );
            }

            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator > ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return IsLess< TB , TA , true >( m2 , m1 );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator > ( const expression_size_type< TA , M , N >& m1 ,
                               const ret_type< TA > s )
            {
                return IsLess< ret_type< TA > , TA , true >( s , m1 );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator > ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return IsLess< TA , ret_type< TA > , true >( m1 , s );
            }


            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator >= ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return IsLess< TA , TB , false >( m1 , m2 );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator >= ( const expression_size_type< TA , M , N >& m1 ,
                               const ret_type< TA > s )
            {
                return IsLess< TA , ret_type< TA > , false >( m1 , s );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator >= ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return IsLess< ret_type< TA > , TA , false >( s , m1 );
            }
        };  // namespace Expression
    };  // namespace Matrix
};  // namespace EH
