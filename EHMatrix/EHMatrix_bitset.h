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
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator == ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return BinaryExp< TA , TB , bool >(
                        m1 , m2 ,
                        []( const auto a , const auto b )->bool
                        {
                            return a == b;
                        }
                    );
            }
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator != ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return BinaryExp< TA , TB , bool >(
                        m1 , m2 ,
                        []( const auto a , const auto b )->bool
                        {
                            return a != b;
                        }
                    );
            }
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator > ( const expression_size_type< TA , M , N >& m1 ,
                              const expression_size_type< TB , M , N >& m2 )
            {
                return BinaryExp< TA , TB , bool >(
                        m1 , m2 ,
                        []( const auto a , const auto b )->bool
                        {
                            return a > b;
                        }
                    );
            }
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator < ( const expression_size_type< TA , M , N >& m1 ,
                              const expression_size_type< TB , M , N >& m2 )
            {
                return BinaryExp< TA , TB , bool >(
                        m1 , m2 ,
                        []( const auto a , const auto b )->bool
                        {
                            return a < b;
                        }
                    );
            }
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator >= ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return BinaryExp< TA , TB , bool >(
                        m1 , m2 ,
                        []( const auto a , const auto b )->bool
                        {
                            return a >= b;
                        }
                    );
            }
            template < typename TA , typename TB ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator <= ( const expression_size_type< TA , M , N >& m1 ,
                               const expression_size_type< TB , M , N >& m2 )
            {
                return BinaryExp< TA , TB , bool >(
                        m1 , m2 ,
                        []( const auto a , const auto b )->bool
                        {
                            return a <= b;
                        }
                    );
            }
            template < typename TA >
            auto operator == ( const Expression< TA >& m1 ,
                               const ret_type< TA > s )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a == s;
                        }
                    );
            }
            template < typename TA >
            auto operator != ( const Expression< TA >& m1 ,
                               const ret_type< TA > s )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a != s;
                        }
                    );
            }
            template < typename TA >
            auto operator >= ( const Expression< TA >& m1 ,
                               const ret_type< TA > s )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a >= s;
                        }
                    );
            }
            template < typename TA >
            auto operator <= ( const Expression< TA >& m1 ,
                               const ret_type< TA > s )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a <= s;
                        }
                    );
            }
            template < typename TA >
            auto operator > ( const Expression< TA >& m1 ,
                              const ret_type< TA > s )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a > s;
                        }
                    );
            }
            template < typename TA >
            auto operator < ( const Expression< TA >& m1 ,
                              const ret_type< TA > s )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a < s;
                        }
                    );
            }

            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator == ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a == s;
                        }
                    );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator != ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a != s;
                        }
                    );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator >= ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a <= s;
                        }
                    );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator <= ( const ret_type< TA > s ,
                               const expression_size_type< TA , M , N >& m1 )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a >= s;
                        }
                    );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator > ( const ret_type< TA > s ,
                              const expression_size_type< TA , M , N >& m1 )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a < s;
                        }
                    );
            }
            template < typename TA ,
                       IndexType M = Traits< TA >::rows ,
                       IndexType N = Traits< TA >::cols >
            auto operator < ( const ret_type< TA > s ,
                              const expression_size_type< TA , M , N >& m1 )
            {
                return UnaryExp< TA , bool >(
                        m1 ,
                        [ s ]( const auto a )->bool
                        {
                            return a > s;
                        }
                    );
            }
        };  // namespace Expression
    };  // namespace Matrix
};  // namespace EH
