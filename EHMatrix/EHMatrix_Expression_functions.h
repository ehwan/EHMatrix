#pragma once

#include "EHMatrix_Global.h"

namespace EH
{
    namespace Matrix
    {
        namespace Expression
        {
            template < typename T >
            using ret_type = typename Traits< T >::return_type;

            template < typename CLS >
            using mat_type = Matrix< typename Traits< CLS >::return_type , Traits< CLS >::rows , Traits< CLS >::cols >;

            template < typename CLS >
            struct is_column_vector
            {
                constexpr const static bool value = Traits< CLS >::cols == 1;
            };
            template < typename CLS >
            struct is_row_vector
            {
                constexpr const static bool value = Traits< CLS >::rows == 1;
            };
            template < typename CLS >
            struct is_vector
            {
                constexpr const static bool value = is_column_vector< CLS >::value || is_row_vector< CLS >::value;
            };
            template < typename CLS >
            struct is_square
            {
                constexpr const static bool value = Traits< CLS >::rows == Traits< CLS >::cols;
            };

            template < typename CLS >
            struct Traits< CLS const > : Traits< CLS >
            {
            };
            template < typename CLS >
            struct Traits< CLS& > : Traits< CLS >
            {
            };
            template < typename CLS >
            struct Traits< CLS const& > : Traits< CLS >
            {
            };

            template < typename CLS >
            struct Traits< Expression< CLS > > : Traits< CLS >
            {
            };

            template < typename ... Ts >
            struct TraitsCombine;
            template < typename TA , typename TB >
            struct TraitsCombine< TA , TB >
            {
                constexpr const static bool is_gettable =
                    Traits< TA >::is_gettable || Traits< TB >::is_gettable;
                constexpr const static bool is_restrict =
                    Traits< TA >::is_restrict && Traits< TB >::is_restrict;
                constexpr const static IndexType cols =
                    std::max( Traits< TA >::cols , Traits< TB >::cols );
                constexpr const static IndexType rows =
                    std::max( Traits< TA >::rows , Traits< TB >::rows );
                constexpr const static int operations = Traits< TA >::operations + Traits< TB >::operations + 1;
                typedef typename std::common_type<
                            typename Traits< TA >::return_type ,
                            typename Traits< TB >::return_type
                        >::type return_type;

                template < typename CLS >
                struct has_same_root
                {
                    constexpr const static bool value =
                        Traits< TA >::template has_same_root< CLS >::value || Traits< TB >::template has_same_root< CLS >::value;
                };
                typedef typename Traits< TA >::root_type root_type;
            };
            template < typename T0 , typename ... Ts >
            struct TraitsCombine< T0 , Ts... >
            {
                constexpr const static bool is_gettable =
                    Traits< T0 >::is_gettable || TraitsCombine< Ts... >::is_gettable;
                constexpr const static bool is_restrict =
                    Traits< T0 >::is_restrict && TraitsCombine< Ts... >::is_restrict;
                constexpr const static IndexType cols =
                    std::max( Traits< T0 >::cols , TraitsCombine< Ts... >::cols );
                constexpr const static IndexType rows =
                    std::max( Traits< T0 >::rows , TraitsCombine< Ts... >::rows );
                constexpr const static int operations =
                    Traits< T0 >::operations + TraitsCombine< Ts... >::operations;
                typedef typename std::common_type<
                            typename Traits< T0 >::return_type ,
                            typename TraitsCombine< Ts... >::return_type
                        >::type return_type;

                template < typename CLS >
                struct has_same_root
                {
                    constexpr const static bool value =
                        Traits< T0 >::template has_same_root< CLS >::value || TraitsCombine< Ts... >::template has_same_root< CLS >::value;
                };

                typedef typename Traits< T0 >::root_type root_type;
            };

            // bit mask of conditional expression
            // choose whick attribute to compare;
            struct exm
            {
                constexpr const static unsigned int type       = 0b1;
                constexpr const static unsigned int rows       = 0b10;
                constexpr const static unsigned int cols       = 0b100;
                constexpr const static unsigned int isgettable = 0b1000;
                constexpr const static unsigned int isrestrict = 0b10000;

                constexpr const static unsigned int size     = rows | cols;
                constexpr const static unsigned int all      = type | size | isgettable | isrestrict;
                constexpr const static unsigned int sizetype = size | type;
            };
            template < typename CLS ,
                       typename T ,
                       IndexType M ,
                       IndexType N ,
                       bool ISGETTABLE ,
                       bool ISRESTRICT ,
                       unsigned int MASK = exm::all >
            struct conditional_expression
            {
                constexpr const static bool type =
                         ( std::is_same< typename Traits< CLS >::return_type , T >::value || std::is_same< void , T >::value );
                constexpr const static bool cols =
                         ( N==0 || Traits< CLS >::cols == N );
                constexpr const static bool rows =
                         ( M==0 || Traits< CLS >::rows == M );
                constexpr const static bool isgettable =
                         ( Traits< CLS >::is_gettable == ISGETTABLE );
                constexpr const static bool isrestrict =
                         ( Traits< CLS >::is_restrict == ISRESTRICT );

                    constexpr const static bool value =
                        ( ( MASK & exm::type       )==0 || type       ) &&
                        ( ( MASK & exm::rows       )==0 || rows       ) &&
                        ( ( MASK & exm::cols       )==0 || cols       ) &&
                        ( ( MASK & exm::isgettable )==0 || isgettable ) &&
                        ( ( MASK & exm::isrestrict )==0 || isrestrict );

                static _ehm_inline void Assert()
                {
                    static_assert( type ,       "Invalid Type" );
                    static_assert( cols ,       "Invalid Width" );
                    static_assert( rows ,       "Invalid Height" );
                    static_assert( isgettable , "Invalid Gettable" );
                    static_assert( isrestrict , "Invalid Restrict" );
                }
            };

            template < typename CLS ,
                       typename T ,
                       IndexType M ,
                       IndexType N ,
                       bool ISGETTABLE ,
                       bool ISRESTRICT ,
                       unsigned int MASK ,
                       typename Enabled = typename std::enable_if<
                            conditional_expression< CLS , T , M , N , ISGETTABLE , ISRESTRICT , MASK >::value >::type
                       >
            using expression_type = Expression< CLS >;
            template < typename CLS ,
                       IndexType M ,
                       IndexType N >
            using expression_size_type = expression_type< CLS , void , M , N , false , false , exm::size >;

            template < typename CLS ,
                       typename T ,
                       IndexType M ,
                       IndexType N >
            using expression_major_type = expression_type< CLS , T , M , N , false , false , exm::sizetype >;

            template < typename CLS , typename THIS >
            using expression_major_same = expression_major_type< CLS , ret_type< THIS > , Traits< THIS >::rows , Traits< THIS >::cols >;

            template < typename CLS ,
                       typename T ,
                       IndexType M ,
                       IndexType N ,
                       bool ISGETTABLE ,
                       bool ISRESTRICT ,
                       unsigned int MASK = exm::all ,
                       typename Enabled = typename std::enable_if<
                           conditional_expression< CLS , T , M , N , ISGETTABLE , ISRESTRICT , MASK >::value == false
                           >::type >
            using expression_error_type = Expression< CLS >;

            /*
             *  GetBy function branching
             *  - scalar
             *  - expression
             *      - is_gettable
             */
            //scalar
            template < typename T , typename = typename std::enable_if< std::is_arithmetic< T >::value >::type >
            constexpr _ehm_inline static T GetBy( const T a , IndexType i )
            {
                return a;
            }
            template < typename T , typename = typename std::enable_if< std::is_arithmetic< T >::value >::type >
            constexpr _ehm_inline static T GetBy(  const T a , IndexType x , IndexType y )
            {
                return a;
            }


            template < typename CLS >
            constexpr _ehm_inline static ret_type< CLS > GetBy( const Expression< CLS >& m , IndexType i )
            {
                return m[i];
            }
            template < typename CLS ,
                     typename = typename std::enable_if<
                         Traits< CLS >::is_gettable == false
                         >::type >
            constexpr _ehm_inline static ret_type< CLS > GetBy( const Expression< CLS >& m , IndexType x , IndexType y )
            {
                return m[ y + x * Traits< CLS >::rows ];
            }
            template < typename CLS , typename = void ,
                     typename = typename std::enable_if<
                         Traits< CLS >::is_gettable == true
                         >::type >
            constexpr _ehm_inline static ret_type< CLS > GetBy( const Expression< CLS >& m , IndexType x , IndexType y )
            {
                return m.Get( x , y );
            }

            template < typename CLS >
            constexpr _ehm_inline static ret_type< CLS >& GetByRef( Expression< CLS >& m , IndexType i )
            {
                return m[i];
            }
            template < typename CLS , typename = void ,
                     typename = typename std::enable_if<
                         Traits< CLS >::is_gettable == false
                         >::type >
            constexpr _ehm_inline static ret_type< CLS >& GetByRef( Expression< CLS >& m , IndexType x , IndexType y )
            {
                return m[ y + x * Traits< CLS >::rows ];
            }
            template < typename CLS , typename = void , typename = void ,
                     typename = typename std::enable_if<
                         Traits< CLS >::is_gettable == true
                         >::type >
            constexpr _ehm_inline static ret_type< CLS >& GetByRef( Expression< CLS >& m , IndexType x , IndexType y )
            {
                return m.Get( x , y );
            }



            template < typename CLS , int ACCESS_FACTOR >
            struct ShouldMakeTemp
            {
                /*
                 * for Expression that has many access to same index; so make it to calculate many times
                 * decide whether to make temp variable to cache the result
                 * expression's size & operations are counting with it; ( how many times they access? )
                 * Access_Factor would be how many times you will access same index;
                 */
                constexpr static const bool value = //Traits< CLS >::operations * ACCESS_FACTOR > 8;
                    Traits< CLS >::operations * (ACCESS_FACTOR-1) > Traits< CLS >::rows * Traits< CLS >::cols;

                // if should make temp object, returns matrix type
                // else returns expression-reference type
                typedef typename std::conditional<
                            value ,
                            typename std::add_const< mat_type< CLS > >::type ,
                            auto_creference< CLS >
                        >::type type;

                constexpr const static bool is_gettable = value ? false : Traits< CLS >::is_gettable;
                constexpr const static bool is_restrict = value ? true : Traits< CLS >::is_restrict;

                constexpr const static int operations = value ? 0 : Traits< CLS >::operations;
            };
        };  // namespace Expression
    };  // namespace Matrix
};  // namespace EH
