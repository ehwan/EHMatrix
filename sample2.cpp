#include "EHMatrix.h"

template < typename T , EH::Matrix::IndexType M , EH::Matrix::IndexType N >
using mat = EH::Matrix::Matrix< T , M , N >;

template < typename T , EH::Matrix::IndexType N >
using vec = EH::Matrix::Vector< T , N >;

int main()
{
    using mat2 = mat< float , 2 , 2 >;

    mat2 m1( 0 , 1 , 2 , 3 );

    m1.Log();

    m1.Column( 1 ) = 0;
    m1.Log();

    m1 = EH::Matrix::Expressions::make_unary( m1 ,
            []( auto a )
            {
                return a + 1;
            }
            );
    m1.Log();


    m1 += 2;

    m1.Log();

    using vec2 = vec< float , 2 >;
    vec2 v1 = { 0 , 1 };
    vec2 v2 = { 1 , 2 };

    m1 = m1 + v1;

    m1.Log();

    m1 = EH::Matrix::Expressions::make_index_filter< false , false , 2 , 2 >( m1 ,
            []( auto& x , auto& y )
            {
                x = 1;
                y = 0;
            }
        );

    m1.Log();
}
