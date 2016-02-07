#include <iostream>
#include <typeinfo>
#include "EHMatrix.h"

// unsigned int
using EH::Matrix::IndexType;
// general matrix type;
template < typename MEMBER_TYPE , IndexType ROWS , IndexType COLS >
using mat = EH::Matrix::Matrix< MEMBER_TYPE , ROWS , COLS >;

// square matrix type
template < typename MEMBER_TYPE , IndexType SIZE >
using mats = mat< MEMBER_TYPE , SIZE , SIZE >;

// vector type
template < typename MEMBER_TYPE , IndexType SIZE >
using vec = mat< MEMBER_TYPE , SIZE  , 1 >;
//can be used as EH::Matrix::Vector< TYPE , SIZE >

int main()
{
    EH::Matrix::Util::EyeMatrix( { std::atan( 1 ) , 0 , 0 } ).Log();

    return 0;
    {
        vec< float , 2 > v1( 0 , 1 );
        vec< float , 2 > v2( 1 , 2 );
        vec< float , 2 > v3 = ( v1 * 4 ) * v2;

        v3.Log();

        //v3.Log();
    }

    // 2 by 3 matrix
    // column-major,
    // ( 0 , 2 , 3 )
    // ( 1 , 3 , 4 )
    mat< float , 2 , 3 > m1 =
    {
        0 , 1 ,
        2 , 3 ,
        3 , 4
    };

    // access by operator []
    // OR
    // Get( x , y ) function;
    m1[3];   // will return 3
    // same as
    m1.Get( 1 , 1 );

    // iterator functions
    std::fill( m1.begin() , m1.end() , 0 );

    // data() to get pointer
    m1.data();

    // vector3 type;
    vec< float , 3 > v1 = { 0 , 1 , 2 };
    v1.Log();
    vec< float , 3 > v3 = { 0 , 2 , 1 };
    v3.Log();

    // aggressive assign;
    vec< float , 5 > v5( v1 , 5.0f , 6.0f );
    v5.Log();
    v5.FillAggressive( 7.0f , v3 , 8.0f );
    v5.Log();
    v5.Convert< unsigned int >();

    /*
     * for Matrix Aggressive assign;
     *
     * it will recursively fill given matrix column-ward
     *  - must fill bottom first to go right
     *
     * example )
     * m1 = 3 x 3
     * v1 = 3 x 1
     * v2 = 1 x 3
     * s1 = scalar
     *
     * m2 = 4 x 4
     *
     * m2( m1 , v2 , v1 , s1 )
     * m2 will be assigned as
     * ( m1 m1 m1 v1 )
     * ( m1 m1 m1 v1 )
     * ( m1 m1 m1 v1 )
     * ( v2 v2 v2 s1 )
     *
     *
     */

    {
        mats< float , 3 > m1 = 3;
        vec< float , 3 > v1 = 4;
        mat< float , 1 , 3 > v2 = { 0 , 1 , 2 };
        mats< float , 4 > m4( m1 , v2 , v1 , 3 );
        m4.Log();
    }



    // EXP::has_same_root( matrix_type )
    // returns true if Expression
    // contains expression of matrix_type
    auto exp = (( v1 * v3 + 2 ) / v3);
    //std::cout << exp.has_same_root( v1 );   // should returns true for v1 & v3

    // general matrix multiply;
    // ( 2x3 multipy 3x2 = 2x2 )
    // the Transpose() function.
    mat< float , 2 , 2 > m2 = m1 * m1.Transpose();

    // there is no 1 by 1 matrix; will return scalar
    // dot product; same as LengthSquared()
    float sq = v1.Transpose() * v1;


    // general matrix-vector multiply
    vec< float , 2 > v2 = m1 * v1;


    v1.x; v1.r;   // aliased-accessing , only for vector
    v1.y; v1.g;   // aliased-accessing , only for vector

    // vector-specialized functions
    v2.Length();
    v2.LengthSquared();
    // return the original length;
    float length = v2.Normalize();


    // for vector type;
    // several special operations

    // vector-scalar add/sub
    // scalar will treated as vector filled with scalar;
    v1 = v1 + 2;        // v1 - vec2( 2 , 2 )
    v1 = 2 + v1;        // vec2( 2 , 2 ) + v1
    v1 = 1 - v3;        // vec2( 1 , 1 ) - v2
    v3 = v3 - 2;        // v2 - vec2( 2 , 2 )

    // also a assign-operator
    v1 += 3;
    // so on...

    // vector-vector multiply
    // apply operations to each component
    v1 = v1 * v3;       // v1[0] = v1[0] * v2[0] , v1[1] = v1[1] * v2[1]
    v1 = v1 / v3;

    // proxy expression Diagonal();
    m2.Diagonal() += vec< float , 2 >{ 0 , 1 };
    //( x += 0 ,  y       )
    //( z      ,  z += 1  )

    // proxy
    // Row & Column
    m2.Column( 1 )  *= vec< float , 2 >{ 2 , 2 };
    m2.Row( 0 )     *= 2;

    // proxy
    //                          M   N    X   Y
    m2 = m2.template SubMatrix< 2 , 2 >( 0 , 0 );

    // some functions implemented
    m2 = EH::Matrix::floor( m2 );
    m2 = EH::Matrix::ceil( m2 );
    m2 = EH::Matrix::round( m2 );


    // compare operator;
    // matrix - matrix compare;
    // matrix - scalar compare;
    //
    // will return bool vector( expression )
    // all() , any() , none() functions
    // can be converted into std::bitset , bool
    // operator bool() will call all() function

    auto bits1 = m1 == 2.0f;                   // bits1 is expression type
    auto bits2 = 5 >= m2;                   // bits2 is expression type
    auto bits3 = m2 == m2.Transpose();      // bits3 is expression type
    // implict-conversion
    std::bitset< 4 > stdbits = bits3;       // implict-conversion to std::bitset
    auto stdbits2 = bits3.bitset();

    mat< bool , 2 , 2 > boolmat = bits3;    // assign to bool-matrix

    bits1.all();                            // expression-specialized
    bits2.any();
    bits3.none();

    if( bits1 ) // same as -> if( (m1 == 2).all() )     // conversion to bool
    {
        // do something if all components of m1 is equals to 2;
    }

}
