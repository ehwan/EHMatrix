#EHMatrix
##expression template Matrix-Vector library
> is aiming to be an glm-like optimized matrix-vector library using [expression template](https://en.wikipedia.org/wiki/Expression_templates)

> written in C++14

> _**sample.cpp**_ contains **Sample code, some functions useful**






##Interface

###Include Header
    put EHMatrix folder and EHMatrix.h header file together.
    #include "EHMatrix.h"
    will automatically include all features



### Main namespace
    namespace EH::Matrix;




###Major Structures
####struct EH::Matrix::Matrix< TYPE , M , N >
* M the number of rows , N the number of columns
* the M by N matrix type filled with M*N TYPE variables
* the major matrix type you will directly use
* derived from **_Expression_** type below

####*Column Matrix* to use *vector type*
* Matrix< TYPE , N , 1 > to use N vector type
* aliased name Vector< TYPE , N >
* several specialized functions

###*Square Matrix*
* Matrix< TYPE , N , N >
* several specialized functions

####struct EH::Matrix::Expression::Expression
* __EVERYTHING__ working with matrix type above, __will return Expression type__
    - _eg.) mat1 + mat2 will return Expression< mat1 + mat2 >_ , not computed mat3
    - Matrix itself is also Expression
* Expression type is **not computed** itself; it's just containing the info of the expression.
* the expression will be computed( which actually takes your runtime resources ) if you **assign to Matrix** type
