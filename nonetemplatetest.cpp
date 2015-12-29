#include <iostream>
#include "EHMatrix.h"

struct testvec
{
    float s[ 2 ];

    testvec(){}
    testvec( float x )
    {
        std::fill( s , s+2 , x );
    }
};

testvec operator * ( const testvec& a , const testvec& b )
{
    testvec ret;
    for( unsigned int i=0; i<2; ++i )
    {
        ret.s[i] = a.s[i] * b.s[i];
    }
    return ret;
}
testvec operator + ( const testvec& a , const testvec& b )
{
    testvec ret;
    for( unsigned int i=0; i<2; ++i )
    {
        ret.s[i] = a.s[i] + b.s[i];
    }
    return ret;
}

int main()
{
}

/*
 *
	movss	(%rdx), %xmm4
	movss	(%rsi), %xmm0
	movaps	%xmm4, %xmm6
	movaps	%xmm0, %xmm7
	movss	(%rdi), %xmm2
	mulss	%xmm2, %xmm7
	movss	4(%rdx), %xmm3
	mulss	%xmm2, %xmm6
	movss	4(%rsi), %xmm5
	movss	4(%rdi), %xmm1
	addss	%xmm7, %xmm6
	addss	%xmm6, %xmm2
	addss	%xmm0, %xmm2
	addss	%xmm4, %xmm2
	mulss	%xmm2, %xmm0
	movaps	%xmm5, %xmm2
	mulss	%xmm1, %xmm2
	addss	%xmm0, %xmm4
	movaps	%xmm3, %xmm0
	mulss	%xmm1, %xmm0
	addss	%xmm2, %xmm0
	addss	%xmm0, %xmm1
	addss	%xmm5, %xmm1
	addss	%xmm3, %xmm1
	mulss	%xmm5, %xmm1
	addss	%xmm1, %xmm3
	addss	%xmm3, %xmm4
	movaps	%xmm4, %xmm0
	ret
*/
float test1( const testvec& v1 , const testvec& v2 , const testvec& v3 )
{
    const testvec vv = ( v1*v2 + v3*v1 + v1 + v2 + v3 ) * v2 + v3;

    return vv.s[0] + vv.s[1];
}

typedef EH::Matrix::Vector< float , 2 > myvec2;

/*
 *
	movss	(%rdi,%rax), %xmm0
	movaps	%xmm0, %xmm4
	movaps	%xmm0, %xmm3
	movss	(%rsi,%rax), %xmm2
	movss	(%rdx,%rax), %xmm1
	mulss	%xmm2, %xmm4
	mulss	%xmm1, %xmm3
	addss	%xmm4, %xmm3
	addss	%xmm3, %xmm0
	addss	%xmm2, %xmm0
	addss	%xmm1, %xmm0
	mulss	%xmm2, %xmm0
	addss	%xmm0, %xmm1
	movss	%xmm1, -24(%rsp,%rax)
	addq	$4, %rax
	cmpq	$8, %rax
	jne	.L6
	movss	-24(%rsp), %xmm0
	addss	-20(%rsp), %xmm0
	ret
*/
float test2( const myvec2& v1 , const myvec2& v2 , const myvec2& v3 )
{
    const myvec2 vv = ( v1*v2 + v3*v1 + v1 + v2 + v3 ) * v2 + v3;

    return vv[0] + vv[1];
}
