#include "ekf_math.h"

#ifdef EMBEDDED	//if set to run on embedded targets with light fp support.

//https://github.com/hcs0/Hackers-Delight/blob/master/rsqrt.c.txt
inline float fast_rsqrtf(float x){	//1/sqrt(x)
	// for Newton iteration
	float xHalf = 0.5f*x;
	// approximation with empirically found "magic number"
	unsigned int *i = (unsigned int*) &x;
	//*i = 0x5F375A86 - (*i>>1);
 	*i = 0x5f37599e - (*i>>1); 	//The constant 0x5f37599e makes the relative error range from 0 to -0.00000463.
	// one Newton iteration, repeating further improves precision
	return x * (1.5f - xHalf*x*x);
	//return x;
}

// https://stackoverflow.com/questions/31031223/fast-approximate-float-division
inline float fast_inv(float x) {
	union { float f; int i; } v;
	float w;//, sx;

	// sx = (x < 0) ? -1:1;
	// x = sx * x;

	v.i = (int)(0x7EF127EA - *(unsigned int *)&x);
	w = x * v.f;

	// Efficient Iterative Approximation Improvement in horner polynomial form.
	v.f = v.f * (2 - w);     // Single iteration, Err = -3.36e-3 * 2^(-flr(log2(x)))
	// v.f = v.f * ( 4 + w * (-6 + w * (4 - w)));  // Second iteration, Err = -1.13e-5 * 2^(-flr(log2(x)))
	// v.f = v.f * (8 + w * (-28 + w * (56 + w * (-70 + w *(56 + w * (-28 + w * (8 - w)))))));  // Third Iteration, Err = +-6.8e-8 *  2^(-flr(log2(x)))
	return v.f;//* sx;
}

/* 	Quick and dirty fast-float-sqrt 
	https://bits.stephan-brumme.com/squareRoot.html 
	https://github.com/hcs0/Hackers-Delight/blob/master/asqrt.c.txt
	https://bits.stephan-brumme.com/invSquareRoot.html 
	https://stackoverflow.com/questions/31117497/fastest-integer-square-root-in-the-least-amount-of-instructions
	*/
inline float fast_sqrtf(float x0) {
	union {int ix; float x;} c;

	c.x = x0;                      		// x can be viewed as int.
	c.ix = 0x1fbb67a8 + (c.ix >> 1); 	// Initial guess.
	c.x = 0.5f*(c.x + x0/c.x);         	// Newton step.
	c.x = 0.5f*(c.x + x0/c.x);  		// 2nd iter (otherwise the estimates would suffer)
	return c.x;
}

#if 0
static inline float fast_sqrtf(float x) {
	unsigned int i = *(unsigned int*) &x;
	// adjust bias
	i  += 127 << 23;
	// approximation of square root
	i >>= 1;
	return *(float*) &i;
}

static inline float fast_sqrtf(float x) {
	unsigned int i = *(unsigned int*) &x;
	i = 0x1fbb4f2e + (i >> 1);
	return *(float*) &i;
}
#endif

#else

#include <math.h>
inline float fast_rsqrtf(float x){
	return(1/sqrtf(x));
}

inline float fast_inv(float x) {
	return(1/x);
}

inline float fast_sqrtf(float x0) {
	return(sqrtf(x0));
}

#endif