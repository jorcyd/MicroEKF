#include "ekf_math.h"

/* References to bit-hack based routines :
 * https://github.com/hcs0/Hackers-Delight/blob/master/rsqrt.c.txt
 * https://bits.stephan-brumme.com/squareRoot.html 
 * https://github.com/hcs0/Hackers-Delight/blob/master/asqrt.c.txt
 * https://bits.stephan-brumme.com/invSquareRoot.html 
 * https://stackoverflow.com/questions/31117497/fastest-integer-square-root-in-the-least-amount-of-instructions
 * https://stackoverflow.com/questions/31031223/fast-approximate-float-division	
 * https://stackoverflow.com/questions/4930307/fastest-way-to-get-the-integer-part-of-sqrtn
 * https://stackoverflow.com/questions/23474796/is-there-a-fast-fabsf-replacement-for-float-in-c
 */

#ifdef EMBEDDED	//if set to run on embedded targets with basic fp support.

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

inline float fast_inv(float x) {
	union { float f; int i; } v;
	float w;//, sx;

	// sx = (x < 0) ? -1:1;
	// x = sx * x;

	v.i = (int)(0x7EF127EA - *(unsigned int *)&x);
	w = x * v.f;

	// Efficient Iterative Approximation Improvement in horner polynomial form.
	//v.f = v.f * (2 - w);     					// Onte iteration, Err = -3.36e-3 * 2^(-flr(log2(x)))
	v.f = v.f * ( 4 + w * (-6 + w * (4 - w)));  // Two iterations, Err = -1.13e-5 * 2^(-flr(log2(x)))
	//v.f = v.f * (8 + w * (-28 + w * (56 + w * (-70 + w *(56 + w * (-28 + w * (8 - w)))))));  // Three Iteration, Err = +-6.8e-8 *  2^(-flr(log2(x)))
	return v.f;//* sx;
}

inline float fast_sqrtf(float x0) {
	union {int ix; float x;} c;

	c.x = x0;                      			// x can be viewed as int.
	c.ix = 0x1fbb67a8 + (c.ix >> 1); 		// Initial guess.
	c.x = 0.5f*(c.x + x0*fast_inv(c.x));    // Newton step.
	c.x = 0.5f*(c.x + x0*fast_inv(c.x));  // 2nd iter (required)
	return c.x;
}

inline float fast_fabsf(float i){
	(*(unsigned int *)&i) &= 0x7fffffff;
    return i;
}

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

inline float fast_fabsf(float x){
	return(fabsf(x));
}

#endif

unsigned int fast_isqrt(unsigned int a) {
	unsigned int rem = 0;
	unsigned int root = 0;
	int i;

	for (i = 0; i < 16; i++) {
		root <<= 1;
		rem <<= 2;
		rem += a >> 30;
		a <<= 2;

		if (root < rem) {
			root++;
			rem -= root;
			root++;
		}
	}

	return (root >> 1);
}
