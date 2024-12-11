/******************************************************************************
 *
 * Filename:    ieeehalfprecision.c
 * Programmer:  James Tursa
 * Version:     2.0
 * Date:        May 18, 2020
 * Copyright:   (c) 2009, 2020 by James Tursa, All Rights Reserved
 *
 *  This code uses the BSD License:
 *
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are 
 *  met:
 *
 *     * Redistributions of source code must retain the above copyright 
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright 
 *       notice, this list of conditions and the following disclaimer in 
 *       the documentation and/or other materials provided with the distribution
 *      
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * This file contains C code to convert between IEEE double, single, and half
 * precision floating point formats. The intended use is for standalone C code
 * that does not rely on MATLAB mex.h. The bit pattern for the half precision
 * floating point format is stored in a 16-bit unsigned int variable. The half
 * precision bit pattern definition is:
 *
 * 1 bit sign bit
 * 5 bits exponent, biased by 15
 * 10 bits mantissa, hidden leading bit, normalized to 1.0
 *
 * Special floating point bit patterns recognized and supported:
 *
 * All exponent bits zero:
 * - If all mantissa bits are zero, then number is zero (possibly signed)
 * - Otherwise, number is a denormalized bit pattern
 *
 * All exponent bits set to 1:
 * - If all mantissa bits are zero, then number is +Infinity or -Infinity
 * - Otherwise, number is NaN (Not a Number)
 *
 * Revision History:
 * 1.0  2009-Mar-03  Original Release
 * 2.0  2020-May-18  Updated for rounding modes, OpenMP, extra functions
 *
 * Functions included:
 *
 * singles2halfp: Converts singles to half precision
 * doubles2halfp: Converts doubles to half precision
 * halfp2singles: Converts half precision to singles
 * halfp2doubles: Converts half precision to doubles
 * halfp2isinf: Returns 1's or 0's if half precision input is inf
 * halfp2isnan: Returns 1's or 0's if half precision input is NaN
 * halfp2isnormal: Returns 1's or 0's if half precision input is normalized
 * halfp2isdenormal: Returns 1's or 0's if half precision input is denormalized
 * halfp2isinf_interleaved_complex: Same as halfp2inf but considers pairs of numbers
 * halfp2isnan_interleaved_complex: Same as halfp2nan but considers pairs of numbers
 * halfp2isnormal_interleaved_complex: Same as halfp2normal but considers pairs of numbers
 * halfp2isdenormal_interleaved_complex: Same as halfp2inf but considers pairs of numbers
 * halfp2eps: Returns the eps of the input (value of least significant mantissa bit)
 *
 ********************************************************************************/

/*

  Roundings to nearest

    Round to nearest, ties to even – rounds to the nearest value; if the number falls midway, it is rounded to the nearest value with an even least significant digit; this is the default for binary floating point and the recommended default for decimal.
    Round to nearest, ties away from zero – rounds to the nearest value; if the number falls midway, it is rounded to the nearest value above (for positive numbers) or below (for negative numbers); this is intended as an option for decimal floating point.

  Directed roundings

    Round toward 0 – directed rounding towards zero (also known as truncation).
    Round toward +inf – directed rounding towards positive infinity (also known as rounding up or ceiling).
    Round toward -inf – directed rounding towards negative infinity (also known as rounding down or floor).

  rounding mode    C name         MSVC name
  -----------------------------------------
  to nearest       FE_TONEAREST   _RC_NEAR
  toward zero      FE_TOWARDZERO  _RC_CHOP
  to +infinity     FE_UPWARD      _RC_UP
  to -infinity     FE_DOWNWARD    _RC_DOWN

  is_quiet = 1 Forces the most significant mantissa bit of NaN results to be set (i.e., quiet). (RECOMMENDED)
               Shifted payload bits are also carried over (those that don't shift off the right).
           = 0 Ambiguous result. Leaves shifted NaN mantissa bits intact if possible, but if the
               shift causes all resulting mantissa bits to be 0, then the quiet bit is set instead.

 */

// Includes -------------------------------------------------------------------

#include <fenv.h>
#include <string.h>
#include <stddef.h>

#include "mex.h"

#include "matlab_version.h"
#include "matlab_version.c"

/* OpenMP ------------------------------------------------------------- */

#ifdef _OPENMP
#include <omp.h>
#endif

// Macros ---------------------------------------------------------------------

#define  INT16_TYPE          short
#define UINT16_TYPE unsigned short
#define  INT32_TYPE          int
#define UINT32_TYPE unsigned int

#define HALFP_PINF ((UINT16_TYPE) 0x7C00u)  // +inf
#define HALFP_MINF ((UINT16_TYPE) 0xFC00u)  // -inf
#define HALFP_PNAN ((UINT16_TYPE) 0x7E00u)  // +nan (only is_quite bit set, no payload)
#define HALFP_MNAN ((UINT16_TYPE) 0xFE00u)  // -nan (only is_quite bit set, no payload)

/* Define our own values for rounding_mode if they aren't already defined */
#ifndef FE_TONEAREST
    #define FE_TONEAREST    0x0000
    #define FE_UPWARD       0x0100
    #define FE_DOWNWARD     0x0200
    #define FE_TOWARDZERO   0x0300
#endif
#define     FE_TONEARESTINF 0xFFFF  /* Round to nearest, ties away from zero (apparently no standard C name for this) */

// Prototypes -----------------------------------------------------------------

int singles2halfp(void *target, void *source, ptrdiff_t numel, int rounding_mode, int is_quiet);
int doubles2halfp(void *target, void *source, ptrdiff_t numel, int rounding_mode, int is_quiet);
int halfp2singles(void *target, void *source, ptrdiff_t numel);
int halfp2doubles(void *target, void *source, ptrdiff_t numel);
void halfp2isinf(void *target, void *source, ptrdiff_t numel);
void halfp2isnan(void *target, void *source, ptrdiff_t numel);
void halfp2isnormal(void *target, void *source, ptrdiff_t numel);
void halfp2isdenormal(void *target, void *source, ptrdiff_t numel);
void halfp2isinf_interleaved_complex(void *target, void *source, ptrdiff_t numel);
void halfp2isnan_interleaved_complex(void *target, void *source, ptrdiff_t numel);
void halfp2isnormal_interleaved_complex(void *target, void *source, ptrdiff_t numel);
void halfp2isdenormal_interleaved_complex(void *target, void *source, ptrdiff_t numel);
void halfp2eps(void *target, void *source, ptrdiff_t numel);

//-----------------------------------------------------------------------------
//
// Routine:  singles2halfp
//
// Input:  source = Address of 32-bit floating point data to convert
//         numel  = Number of values at that address to convert
//         rounding_mode = (see above)
//         is_quiet = (see above)
//
// Output: target = Address of 16-bit data to hold output (numel values)
//         return value = 0 if native floating point format is IEEE
//                      = 1 if native floating point format is not IEEE
//
// Programmer:  James Tursa
//
//-----------------------------------------------------------------------------

int singles2halfp(void *target, void *source, ptrdiff_t numel, int rounding_mode, int is_quiet)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) target; // Type pun output as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) source; // Type pun input as an unsigned 32-bit int
    UINT16_TYPE    hs, he, hm, hr, hq;
    UINT32_TYPE x, xs, xe, xm, xt, zm, zt, z1;
    ptrdiff_t i;
    int hes, N;
    static int next;  // Little Endian adjustment
    static int checkieee = 1;  // Flag to check for IEEE754, Endian, and word size
    double one = 1.0; // Used for checking IEEE754 floating point format
    UINT32_TYPE *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (UINT32_TYPE *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            return 1;  // Floating point bit pattern is not IEEE 754
        }
        if( sizeof(INT16_TYPE) != 2 || sizeof(INT32_TYPE) != 4 ) {
            return 1;  // short is not 16-bits, or long is not 32-bits.
        }
        checkieee = 0; // Everything checks out OK
    }
    
    if( source == NULL || target == NULL ) { // Nothing to convert (e.g., imag part of pure real)
        return 0;
    }
    
    // Make sure rounding mode is valid
    if( rounding_mode < 0 )
        rounding_mode = fegetround( );
    if( rounding_mode != FE_TONEAREST &&
        rounding_mode != FE_UPWARD &&
        rounding_mode != FE_DOWNWARD &&
        rounding_mode != FE_TOWARDZERO &&
        rounding_mode != FE_TONEARESTINF ) {
        rounding_mode = FE_TONEAREST;
    }
    hq = is_quiet ? (UINT16_TYPE) 0x0200u: (UINT16_TYPE) 0x0000u;  // Force NaN results to be quiet?
    
    // Loop through the elements
#pragma omp parallel for private(x,xs,xe,xm,xt,zm,zt,z1,hs,he,hm,hr,hes,N)
    for( i=0; i<numel; i++ ) {
        x = xp[i];
        if( (x & 0x7FFFFFFFu) == 0 ) {  // Signed zero
            hp[i] = (UINT16_TYPE) (x >> 16);  // Return the signed zero
        } else { // Not zero
            xs = x & 0x80000000u;  // Pick off sign bit
            xe = x & 0x7F800000u;  // Pick off exponent bits
            xm = x & 0x007FFFFFu;  // Pick off mantissa bits
            xt = x & 0x00001FFFu;  // Pick off trailing 13 mantissa bits beyond the shift (used for rounding normalized determination)
            if( xe == 0 ) {  // Denormal will underflow, return a signed zero or smallest denormal depending on rounding_mode
                if( (rounding_mode == FE_UPWARD   && xm && !xs) ||  // Mantissa bits are non-zero and sign bit is 0
                    (rounding_mode == FE_DOWNWARD && xm &&  xs) ) { // Mantissa bits are non-zero and sigh bit is 1
                    hp[i] = (UINT16_TYPE) (xs >> 16) | (UINT16_TYPE) 1u;  // Signed smallest denormal
                } else {
                    hp[i] = (UINT16_TYPE) (xs >> 16);  // Signed zero
                }
            } else if( xe == 0x7F800000u ) {  // Inf or NaN (all the exponent bits are set)
                if( xm == 0 ) { // If mantissa is zero ...
                    hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else {
                    hm = (UINT16_TYPE) (xm >> 13); // Shift mantissa over
                    if( hm ) { // If we still have some non-zero bits (payload) after the shift ...
                        hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7C00u | hq | hm); // Signed NaN, shifted mantissa bits set
                    } else {
                        hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7E00u); // Signed NaN, only 1st mantissa bit set (quiet)
                    }
                }
            } else { // Normalized number
                hs = (UINT16_TYPE) (xs >> 16); // Sign bit
                hes = ((int)(xe >> 23)) - 127 + 15; // Exponent unbias the single, then bias the halfp
                if( hes >= 0x1F ) {  // Overflow
                    hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else if( hes <= 0 ) {  // Underflow exponent, so halfp will be denormal
                    xm |= 0x00800000u;  // Add the hidden leading bit
                    N = (14 - hes);  // Number of bits to shift mantissa to get it into halfp word
                    hm = (N < 32) ? (UINT16_TYPE) (xm >> N) : (UINT16_TYPE) 0u; // Halfp mantissa
                    hr = (UINT16_TYPE) 0u; // Rounding bit, default to 0 for now (this will catch FE_TOWARDZERO and other cases)
                    if( N <= 24 ) {  // Mantissa bits have not shifted away from the end
                        zm = (0x00FFFFFFu >> N) << N;  // Halfp denormal mantissa bit mask
                        zt = 0x00FFFFFFu & ~zm;  // Halfp denormal trailing manntissa bits mask
                        z1 = (zt >> (N-1)) << (N-1);  // First bit of trailing bit mask
                        xt = xm & zt;  // Trailing mantissa bits
                        if( rounding_mode == FE_TONEAREST ) {
                            if( xt > z1 || xt == z1 && (hm & 1u) ) { // Trailing bits are more than tie, or tie and mantissa is currently odd
                                hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                            }
                        } else if( rounding_mode == FE_TONEARESTINF ) {
                            if( xt >= z1  ) { // Trailing bits are more than or equal to tie
                                hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                            }
                        } else if( (rounding_mode == FE_UPWARD   && xt && !xs) ||  // Trailing bits are non-zero and sign bit is 0
                                   (rounding_mode == FE_DOWNWARD && xt &&  xs) ) { // Trailing bits are non-zero and sign bit is 1
                            hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                        }                    
                    } else {  // Mantissa bits have shifted at least one bit beyond the end (ties not possible)
                        if( (rounding_mode == FE_UPWARD   && xm && !xs) ||  // Trailing bits are non-zero and sign bit is 0
                            (rounding_mode == FE_DOWNWARD && xm &&  xs) ) { // Trailing bits are non-zero and sign bit is 1
                            hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                        }                    
                    }
                    hp[i] = (hs | hm) + hr; // Combine sign bit and mantissa bits and rounding bit, biased exponent is zero
                } else {
                    he = (UINT16_TYPE) (hes << 10); // Exponent
                    hm = (UINT16_TYPE) (xm >> 13); // Mantissa
                    hr = (UINT16_TYPE) 0u; // Rounding bit, default to 0 for now
                    if( rounding_mode == FE_TONEAREST ) {
                        if( xt > 0x00001000u || xt == 0x00001000u && (hm & 1u) ) { // Trailing bits are more than tie, or tie and mantissa is currently odd
                            hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                        }
                    } else if( rounding_mode == FE_TONEARESTINF ) {
                        if( xt >= 0x00001000u  ) { // Trailing bits are more than or equal to tie
                            hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                        }
                    } else if( (rounding_mode == FE_UPWARD   && xt && !xs) ||  // Trailing bits are non-zero and sign bit is 0
                               (rounding_mode == FE_DOWNWARD && xt &&  xs) ) { // Trailing bits are non-zero and sign bit is 1
                        hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                    }
                    hp[i] = (hs | he | hm) + hr;  // Adding rounding bit might overflow into exp bits, but that is OK
                }
            }
        }
    }
    return 0;
}

//-----------------------------------------------------------------------------
//
// Routine:  doubles2halfp
//
// Input:  source = Address of 64-bit floating point data to convert
//         numel  = Number of values at that address to convert
//         rounding_mode = (see above)
//         is_quiet = (see above)
//
// Output: target = Address of 16-bit data to hold output (numel values)
//         return value = 0 if native floating point format is IEEE
//                      = 1 if native floating point format is not IEEE
//
// Programmer:  James Tursa
//
//-----------------------------------------------------------------------------

int doubles2halfp(void *target, void *source, ptrdiff_t numel, int rounding_mode, int is_quiet)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) target; // Type pun output as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) source; // Type pun input as an unsigned 32-bit int
    UINT16_TYPE    hs, he, hm, hq, hr;
    UINT32_TYPE x, xs, xe, xm, xn, xt, zm, zt, z1;
    int hes, N;
    ptrdiff_t i;
    static int next;  // Little Endian adjustment
    static int checkieee = 1;  // Flag to check for IEEE754, Endian, and word size
    double one = 1.0; // Used for checking IEEE754 floating point format
    UINT32_TYPE *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (UINT32_TYPE *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            return 1;  // Floating point bit pattern is not IEEE 754
        }
        if( sizeof(INT16_TYPE) != 2 || sizeof(INT32_TYPE) != 4 ) {
            return 1;  // short is not 16-bits, or long is not 32-bits.
        }
        checkieee = 0; // Everything checks out OK
    }
    
    if( source == NULL || target == NULL ) { // Nothing to convert (e.g., imag part of pure real)
        return 0;
    }
    
    hq = is_quiet ? (UINT16_TYPE) 0x0200u: (UINT16_TYPE) 0x0000u;  // Force NaN results to be quiet?
    
#pragma omp parallel for private(x,xs,xe,xm,xt,zm,zt,z1,hs,he,hm,hr,hes,N)
    for( i=0; i<numel; i++ ) {
        if( next ) { // Little Endian
            xn = xp[2*i];  // Lower mantissa
            x  = xp[2*i+1];  // Sign, exponent, upper mantissa
        } else { // Big Endian
            x  = xp[2*i];  // Sign, exponent, upper mantissa
            xn = xp[2*i+1];  // Lower mantissa
        }
        if( (x & 0x7FFFFFFFu) == 0 ) {  // Signed zero
            hp[i]= (UINT16_TYPE) (x >> 16);  // Return the signed zero
        } else { // Not zero
            xs = x & 0x80000000u;  // Pick off sign bit
            xe = x & 0x7FF00000u;  // Pick off exponent bits
            xm = x & 0x000FFFFFu;  // Pick off mantissa bits
            xt = x & 0x000003FFu;  // Pick off trailing 10 mantissa bits beyond the shift (used for rounding normalized determination)
            if( xe == 0 ) {  // Denormal will underflow, return a signed zero or signed smallest denormal depending on rounding_mode
                if( (rounding_mode == FE_UPWARD   && (xm || xn) && !xs) ||  // Mantissa bits are non-zero and sign bit is 0
                    (rounding_mode == FE_DOWNWARD && (xm || xn) &&  xs) ) { // Mantissa bits are non-zero and sigh bit is 1
                    hp[i] = (UINT16_TYPE) (xs >> 16) | (UINT16_TYPE) 1u;  // Signed smallest denormal
                } else {
                    hp[i] = (UINT16_TYPE) (xs >> 16);  // Signed zero
                }
            } else if( xe == 0x7FF00000u ) {  // Inf or NaN (all the exponent bits are set)
                if( xm == 0 && xn == 0 ) { // If mantissa is zero ...
                    hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else {
                    hm = (UINT16_TYPE) (xm >> 10); // Shift mantissa over
                    if( hm ) { // If we still have some non-zero bits (payload) after the shift ...
                        hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7C00u | hq | hm); // Signed NaN, shifted mantissa bits set
                    } else {
                        hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7E00u); // Signed NaN, only 1st mantissa bit set (quiet)
                    }
                }
            } else { // Normalized number
                hs = (UINT16_TYPE) (xs >> 16); // Sign bit
                hes = ((int)(xe >> 20)) - 1023 + 15; // Exponent unbias the double, then bias the halfp
                if( hes >= 0x1F ) {  // Overflow
                    hp[i] = (UINT16_TYPE) ((xs >> 16) | 0x7C00u); // Signed Inf
                } else if( hes <= 0 ) {  // Underflow exponent, so halfp will be denormal
                    xm |= 0x00100000u;  // Add the hidden leading bit
                    N = (11 - hes);  // Number of bits to shift mantissa to get it into halfp word
                    hm = (N < 32) ? (UINT16_TYPE) (xm >> N) : (UINT16_TYPE) 0u; // Halfp mantissa
                    hr = (UINT16_TYPE) 0u; // Rounding bit, default to 0 for now (this will catch FE_TOWARDZERO and other cases)
                    if( N <= 21 ) {  // Mantissa bits have not shifted away from the end
                        zm = (0x001FFFFFu >> N) << N;  // Halfp denormal mantissa bit mask
                        zt = 0x001FFFFFu & ~zm;  // Halfp denormal trailing manntissa bits mask
                        z1 = (zt >> (N-1)) << (N-1);  // First bit of trailing bit mask
                        xt = xm & zt;  // Trailing mantissa bits
                        if( rounding_mode == FE_TONEAREST ) {
                            if( (xt > z1 || xt == z1 && xn) || (xt == z1 && !xn) && (hm & 1u) ) { // Trailing bits are more than tie, or tie and mantissa is currently odd
                                hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                            }
                        } else if( rounding_mode == FE_TONEARESTINF ) {
                            if( xt >= z1 ) { // Trailing bits are more than or equal to tie
                                hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                            }
                        } else if( (rounding_mode == FE_UPWARD   && (xt || xn) && !xs) ||  // Trailing bits are non-zero and sign bit is 0
                                   (rounding_mode == FE_DOWNWARD && (xt || xn) &&  xs) ) { // Trailing bits are non-zero and sign bit is 1
                            hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                        }                    
                    } else {  // Mantissa bits have shifted at least one bit beyond the end (ties not possible)
                        if( (rounding_mode == FE_UPWARD   && (xm || xn) && !xs) ||  // Trailing bits are non-zero and sign bit is 0
                            (rounding_mode == FE_DOWNWARD && (xm || xn) &&  xs) ) { // Trailing bits are non-zero and sign bit is 1
                            hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                        }                    
                    }
                    hp[i] = (hs | hm) + hr; // Combine sign bit and mantissa bits and rounding bit, biased exponent is zero
                } else {
                    he = (UINT16_TYPE) (hes << 10); // Exponent
                    hm = (UINT16_TYPE) (xm >> 10); // Mantissa
                    hr = (UINT16_TYPE) 0u; // Rounding bit, default to 0 for now
                    if( rounding_mode == FE_TONEAREST ) {
                        if( (xt > 0x00000200u || xt == 0x00000200u && xn) || (xt == 0x00000200u && !xn) && (hm & 1u) ) { // Trailing bits are more than tie, or tie and mantissa is currently odd
                            hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                        }
                    } else if( rounding_mode == FE_TONEARESTINF ) {
                        if( xt >= 0x00000200u  ) { // Trailing bits are more than or equal to tie
                            hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                        }
                    } else if( (rounding_mode == FE_UPWARD   && (xt || xn) && !xs) ||  // Trailing bits are non-zero and sign bit is 0
                               (rounding_mode == FE_DOWNWARD && (xt || xn) &&  xs) ) { // Trailing bits are non-zero and sign bit is 1
                        hr = (UINT16_TYPE) 1u; // Rounding bit set to 1
                    }
                    hp[i] = (hs | he | hm) + hr;  // Adding rounding bit might overflow into exp bits, but that is OK
                }
            }
        }
    }
    return 0;
}

//-----------------------------------------------------------------------------
//
// Routine:  halfp2singles
//
// Input:  source = address of 16-bit data to convert
//         numel  = Number of values at that address to convert
//
// Output: target = Address of 32-bit floating point data to hold output (numel values)
//         return value = 0 if native floating point format is IEEE
//                      = 1 if native floating point format is not IEEE
//
// Programmer:  James Tursa
//
//-----------------------------------------------------------------------------

int halfp2singles(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) target; // Type pun output as an unsigned 32-bit int
    UINT16_TYPE h, hs, he, hm;
    UINT32_TYPE xs, xe, xm;
    INT32_TYPE xes;
    ptrdiff_t i;
    int e;
    static int next;  // Little Endian adjustment
    static int checkieee = 1;  // Flag to check for IEEE754, Endian, and word size
    double one = 1.0; // Used for checking IEEE754 floating point format
    UINT32_TYPE *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (UINT32_TYPE *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            return 1;  // Floating point bit pattern is not IEEE 754
        }
        if( sizeof(INT16_TYPE) != 2 || sizeof(INT32_TYPE) != 4 ) {
            return 1;  // short is not 16-bits, or long is not 32-bits.
        }
        checkieee = 0; // Everything checks out OK
    }
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return 0;
    
#pragma omp parallel for private(xs,xe,xm,h,hs,he,hm,xes,e)
    for( i=0; i<numel; i++ ) {
        h = hp[i];
        if( (h & 0x7FFFu) == 0 ) {  // Signed zero
            xp[i] = ((UINT32_TYPE) h) << 16;  // Return the signed zero
        } else { // Not zero
            hs = h & 0x8000u;  // Pick off sign bit
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( he == 0 ) {  // Denormal will convert to normalized
                e = -1; // The following loop figures out how much extra to adjust the exponent
                do {
                    e++;
                    hm <<= 1;
                } while( (hm & 0x0400u) == 0 ); // Shift until leading bit overflows into exponent bit
                xs = ((UINT32_TYPE) hs) << 16; // Sign bit
                xes = ((INT32_TYPE) (he >> 10)) - 15 + 127 - e; // Exponent unbias the halfp, then bias the single
                xe = (UINT32_TYPE) (xes << 23); // Exponent
                xm = ((UINT32_TYPE) (hm & 0x03FFu)) << 13; // Mantissa
                xp[i] = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            } else if( he == 0x7C00u ) {  // Inf or NaN (all the exponent bits are set)
                xp[i] = (((UINT32_TYPE) hs) << 16) | ((UINT32_TYPE) 0x7F800000u) | (((UINT32_TYPE) hm) << 13); // Signed Inf or NaN
            } else { // Normalized number
                xs = ((UINT32_TYPE) hs) << 16; // Sign bit
                xes = ((INT32_TYPE) (he >> 10)) - 15 + 127; // Exponent unbias the halfp, then bias the single
                xe = (UINT32_TYPE) (xes << 23); // Exponent
                xm = ((UINT32_TYPE) hm) << 13; // Mantissa
                xp[i] = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            }
        }
    }
    return 0;
}

//-----------------------------------------------------------------------------
//
// Routine:  halfp2singles
//
// Input:  source = address of 16-bit data to convert
//         numel  = Number of values at that address to convert
//
// Output: target = Address of 32-bit floating point data to hold output (numel values)
//         return value = 0 if native floating point format is IEEE
//                      = 1 if native floating point format is not IEEE
//
// Programmer:  James Tursa
//
//-----------------------------------------------------------------------------

int halfp2doubles(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) target; // Type pun output as an unsigned 32-bit int
    UINT16_TYPE h, hs, he, hm;
    UINT32_TYPE xs, xe, xm, xn;
    INT32_TYPE xes;
    ptrdiff_t i, j;
    int e;
    static int next;  // Little Endian adjustment
    static int checkieee = 1;  // Flag to check for IEEE754, Endian, and word size
    double one = 1.0; // Used for checking IEEE754 floating point format
    UINT32_TYPE *ip; // Used for checking IEEE754 floating point format
    
    if( checkieee ) { // 1st call, so check for IEEE754, Endian, and word size
        ip = (UINT32_TYPE *) &one;
        if( *ip ) { // If Big Endian, then no adjustment
            next = 0;
        } else { // If Little Endian, then adjustment will be necessary
            next = 1;
            ip++;
        }
        if( *ip != 0x3FF00000u ) { // Check for exact IEEE 754 bit pattern of 1.0
            return 1;  // Floating point bit pattern is not IEEE 754
        }
        if( sizeof(INT16_TYPE) != 2 || sizeof(INT32_TYPE) != 4 ) {
            return 1;  // short is not 16-bits, or long is not 32-bits.
        }
        checkieee = 0; // Everything checks out OK
    }
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return 0;
    
    xp += next;  // Little Endian adjustment if necessary
    
#pragma omp parallel for private(xs,xe,xm,h,hs,he,hm,xes,e,j)
    for( i=0; i<numel; i++ ) {
        j = 2*i;
        if( next ) {
            xp[j-1] = 0; // Set lower mantissa bits, Little Endian
        } else {
            xp[j+1] = 0; // Set lower mantissa bits, Big Endian
        }
        h = hp[i];
        if( (h & 0x7FFFu) == 0 ) {  // Signed zero
            xp[j] = ((UINT32_TYPE) h) << 16;  // Return the signed zero
        } else { // Not zero
            hs = h & 0x8000u;  // Pick off sign bit
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( he == 0 ) {  // Denormal will convert to normalized
                e = -1; // The following loop figures out how much extra to adjust the exponent
                do {
                    e++;
                    hm <<= 1;
                } while( (hm & 0x0400u) == 0 ); // Shift until leading bit overflows into exponent bit
                xs = ((UINT32_TYPE) hs) << 16; // Sign bit
                xes = ((INT32_TYPE) (he >> 10)) - 15 + 1023 - e; // Exponent unbias the halfp, then bias the double
                xe = (UINT32_TYPE) (xes << 20); // Exponent
                xm = ((UINT32_TYPE) (hm & 0x03FFu)) << 10; // Mantissa
                xp[j] = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            } else if( he == 0x7C00u ) {  // Inf or NaN (all the exponent bits are set)
                xp[j] = (((UINT32_TYPE) hs) << 16) | ((UINT32_TYPE) 0x7FF00000u) | (((UINT32_TYPE) hm) << 10); // Signed Inf or NaN
            } else {
                xs = ((UINT32_TYPE) hs) << 16; // Sign bit
                xes = ((INT32_TYPE) (he >> 10)) - 15 + 1023; // Exponent unbias the halfp, then bias the double
                xe = (UINT32_TYPE) (xes << 20); // Exponent
                xm = ((UINT32_TYPE) hm) << 10; // Mantissa
                xp[j] = (xs | xe | xm); // Combine sign bit, exponent bits, and mantissa bits
            }
        }
    }
    return 0;
}

//-----------------------------------------------------------------------------

void halfp2isinf(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT08_TYPE *xp = (UINT08_TYPE *) target; // Type pun output as an unsigned 8-bit int
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
#pragma omp parallel for
    for( i=0; i<numel; i++ ) {
        xp[i] = (hp[i] & 0x7FFFu) == 0x7C00u;  // All exponent bits set and mantissa is all 0's
    }
}

//-----------------------------------------------------------------------------

void halfp2isnan(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT08_TYPE *xp = (UINT08_TYPE *) target; // Type pun output as an unsigned 8-bit int
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
#pragma omp parallel for
    for( i=0; i<numel; i++ ) {
        xp[i] = (hp[i] & 0x7FFFu) >= 0x7C01u;  // All exponent bits set and mantissa is not all 0's
    }
}

//-----------------------------------------------------------------------------

void halfp2isnormal(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT08_TYPE *xp = (UINT08_TYPE *) target; // Type pun output as an unsigned 8-bit int
    UINT16_TYPE h, he, hm;
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
#pragma omp parallel for private(h,he,hm)
    for( i=0; i<numel; i++ ) {
        h = hp[i];
        he = h & 0x7C00u;  // Pick off exponent bits
        hm = h & 0x03FFu;  // Pick off mantissa bits
        if( he && he < 0x7C00u || !he && !hm ) {  // Exponent bits non-zero or mantissa is zero
            xp[i] = 1u;
        }
    }
}

//-----------------------------------------------------------------------------

void halfp2isdenormal(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT08_TYPE *xp = (UINT08_TYPE *) target; // Type pun output as an unsigned 8-bit int
    UINT16_TYPE h, he, hm;
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
#pragma omp parallel for private(h,he,hm)
    for( i=0; i<numel; i++ ) {
        h = hp[i];
        he = h & 0x7C00u;  // Pick off exponent bits
        hm = h & 0x03FFu;  // Pick off mantissa bits
        if( !he && hm ) {  // All exponent bits 0's and mantissa is non-zero
            xp[i] = 1u;
        }
    }
}

//-----------------------------------------------------------------------------

void halfp2isinf_interleaved_complex(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT08_TYPE *xp = (UINT08_TYPE *) target; // Type pun output as an unsigned 8-bit int
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
#pragma omp parallel
    for( i=0; i<numel; i++ ) {
        xp[i] = (hp[2*i] & 0x7FFFu) == 0x7C00u || (hp[2*i+1] & 0x7FFFu) == 0x7C00u;
    }
}

//-----------------------------------------------------------------------------

void halfp2isnan_interleaved_complex(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT08_TYPE *xp = (UINT08_TYPE *) target; // Type pun output as an unsigned 8-bit int
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
#pragma omp parallel for
    for( i=0; i<numel; i++ ) {
        xp[i] = (hp[2*i] & 0x7FFFu) >= 0x7C01u || (hp[2*i+1] & 0x7FFFu) >= 0x7C01u;
    }
}

//-----------------------------------------------------------------------------

void halfp2isnormal_interleaved_complex(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT08_TYPE *xp = (UINT08_TYPE *) target; // Type pun output as an unsigned 8-bit int
    UINT16_TYPE h, he, hm;
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
#pragma omp parallel for private(h,he,hm)
    for( i=0; i<numel; i++ ) {
        h = hp[2*i];
        he = h & 0x7C00u;  // Pick off exponent bits
        hm = h & 0x03FFu;  // Pick off mantissa bits
        if( he && he < 0x7C00u || !he && !hm ) {  // Exponent bits non-zero or mantissa is zero
            xp[i] = 1u;
        } else {
            h = hp[2*i+1];
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( he && he < 0x7C00u || !he && !hm ) {  // Exponent bits non-zero or mantissa is zero
                xp[i] = 1u;
            }
        }
    }
}

//-----------------------------------------------------------------------------

void halfp2isdenormal_interleaved_complex(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT08_TYPE *xp = (UINT08_TYPE *) target; // Type pun output as an unsigned 8-bit int
    UINT16_TYPE h, he, hm;
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
#pragma omp parallel for private(h,he,hm)
    for( i=0; i<numel; i++ ) {
        h = hp[2*i];
        he = h & 0x7C00u;  // Pick off exponent bits
        hm = h & 0x03FFu;  // Pick off mantissa bits
        if( !he && hm ) {  // All exponent bits 0's and mantissa is non-zero
            xp[i] = 1u;
        } else {
            h = hp[2*i+1];
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( !he && hm ) {  // All exponent bits 0's and mantissa is non-zero
                xp[i] = 1u;
            }
        }
    }
}

//-----------------------------------------------------------------------------

void halfp2eps(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT16_TYPE *xp = (UINT16_TYPE *) target; // Type pun output as an unsigned 32-bit int
    UINT16_TYPE h, he, hm;
    ptrdiff_t i;
    int e;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
    
#pragma omp parallel for private(h,he,hm,e)
    for( i=0; i<numel; i++ ) {
        h = hp[i];
        if( (h & 0x7FFFu) == 0 ) {  // Signed zero
            xp[i] = 0x0001u;  // Return the smallest denormal
        } else { // Not zero
            he = h & 0x7C00u;  // Pick off exponent bits
            hm = h & 0x03FFu;  // Pick off mantissa bits
            if( he == 0 ) {  // Denormal just returns smallest denormal
                xp[i] = 0x0001u;  // Return the smallest denormal
            } else if( he == 0x7C00u ) {  // Inf or NaN (all the exponent bits are set)
                xp[i] = 0x7E00u; // Returns +NaN
            } else { // Normalized number
                e = (int)(he >> 10) - 10;
                if( e <= 0 ) {  // Result will be denormal
                    xp[i] = (UINT16_TYPE) (0x0400u >> (1-e)); // Hidden leading bit shifted into denormal
                } else {  // Result will be normalized
                    xp[i] = (UINT16_TYPE) (e << 10);
                }
            }
        }
    }
}
