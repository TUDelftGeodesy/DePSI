/******************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Function:    halfprecision
 * Filename:    halfprecision.c
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
 * halfprecision converts the input argument to/from a half precision floating
 * point bit pattern corresponding to IEEE 754r. The bit pattern is stored in a
 * uint16 class variable. Please note that halfprecision is *not* a class. That
 * is, you cannot do any arithmetic with the half precision bit patterns.
 * halfprecision is simply a function that converts the IEEE 754r half precision
 * bit pattern to/from other numeric MATLAB variables. You can, however, take
 * the half precision bit patterns, convert them to single or double, do the
 * operation, and then convert the result back manually.
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
 * Building:
 *
 * halfprecision requires that a mex routine be built (one time only). This
 * process is typically self-building the first time you call the function
 * as long as you have the files halfprecision.m and halfprecision.c in the
 * same directory somewhere on the MATLAB path. If you need to manually build
 * the mex function, here are the commands:
 *
 * On older versions of MATLAB:
 * >> mex -setup
 *   (then follow instructions to select a C / C++ compiler of your choice)
 * >> mex halfprecision.c
 *
 * On newer versions of MATLAB R2017b or earlier:
 * >> mex halfprecision.c
 *
 * On later versions of MATLAB R2018a or later:
 * >> mex halfprecision.c -R2018a
 *
 * If you have an OPENMP compliant compiler, you should compile with the
 * appropriate compiler flags so that halfprecision can multi-thread.
 *
 * E.g., for Windows MSVC:
 * >> mex('COMPFLAGS="$COMPFLAGS /openmp"','halfprecision.c')
 * or
 * >> mex('COMPFLAGS="$COMPFLAGS /openmp"','halfprecision.c','-R2018a')
 *
 * E.g., for Windows MinGW64:
 * >> mex('CFLAGS="$CFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"','halfprecision.c')
 * or
 * >> mex('CFLAGS="$CFLAGS -fopenmp"','LDFLAGS="$LDFLAGS -fopenmp"','halfprecision.c','-R2018a')
 *
 * Syntax
 *
 * B = halfprecision(A)
 * C = halfprecision(B,classname or function)
 * L = halfprecision(directive)
 *     halfprecision(B,'disp')
 *
 * Description
 *
 * A = a MATLAB numeric array, char array, or logical array.
 *
 * B = the variable A converted into half precision floating point bit pattern.
 *     The bit pattern will be returned as a uint16 class variable. The values
 *     displayed are simply the bit pattern interpreted as if it were an unsigned
 *     16-bit integer. To see the halfprecision values, use the 'disp' option, which
 *     simply converts the bit patterns into a single class and then displays them.
 *
 * C = the half precision floating point bit pattern in B converted into class S.
 *     B must be a uint16 or int16 class variable.
 *
 * classname = char string naming the desired class (e.g., 'single', 'int32', etc.)
 *
 * function = char string giving one of the following functions:
 *            'isinf' = returns a logical variable, true where B is inf
 *            'isnan' = returns a logical variable, true where B is nan
 *            'isnormal' = returns a logical variable, true where B is normalized
 *            'isdenormal' = returns a logical variable, true where B is denormalized
 *            'eps' = returns eps of the half precision values
 *
 * directive = char string giving one of the following directives:
 *             'openmp' = returns a logical variable, true when compiled with OpenMP
 *             'realmax' = returns max half precision value
 *             'realmin' = returns min half precision normalized value
 *             'realmindenormal' = returns min half precision denormalized value
 *             'version' = returns a string with compilation memory model
 *
 *     'disp' = The floating point bit values are simply displayed.
 *
 * Examples
 * 
 * >> a = [-inf -1e30 -1.2 NaN 1.2 1e30 inf]
 * a =
 * 1.0e+030 *
 *     -Inf   -1.0000   -0.0000       NaN    0.0000    1.0000       Inf
 *
 * >> b = halfprecision(a)
 * b =
 * 64512  64512  48333  65024  15565  31744  31744
 *
 * >> halfprecision(b,'disp')
 *     -Inf      -Inf   -1.2002       NaN    1.2002       Inf       Inf
 *
 * >> halfprecision(b,'double')
 * ans =
 *     -Inf      -Inf   -1.2002       NaN    1.2002       Inf       Inf
 *
 * >> 2^(-24)
 * ans =
 * 5.9605e-008
 *
 * >> halfprecision(ans)
 * ans =
 *     1
 *
 * >> halfprecision(ans,'disp')
 * 5.9605e-008
 *
 * >> 2^(-25)
 * ans =
 * 2.9802e-008
 *
 * >> halfprecision(ans)
 * ans =
 *     1
 *
 * >> halfprecision(ans,'disp')
 * 5.9605e-008
 *
 * >> 2^(-26)
 * ans =
 *  1.4901e-008
 *
 * >> halfprecision(ans)
 * ans =
 *     0
 *
 * >> halfprecision(ans,'disp')
 *    0
 *
 * Note that the special cases of -Inf, +Inf, and NaN are handled correctly.
 * Also, note that the -1e30 and 1e30 values overflow the half precision format
 * and are converted into half precision -Inf and +Inf values, and stay that
 * way when they are converted back into doubles.
 *
 * For the denormalized cases, note that 2^(-24) is the smallest number that can
 * be represented in half precision exactly. 2^(-25) will convert to 2^(-24)
 * because of the rounding algorithm used, and 2^(-26) is too small and underflows
 * to zero.
 *
 * Revision History:
 * 1.0  2009-Mar-03  Original Release
 * 2.0  2020-May-16  Updated for rounding modes, directives, R2018a, R2019b
 *
 ********************************************************************************/

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

#define UINT08_TYPE unsigned char
#define  INT16_TYPE          short
#define UINT16_TYPE unsigned short
#define  INT32_TYPE          int
#define UINT32_TYPE unsigned int

#define HALFP_PINF ((UINT16_TYPE) 0x7C00u)  // +inf
#define HALFP_MINF ((UINT16_TYPE) 0xFC00u)  // -inf
#define HALFP_PNAN ((UINT16_TYPE) 0x7E00u)  // +nan (only is_quite bit set, no payload)
#define HALFP_MNAN ((UINT16_TYPE) 0xFE00u)  // -nan (only is_quite bit set, no payload)

// Global ---------------------------------------------------------------------

int next;  // Little Endian adjustment
int checkieee = 1;  // Flag to check for IEEE754, Endian, and word size
int warn = 1;  // Warning message if complex is used in R2018a+ but compiled in R2017b-

// Prototypes -----------------------------------------------------------------

void singles2halfp(void *target, void *source, ptrdiff_t numel, int rounding_mode, int is_quiet);
void doubles2halfp(void *target, void *source, ptrdiff_t numel, int rounding_mode, int is_quiet);

void halfp2singles(void *target, void *source, ptrdiff_t numel);
void halfp2doubles(void *target, void *source, ptrdiff_t numel);

void halfp2isinf(void *target, void *source, ptrdiff_t numel);
void halfp2isnan(void *target, void *source, ptrdiff_t numel);
void halfp2isnormal(void *target, void *source, ptrdiff_t numel);
void halfp2isdenormal(void *target, void *source, ptrdiff_t numel);

void halfp2isinf_interleaved_complex(void *target, void *source, ptrdiff_t numel);
void halfp2isnan_interleaved_complex(void *target, void *source, ptrdiff_t numel);
void halfp2isnormal_interleaved_complex(void *target, void *source, ptrdiff_t numel);
void halfp2isdenormal_interleaved_complex(void *target, void *source, ptrdiff_t numel);
void halfp2eps(void *target, void *source, ptrdiff_t numel);

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

/* Define our own values for rounding_mode if they aren't already defined */
#ifndef FE_TONEAREST
    #define FE_TONEAREST    0x0000
    #define FE_UPWARD       0x0100
    #define FE_DOWNWARD     0x0200
    #define FE_TOWARDZERO   0x0300
#endif
#define     FE_TONEARESTINF 0xFFFF  /* Round to nearest, ties away from zero (apparently no standard C name for this) */

// Gateway Function -----------------------------------------------------------

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    mwSize ndim; // Number of dimensions of input
    mwSize numel; // Number of elements of input
    const mwSize *dims; // Pointer to dimensions array
    mxClassID classid; // Class id of input or desired output
    mxComplexity complexity; // Complexity of input
    mxArray *rhs[2], *lhs[1], *currentformat[1]; // Used for callbacks into MATLAB
    char *classname; // Class name of desired output
    int disp = 0; // Display flag
    int isinf = 0; // Is inf flag
    int isnan = 0; // Is NaN flag
    int isnormal = 0; // Is normal flag
    int isdenormal = 0; // Is denormal flag
    int eps = 0; // Eps flag
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
            mexErrMsgTxt("Floating point bit pattern is not IEEE 754");
        }
        if( sizeof(INT16_TYPE) != 2 || sizeof(INT32_TYPE) != 4 ) {
            mexErrMsgTxt("Internal error. short is not 16-bits, or int is not 32-bits.");
        }
        checkieee = 0; // Everything checks out OK
    }

// Check input arguments for number and type
    
    if( nlhs > 1 ) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if( nrhs != 1 && nrhs != 2 ) {
        mexErrMsgTxt("Need one or two input arguments.");
    }
    if( mxIsSparse(prhs[0]) ) {
        mexErrMsgTxt("Sparse matrices not supported.");
    }
    if( nrhs == 2 ) {
        if( !mxIsChar(prhs[1]) )
            mexErrMsgTxt("2nd input argument must be char string naming desired class, or 'disp'.");
        if( !mxIsInt16(prhs[0]) && !mxIsUint16(prhs[0]) )
            mexErrMsgTxt("1st input argument must be uint16 or int16 class.");
        classname = mxArrayToString(prhs[1]);
        if( strcmp(classname,"double") == 0 ) {  // Check 2nd input for proper class name string
            classid = mxDOUBLE_CLASS;
        } else if( strcmp(classname,"single") == 0 ) {
            classid = mxSINGLE_CLASS;
        } else if( strcmp(classname,"int8") == 0 ) {
            classid = mxINT8_CLASS;
        } else if( strcmp(classname,"uint8") == 0 ) {
            classid = mxUINT8_CLASS;
        } else if( strcmp(classname,"int16") == 0 ) {
            classid = mxINT16_CLASS;
        } else if( strcmp(classname,"uint16") == 0 ) {
            classid = mxUINT16_CLASS;
        } else if( strcmp(classname,"int32") == 0 ) {
            classid = mxINT32_CLASS;
        } else if( strcmp(classname,"uint32") == 0 ) {
            classid = mxUINT32_CLASS;
        } else if( strcmp(classname,"int64") == 0 ) {
            classid = mxINT64_CLASS;
        } else if( strcmp(classname,"uint64") == 0 ) {
            classid = mxUINT64_CLASS;
        } else if( strcmp(classname,"char") == 0 ) {
            classid = mxCHAR_CLASS;
        } else if( strcmp(classname,"logical") == 0 ) {
            classid = mxLOGICAL_CLASS;
        } else if( strcmp(classname,"disp") == 0 ) {
            disp = 1;
        } else if( strcmp(classname,"isinf") == 0 ) {
            isinf = 1;
        } else if( strcmp(classname,"isnan") == 0 ) {
            isnan = 1;
        } else if( strcmp(classname,"isnormal") == 0 ) {
            isnormal = 1;
        } else if( strcmp(classname,"isdenormal") == 0 ) {
            isdenormal = 1;
        } else if( strcmp(classname,"eps") == 0 ) {
            eps = 1;
        } else {
            mexErrMsgTxt("2nd input argument must be char string naming desired class, or directive.");
        }
    }
    
    if( nrhs == 1 ) {  // Convert from input MATLAB variable to halfprecision ---------------------

        if( mxIsChar(prhs[0]) ) { // Check for directive
            classname = mxArrayToString(prhs[0]);
            if(        strcmp(classname,"openmp") == 0 ) {
#ifdef _OPENMP
                plhs[0] = mxCreateLogicalScalar(1);
#else
                plhs[0] = mxCreateLogicalScalar(0);
#endif
                return;
                
            } else if( strcmp(classname,"realmax") == 0 ) {
                plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL); // halfprecision stored as uint16
                *(UINT16_TYPE *)mxGetData(plhs[0]) = 0x7BFFu;
                return;
                
            } else if( strcmp(classname,"realmin") == 0 ) {
                plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL); // halfprecision stored as uint16
                *(UINT16_TYPE *)mxGetData(plhs[0]) = 0x0400u;
                return;
                
            } else if( strcmp(classname,"realmindenormal") == 0 ) {
                plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL); // halfprecision stored as uint16
                *(UINT16_TYPE *)mxGetData(plhs[0]) = 0x0001u;
                return;
                
            } else if( strcmp(classname,"version") == 0 ) {
#if TARGET_API_VERSION == R2017b
                plhs[0] = mxCreateString("2017b");
#else
                plhs[0] = mxCreateString("2018a");
#endif
                return;
            }
        }
        classid = mxGetClassID(prhs[0]);  // Check for supported class
        switch( classid ) {
        case mxDOUBLE_CLASS:
        case mxSINGLE_CLASS:
        case mxUINT64_CLASS:
        case mxINT64_CLASS:
        case mxUINT32_CLASS:
        case mxINT32_CLASS:
        case mxUINT16_CLASS:
        case mxINT16_CLASS:
        case mxUINT8_CLASS:
        case mxINT8_CLASS:
        case mxLOGICAL_CLASS:
        case mxCHAR_CLASS:
            break;
        default:
            mexErrMsgTxt("Unable to convert input argument to halfprecision.");
        }
        complexity = mxIsComplex(prhs[0]) ? mxCOMPLEX : mxREAL; // Get stats of input variable
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT16_CLASS, complexity); // halfprecision stored as uint16
        switch( classid ) {
        case mxDOUBLE_CLASS:  // Custom code for double class
#if TARGET_API_VERSION == R2017b
            doubles2halfp(mxGetData(plhs[0]),mxGetData(prhs[0]),numel,FE_TONEAREST,1);
            doubles2halfp(mxGetImagData(plhs[0]),mxGetImagData(prhs[0]),numel,FE_TONEAREST,1);
            if( complexity == mxCOMPLEX && matlab_version() >= 0x2018a && warn ) {
                mexWarnMsgTxt("Using complex variables in R2018a+ will cause extra deep copies");
                mexWarnMsgTxt("To avoid this recompile with -R2018a option");
                warn = 0;
            }
#else
            if( complexity == mxREAL ) {
                doubles2halfp(mxGetData(plhs[0]),mxGetData(prhs[0]),numel,FE_TONEAREST,1);
            } else {
                doubles2halfp(mxGetData(plhs[0]),mxGetData(prhs[0]),2*numel,FE_TONEAREST,1);
            }
#endif
            break;
        case mxSINGLE_CLASS:  // Custom code for single class
#if TARGET_API_VERSION == R2017b
            singles2halfp(mxGetData(plhs[0]),mxGetData(prhs[0]),numel,FE_TONEAREST,1);
            singles2halfp(mxGetImagData(plhs[0]),mxGetImagData(prhs[0]),numel,FE_TONEAREST,1);
            if( complexity == mxCOMPLEX && matlab_version() >= 0x2018a && warn ) {
                mexWarnMsgTxt("Using complex variables in R2018a+ will cause extra deep copies");
                mexWarnMsgTxt("To avoid this recompile with -R2018a option");
                warn = 0;
            }
#else
            if( complexity == mxREAL ) {
                singles2halfp(mxGetData(plhs[0]),mxGetData(prhs[0]),numel,FE_TONEAREST,1);
            } else {
                singles2halfp(mxGetData(plhs[0]),mxGetData(prhs[0]),2*numel,FE_TONEAREST,1);
            }
#endif
            break;
        case mxUINT64_CLASS:
        case mxINT64_CLASS:
        case mxUINT32_CLASS:
        case mxINT32_CLASS:
        case mxUINT16_CLASS:
        case mxINT16_CLASS:
        case mxUINT8_CLASS:
        case mxINT8_CLASS:
        case mxLOGICAL_CLASS:
        case mxCHAR_CLASS:  // All other classes get converted to single first
            mexCallMATLAB(1,lhs,1,(mxArray **)prhs,"single");
#if TARGET_API_VERSION == R2017b
            singles2halfp(mxGetData(plhs[0]),mxGetData(lhs[0]),numel,FE_TONEAREST,1);
            singles2halfp(mxGetImagData(plhs[0]),mxGetImagData(lhs[0]),numel,FE_TONEAREST,1);
#else
            if( complexity == mxREAL ) {
                singles2halfp(mxGetData(plhs[0]),mxGetData(lhs[0]),2*numel,FE_TONEAREST,1);
            } else {
                singles2halfp(mxGetData(plhs[0]),mxGetData(lhs[0]),numel,FE_TONEAREST,1);
            }
#endif
            mxDestroyArray(lhs[0]);
            break;
        }
    
    } else if ( disp ) {  // Display halfprecision input values to screen -------------------------
        rhs[0] = mxCreateDoubleScalar(0.0);
        rhs[1] = mxCreateString("Format");
        mexCallMATLAB(1,currentformat,2,rhs,"get");  // Get the current display format
        mxDestroyArray(rhs[1]);
        mxDestroyArray(rhs[0]);
        rhs[0] = mxCreateString("short");
        mexCallMATLAB(0,lhs,1,rhs,"format");  // Set the display format to short
        mxDestroyArray(rhs[0]);
        complexity = mxIsComplex(prhs[0]) ? mxCOMPLEX : mxREAL; // Get stats of input variable
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        rhs[0] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, complexity);
#if TARGET_API_VERSION == R2017b
        halfp2singles(mxGetData(rhs[0]), mxGetData(prhs[0]), numel);  // Convert input to single
        halfp2singles(mxGetImagData(rhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
        if( mxIsComplex(prhs[0]) ) {
            halfp2singles(mxGetData(rhs[0]), mxGetData(prhs[0]), 2*numel);  // Convert input to single
        } else {
            halfp2singles(mxGetData(rhs[0]), mxGetData(prhs[0]), numel);  // Convert input to single
        }
#endif
        mexCallMATLAB(0,lhs,1,rhs,"disp");  // Display it as single
        mxDestroyArray(rhs[0]);
        mexCallMATLAB(0,lhs,1,currentformat,"format");  // Restore current display format
        mxDestroyArray(currentformat[0]);
        
    } else if ( isinf ) {  // Return halfprecision isinf result -------------------------
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        plhs[0] = mxCreateLogicalArray(ndim, dims);
#if TARGET_API_VERSION == R2017b
        halfp2isinf(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isinf
        halfp2isinf(mxGetData(plhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
        if( mxIsComplex(prhs[0]) ) {
            halfp2isinf_interleaved_complex(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isinf
        } else {
            halfp2isinf(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isinf
        }
#endif
        
    } else if ( isnan ) {  // Return halfprecision isnan result -------------------------
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        plhs[0] = mxCreateLogicalArray(ndim, dims);
#if TARGET_API_VERSION == R2017b
        halfp2isnan(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isnan
        halfp2isnan(mxGetData(plhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
        if( mxIsComplex(prhs[0]) ) {
            halfp2isnan_interleaved_complex(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isnan
        } else {
            halfp2isnan(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isinf
        }
#endif
        
    } else if ( isnormal ) {  // Return halfprecision isdenormal result -------------------------
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        plhs[0] = mxCreateLogicalArray(ndim, dims);
#if TARGET_API_VERSION == R2017b
        halfp2isnormal(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isnormal
        halfp2isnormal(mxGetData(plhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
        if( mxIsComplex(prhs[0]) ) {
            halfp2isnormal_interleaved_complex(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isnormal
        } else {
            halfp2isnormal(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isinf
        }
#endif
        
    } else if ( isdenormal ) {  // Return halfprecision isdenormal result -------------------------
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        plhs[0] = mxCreateLogicalArray(ndim, dims);
#if TARGET_API_VERSION == R2017b
        halfp2isdenormal(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isdenormal
        halfp2isdenormal(mxGetData(plhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
        if( mxIsComplex(prhs[0]) ) {
            halfp2isdenormal_interleaved_complex(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isdenormal
        } else {
            halfp2isdenormal(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp isinf
        }
#endif
        
    } else if ( eps ) {  // Return halfp eps values -------------------------
        complexity = mxIsComplex(prhs[0]) ? mxCOMPLEX : mxREAL; // Get stats of input variable
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        plhs[0] = mxCreateNumericArray(ndim, dims, mxUINT16_CLASS, complexity); // halfprecision stored as uint16
#if TARGET_API_VERSION == R2017b
        halfp2eps(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp eps
        halfp2eps(mxGetData(plhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
        if( mxIsComplex(prhs[0]) ) {
            halfp2eps(mxGetData(plhs[0]), mxGetData(prhs[0]), 2*numel);  // halfp eps
        } else {
            halfp2eps(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // halfp eps
        }
#endif
        
    } else {  // Convert halfprecision to desired class ----------------------------------
        complexity = mxIsComplex(prhs[0]) ? mxCOMPLEX : mxREAL; // Get stats of input variable
        if( complexity == mxCOMPLEX && (classid == mxLOGICAL_CLASS || classid == mxCHAR_CLASS) ) {
            mexErrMsgTxt("Cannot convert complex into logical or char");
        }
        numel = mxGetNumberOfElements(prhs[0]);
        ndim = mxGetNumberOfDimensions(prhs[0]);
        dims = mxGetDimensions(prhs[0]);
        switch( classid ) {
        case mxDOUBLE_CLASS: // Custom code for double class
            plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, complexity);
#if TARGET_API_VERSION == R2017b
            halfp2doubles(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // Convert input to double
            halfp2doubles(mxGetImagData(plhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
            if( mxIsComplex(prhs[0]) ) {
                halfp2doubles(mxGetData(plhs[0]), mxGetData(prhs[0]), 2*numel);  // Convert input to double
            } else {
                halfp2doubles(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // Convert input to double
            }
#endif
            break;
        case mxSINGLE_CLASS: // Custom code for single class
            plhs[0] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, complexity);
#if TARGET_API_VERSION == R2017b
            halfp2singles(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // Convert input to single
            halfp2singles(mxGetImagData(plhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
            if( mxIsComplex(prhs[0]) ) {
                halfp2singles(mxGetData(plhs[0]), mxGetData(prhs[0]), 2*numel);  // Convert input to double
            } else {
                halfp2singles(mxGetData(plhs[0]), mxGetData(prhs[0]), numel);  // Convert input to double
            }
#endif
            break;
        case mxUINT64_CLASS:
        case mxINT64_CLASS:
        case mxUINT32_CLASS:
        case mxINT32_CLASS:
        case mxUINT16_CLASS:
        case mxINT16_CLASS:
        case mxUINT8_CLASS:
        case mxINT8_CLASS:
        case mxLOGICAL_CLASS:
        case mxCHAR_CLASS: // All other classes get converted to single first
            rhs[0] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, complexity);
#if TARGET_API_VERSION == R2017b
            halfp2singles(mxGetData(rhs[0]), mxGetData(prhs[0]), numel);  // Convert input to single
            halfp2singles(mxGetImagData(rhs[0]), mxGetImagData(prhs[0]), numel);  // ditto for imag
#else
            if( mxIsComplex(prhs[0]) ) {
                halfp2singles(mxGetData(rhs[0]), mxGetData(prhs[0]), 2*numel);  // Convert input to single
            } else {
                halfp2singles(mxGetData(rhs[0]), mxGetData(prhs[0]), numel);  // Convert input to single
            }
#endif
            mexCallMATLAB(1,plhs,1,rhs,classname);  // Now convert the single to desired class
            mxDestroyArray(rhs[0]);
            break;
        default:
            mexErrMsgTxt("Unable to convert input argument to halfprecision.");
        }
        mxFree(classname);
    }
    
}

//-----------------------------------------------------------------------------

void singles2halfp(void *target, void *source, ptrdiff_t numel, int rounding_mode, int is_quiet)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) target; // Type pun output as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) source; // Type pun input as an unsigned 32-bit int
    UINT16_TYPE    hs, he, hm, hr, hq;
    UINT32_TYPE x, xs, xe, xm, xt, zm, zt, z1;
    int hes, N;
    ptrdiff_t i;
    
    if( source == NULL || target == NULL ) { // Nothing to convert (e.g., imag part of pure real)
        return;
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
}

//-----------------------------------------------------------------------------

void doubles2halfp(void *target, void *source, ptrdiff_t numel, int rounding_mode, int is_quiet)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) target; // Type pun output as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) source; // Type pun input as an unsigned 32-bit int
    UINT16_TYPE    hs, he, hm, hq, hr;
    UINT32_TYPE x, xs, xe, xm, xn, xt, zm, zt, z1;
    int hes, N;
    ptrdiff_t i;

    if( source == NULL || target == NULL ) { // Nothing to convert (e.g., imag part of pure real)
        return;
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
}

//-----------------------------------------------------------------------------

void halfp2singles(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) target; // Type pun output as an unsigned 32-bit int
    UINT16_TYPE h, hs, he, hm;
    UINT32_TYPE xs, xe, xm;
    INT32_TYPE xes;
    ptrdiff_t i;
    int e;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
    
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
}

//-----------------------------------------------------------------------------

void halfp2doubles(void *target, void *source, ptrdiff_t numel)
{
    UINT16_TYPE *hp = (UINT16_TYPE *) source; // Type pun input as an unsigned 16-bit int
    UINT32_TYPE *xp = (UINT32_TYPE *) target; // Type pun output as an unsigned 32-bit int
    UINT16_TYPE h, hs, he, hm;
    UINT32_TYPE xs, xe, xm;
    INT32_TYPE xes;
    ptrdiff_t i, j;
    int e;
    
    if( source == NULL || target == NULL ) // Nothing to convert (e.g., imag part of pure real)
        return;
    
    xp += next;  // Little Endian adjustment if necessary
    
#pragma omp parallel for private(xs,xe,xm,h,hs,he,hm,xes,e,j)
    for( i=0; i<numel; i++ ) {
        j = 2*i;
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

#if TARGET_API_VERSION == R2018a

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

#endif
