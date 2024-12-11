/* matlab_version.h */
/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Filename:    matlab_version.h
 * Programmer:  James Tursa
 * Version:     5.00
 * Date:        March 13, 2020
 * Copyright:   (c) 2018, 2019, 2020 by James Tursa, All Rights Reserved
 *
 * Change Log:
 * 2012/Mar/21 --> 1.00, Initial Release
 * 2018/Apr/20 --> 2.00, Updated the comments for more information
 *                       Updated use of MATLAB_VERSION EMLRT_VERSION_INFO
 *                       Updated conditional use of MX_API_VER
 * 2018/May/16 --> 3.00, Updated comments, mostly additional R2018a info
 *                       Added prototype for matlab_version() function
 * 2018/Aug/15 --> 4.00, Updated comments for mxArrayToString and
 *                       mxArrayToUTF8String
 * 2019/May/20 --> 4.10, Updated comments for R2019a CrossLink field
 * 2020/Mar/13 --> 5.00, Updated comments for R2019b ndim field
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
 * matlab_version.h is an include file that determines, at compile time, what
 * version of MATLAB is being used for the compilation. It has been tested to
 * work under R2006a - R2011a, 32-bit Windows, and some R2014a - R2017b 64-bit
 * Windows versions. It may work under other versions but I don't have access
 * to test them. If you do run into any problems on any versions, please e-mail
 * me and let me know so that I can update this file. The ultimate goal is to
 * have one file that works for all platforms and MATLAB versions (within reason).
 *
 * There are two macros defined with this pre-processor code discussed below:
 *
 *   MATLAB_VERSION
 *   TARGET_API_VERSION
 *
 * And the following function prototype:
 *
 *   int matlab_version(void);  <-- The function code is in matlab_version.c
 *
 * Usage: Simply include this file at the top of your code (usually after you
 * have already included mex.h or engine.h) and you will get a macro defined
 * called MATLAB_VERSION that will contain an integer representation of the
 * MATLAB version that is being used for the compile:
 * (version numbering uses hex values to support the "a" and "b" suffix)
 *
 *   0x00000 = UNKNOWN
 *   0x2006a = R2006a or earlier
 *   0x2006b = R2006b
 *   0x2007a = R2007a
 *           :
 *          etc
 *
 * I don't have access to versions earlier than R2006a, so I don't yet have
 * logic to differentiate those versions. If you can come up with something
 * please let me know. Thanks.
 *
 * Relies on macros that are defined in matrix.h, tmwtypes.h, blas.h, emlrt.h,
 * lapack.h (you don't have to explicitly include them since this header will
 * include them if necessary). Also, this macro can essentially be disabled
 * if the macro MATLAB_VERSION has already been defined (e.g., by using the
 * -DMATLAB_VERSION=blahblah option on the mex command line).
 *
 * The purpose of knowing the MATLAB version at compile time is so that your
 * code can make intelligent decisions with regards to which API functions
 * are available (e.g., mxGetProperty and mxSetProperty) for linking. That
 * way you can either write substitude functions, dummy stub functions, or
 * whatever makes sense for your application, and use #if tests on the value
 * of MATLAB_VERSION to determine which code to include at compile time.
 *
 * Important version specific changes to be aware of:
 *
 * R2006b - First introductin of specific types for integer sizes, etc:
 *              mwSize, mwIndex, mwSignedIndex
 *          And the related macros:
 *              MWSIZE_MAX, MWSIZE_MIN, MWINDEX_MAX, MWINDEX_MIN,
 *              MWSINDEX_MAX, MWSINDEX_MIN
 *
 * R2007b - The combination blas/lapack library is split into two libraries.
 *          First introduction of the blas.h and lapack.h header files. However,
 *          all of the functions that return complex scalar values have either
 *          kludged double (for complex float return) or incorrect double
 *          (for complex double return). The complex float functions *might*
 *          be coaxed into extracting the two individual real & imag floats
 *          from the double return value as long as the bit pattern doesn't
 *          inadvertantly match a signaling NaN. Very suspect. The double
 *          complex return values I don't see how could ever work properly.
 *              cdotc, cdotu, zdotc, zdotu
 *
 * R2008a - First introduction of the new classdef OOP. So opaque objects
 *          can now be used, and the API functions mxGetProperty and
 *          mxSetProperty first appear (although these function only work
 *          with *deep* copies, not pointers to original variables). The
 *          API functions mxGetPropertyShared and mxSetPropertyShared
 *          (both undocumented) are available for C++ only (no C export).
 *
 * R2008b - First version where the Ir and Jc arrays are shared in shared
 *          data copies of sparse matrices. Prior to this version only the
 *          Pr and Pi arrays were shared (each sparse shared data copy
 *          variable used to have deep copies of the Ir and Jc arrays).
 *
 *          First version where mxCreateUninitNumericArray is available.
 *
 * R2009a - Change in the mxArray header definition. The first item used
 *          to be a pointer to a C-style null terminated string containing
 *          the name of the variable. Now this item is NULL until R2011a.
 *          Thus all of the older deprecated API functions that depended
 *          on this become obsolete and are no longer included in the API:
 *            matPutArray, matPutNextArray, matPutArrayAsGlobal,
 *            matGetArray, matGetArrayHeader, matGetNextArray,
 *            matGetNextArrayHeader, matDeleteArray
 *
 * R2009b - Fixes the blas incorrect complex scalar return values in the
 *          blas.h file so that they correctly reflect the paired float
 *          or double return values. However, apparently on some platforms
 *          these functions do *not* return any value. Instead, the first
 *          argument is a pointer to the value pair to be returned.
 *
 * R2010b - The struct mxArray_tag definition is removed from matrix.h
 *
 * R2011a - Change in the mxArray header definition. The first item used
 *          to be NULL (since R2009a) but is now a reverse CrossLink pointer.
 *          I.e., the CrossLink that is used to connect shared data copies
 *          of variables is now a double linked list (previously was only a
 *          single forward linked list).
 *
 * R2011b - Change in cell array Pi data pointer. This used to always be
 *          NULL since nothing else was needed at the data level to keep
 *          track of the mxArray data addresses. For R2011b, the Pi data
 *          pointer now contains a value. Unknown what this is used for.
 *          (maybe to make it look internally more like a struct array?)
 *
 *          Change in struct array Pi data pointer. The location of the
 *          fieldname information behind this pointer has changed. E.g.,
 *          The old 32-bit format was:
 *            Pi[0] = Number of fields
 *            Pi[1] = Number of variables sharing the Pi memory block
 *            Pi[2] = 0
 *            Pi[3] = 1st fieldname pointer (char *) C-style null terminated
 *            Pi[4] = 2nd fieldname pointer (char *) etc.
 *          The new 32-bit format is:
 *            Pi[0] = ?
 *            Pi[1] = Pointer to the fieldname information
 *            Pi[1][0] = Number of fields
 *            Pi[1][1] = Number of variables sharing the Pi memory block
 *            Pi[1][2] = 0
 *            Pi[1][3] = 1st fieldname pointer (char *) C-style null terminated
 *            Pi[1][4] = 2nd fieldname pointer (char *) etc.
 *          The new 64-bit format is:
 *            Pi[0] = ?
 *            Pi[1] = Pointer to the fieldname information
 *            Pi[1][0] = Number of fields
 *            Pi[1][1] = Number of variables sharing the Pi memory block
 *            Pi[1][2] = 1st fieldname pointer (char *) C-style null terminated
 *            Pi[1][3] = 2nd fieldname pointer (char *) etc.
 *          Shared data copies of struct variables no longer share Pi.
 *          Note that fieldname pointers are shared among struct variables.
 *          It is as if there is a general pool of fieldnames strings that
 *          MATLAB keeps in the background. If you create a new struct with
 *          any fieldnames that match strings in this existing pool, then
 *          the fieldname string pointer for that particular fieldname will
 *          be for that existing string. This fieldname string pool is
 *          maintained in the background separate of variables. I.e., even
 *          if you clear all the variables, the pool remains intact and if
 *          you subsequently create new struct variables with fieldnames
 *          that match anything in the pool, then again pointers to those
 *          pre-existing strings will be used for this new struct variable.
 *          The Pi[1] (or Pi[1][1]) value seems to behave kind of like a
 *          reference counter for the Pi memory block. I.e., it keeps track
 *          of how many variables are sharing the Pi memory block.
 *
 * R2012a - Change in struct array Pi data pointer. The number of variables
 *          sharing the Pi memory block is no longer used. It is 0. E.g.,
 *          The new 32-bit format is:
 *            Pi[0] = ?
 *            Pi[1] = Pointer to the fieldname information
 *            Pi[1][0] = Number of fields
 *            Pi[1][1] = 0
 *            Pi[1][2] = 0
 *            Pi[1][3] = 1st fieldname pointer (char *) C-style null terminated
 *            Pi[1][4] = 2nd fieldname pointer (char *) etc.
 *          The new 64-bit format is:
 *            Pi[0] = ?
 *            Pi[1] = Pointer to the fieldname information
 *            Pi[1][0] = Number of fields
 *            Pi[1][1] = 0
 *            Pi[1][2] = 1st fieldname pointer (char *) C-style null terminated
 *            Pi[1][3] = 2nd fieldname pointer (char *) etc.
 *
 * R2012b - Change in struct array Pi data pointer. Something else is now
 *          stored in the Pi[1][1] area. Looks like a pointer. Unknown what
 *          it points to.
 *          The new 32-bit format is:
 *            Pi[0] = ?
 *            Pi[1] = Pointer to the fieldname information
 *            Pi[1][0] = Number of fields
 *            Pi[1][1] = UNKNOWN POINTER
 *            Pi[1][2] = 0
 *            Pi[1][3] = 1st fieldname pointer (char *) C-style null terminated
 *            Pi[1][4] = 2nd fieldname pointer (char *) etc.
 *          The new 64-bit format is:
 *            Pi[0] = ?
 *            Pi[1] = Pointer to the fieldname information
 *            Pi[1][0] = Number of fields
 *            Pi[1][1] = UNKNOWN POINTER
 *            Pi[1][2] = 1st fieldname pointer (char *) C-style null terminated
 *            Pi[1][3] = 2nd fieldname pointer (char *) etc.
 *
 * R2014a - First version where the undocumented C mxCreateReference function
 *          is removed from the API library. 
 *
 * R2015a - First version where the previously undocumented API functions
 *          mxCreateUninitNumericMatrix and mxCreateUninitNumericArray
 *          appear in the official documentation (they have been part of
 *          the API library for several years, but now become official).
 *
 * R2015b - First version where reference copies (sharing the actual mxArray
 *          header struct) are used among top level workspace variables and
 *          the method of passing variables to mex routines changes. There
 *          appears to have been a major change in the way that MATLAB works
 *          with shared data copy and reference copy variables with this
 *          release, and also a major change in the way that MATLAB passes
 *          input arguments to mex functions and through the mexCallMATLAB
 *          API function. For instance, normal workspace variables created
 *          with literal assignments (e.g., x = 1:5) can have a reference
 *          copy of the 1:5 saved in the background, which MATLAB will use
 *          in the future if the same pattern (x = 1:5) is detected. That
 *          is, the 1:5 is not created on the fly anymore ... it is only
 *          created once, saved in the background, and then re-used in later
 *          assignments when the assignment pattern is the same. This is
 *          just one example of the R2015b changes. I don't know what all
 *          of the changes are. E.g.,
 *
 *          R2015a- Top level workspace variables only shared data through
 *                  the shared data copy mechanism. Reference copies, where
 *                  the mxArray struct header is also shared, were only used
 *                  for cell and struct field elements. E.g.,
 *                    y = x; % y is a shared data copy of x
 *          R2015b+ Top level workspace variables can be reference copies of
 *                  other variables (e.g., for direct assignments). E.g.,
 *                    y = x; % y is a reference copy of x
 *
 *          R2015a- Top level workspace variables were always passed to mex
 *                  routines by original mxArray address. I.e., you always
 *                  got the address of the original workspace variable being
 *                  passed in. Cell and struct field elements could be passed
 *                  in as shared data copies. E.g.,
 *                    mymex(x); % prhs[0] is the address of x in the workspace
 *                    mymex(c{1}); % prhs[0] could be the address of shared data copy of c{1}
 *          R2015b+ Top level workspace variables are now passed to mex
 *                  routines as shared data copies. I.e., you don't get the
 *                  address of the original workspace mxArray header anymore.
 *                  Cell and struct field element passing is the same, and
 *                  could be passed in as shared data copies. E.g.,
 *                    mymex(x); % prhs[0] is the address of shared data copy of x
 *                    mymex(c{1}); % prhs[0] could be the address of shared data copy of c{1}
 *                  Also, the act of calling a mex function seems to create
 *                  extra reference copies of the input arguments for some
 *                  reason. I.e., the RefCount of the input arguments gets
 *                  bumped up by 1, at least for the duration of the call.
 *                  I don't know the purpose of these changes (maybe to protect
 *                  against the user inadvertenly calling mxDestroyArray on
 *                  the prhs[ ] variables?)
 *
 *          Last version that has 32-bit MATLAB. All later versions are 64-bit.
 *
 * R2016b - First version where the string class is introduced.
 *          Last version where the allocated memory from the mxArrayToString
 *          and mxArrayToUTF8String functions are persistent and are *not*
 *          on the garbage collection list. You must manually mxFree this
 *          memory or you will get a permanent memory leak that can only be
 *          recovered by restarting MATLAB.
 *
 * R2017a - First version where the macro TARGET_API_VERSION is defined.
 *          First version where the allocated memory from the mxArrayToString
 *          and mxArrayToUTF8String functions are non-persistent and *are*
 *          on the garbage collection list. This memory will be automatically
 *          freed by MATLAB when the mex function exits.
 *
 * R2017b - The macro TARGET_API_VERSION is no longer defined.
 *
 * R2018a - First version where interleaved complex memory storage is used.
 *          The Pi data pointer has been removed from the mxArray header.
 *          You can either link with a library that supports the old separate
 *          complex data model interface (but will do a copy-in/copy-out of
 *          the complex data), or you can link with a library that supports
 *          the new interleaved complex memory model (to avoid the copy-in/
 *          copy-out of complex data). Any hacks (such as modifying prhs[]
 *          input variables inplace) using the separate complex library
 *          are likely to have major problems (will not work or will crash
 *          MATLAB). The default mex compiling option is -R2017b, which
 *          means you are linking against the API that does copy-in/copy-out
 *          for complex variables. You have to explicitly use the -R2018a
 *          option in the mex command to link against the interleaved
 *          complex API. API function changes for the -R2018a option:
 *
 *          C Matrix API functions not supported in -R2018a API:
 *              mxGetPi
 *              mxGetImagData
 *              mxSetPi
 *              mxSetImagData
 *          The behavior of these C Matrix API functions changes in -R2018a API:
 *              mxGetData
 *              mxSetData
 *              mxGetElementSize
 *              mexCallMATLAB
 *              mexGetVariable
 *              mexGetVariablePtr
 *              mexPutVariable
 *          C Matrix API functions to be phased out:
 *              mxGetPr
 *              mxSetPr
 *
 *          First version where a whole slew of type specific data pointer
 *          routines are available. E.g., mxGetDoubles, mxGetComplexDoubles,
 *          mxGetInt8s, mxGetComplexInt8s, etc. There are real and imag
 *          Get and Set routines for every numeric data type introduced.
 *
 *          First version where the udocumented mxGetPropertyShared and
 *          mxSetPropertyShared functions appear in the C API library (they
 *          have been part of the C++ API library for years). To link with
 *          them you will need to undefine their override macros. E.g.,
 *
 *              #include "mex.h"
 *              #undef mxGetPropertyShared
 *              #undef mxSetPropertyShared
 *              mxArray *mxGetPropertyShared(const mxArray *pa, mwIndex index, const char *propname);
 *              void mxSetPropertyShared(mxArray *pa, mwIndex index, const char *propname, const mxArray *value);
 *
 *          First version where the Pi data pointer is removed from mxArray.
 *          Note that the mxArray of all variables is the new R2018a version
 *          with the Pi data pointer removed, even with the -R2017b option.
 *          That is, with the -R2017b option the separate Pr and Pi data
 *          areas are actually kept off to the side somewhere as deep data
 *          copies and the API library functions deal with that in the
 *          background for you. For complex variables this means a deep data
 *          copy-in/copy-out when entering and exiting the mex routine. But
 *          the underlying mxArray for every variable is still the R2018a
 *          version without the Pi pointer even when using -R2017b option.
 *          
 *          First version where there are major changes to the meaning of
 *          the mxArray header bit flags:
 *          R2017b and earlier:
 *              bit  0 = is scalar double full
 *              bit  2 = is empty double full
 *              bit  4 = is temporary
 *              bit  5 = is sparse
 *              bit  9 = is numeric
 *              bits 24 - 31 = User Bits
 *          R2018a and later:
 *	            bit  0 = is scalar double full
 *              bit  4 = is sparse
 *              bit  7 = is numeric
 *              bit 11 = is complex
 *              bits 12 - 17 = Various values related to byte usage of non-complex variables
 *              Number of Data Bytes      Bits
 *              flags bits 17-12 values for the number of data bytes
 *                      0                000000
 *                      1 -    32        000001
 *                     33 -    64        000011
 *                     65 -    96        000101
 *                     97 -   128        000111
 *                    129 -   160        001001
 *                    161 -   192        001011
 *                    193 -   224        001101
 *                    225 -   256        001111
 *                    257 -   352        010001
 *                    353 -   960        010011
 *                    961 -  2048        010101
 *                   2049 -  4096        010111
 *                   4097 -  6144        011001
 *                   6145 -  8192        011011
 *                   8193 -  9216        011101
 *                   9217 - 10240        011111
 *                  10241 - 13312        100001
 *                  13313 - 16384        100011
 *                  16385 -              000000
 *              bits 24 - 31 = User Bits
 *          So it used to be that you only had to check if Pi is not NULL to
 *          see if a numeric variable is complex. Now you have to check bit
 *          flag 11 for 1 or 0 to see if a numeric variable is complex.
 *
 *          First version where the following undocumented API functions
 *          are prevented by the mex command from linking with mex code:
 *          (the functions still appear to be in the API library, they are
 *           just deliberately prevented from linking with mex code).
 *
 *              mxCreateSharedDataCopy
 *              mxCreateUninitDoubleMatrix
 *              mxFastZeros
 *              mxUnreference
 *              mxUnshareArray
 *              mxGetPropertyShared
 *              mxSetPropertyShared
 *
 *          The above functions are apparently slated to be removed from
 *          the API library in future versions of MATLAB. The mechanism
 *          for preventing the linking is to define macros with the same
 *          names as these functions. To re-enable the ability to link,
 *          simply undefine the macros. E.g., in your C source code:
 *
 *              #include "mex.h"
 *              #undef mxCreateSharedDataCopy
 *              mxArray *mxCreateSharedDataCopy(const mxArray *);
 *
 *          Even though the above functions are still in the API library and
 *          you can link with them using the above hack, it is unknown if
 *          they actually work correctly in R2018a with either the -R2017b
 *          or -R2018a API library compile option, so caveat emptor. And
 *          keep in mind that even this hack will break in future versions
 *          if the functions are physically removed from the library.
 *
 *          First version where the macros R2017b, R2018a are available, and
 *          the macro TARGET_API_VERSION appears again. They are defined as:
 *
 *            #define R2017b 700
 *            #define R2018a 800
 *
 *          The macro TARGET_API_VERSION is then apparently set to either
 *          R2017b or R2018a by the mex command as follows
 *
 *            TARGET_API_VERSION = 700 if mex option -R2017b is present, or
 *                                     if neither -R2017b or -R2018a is used
 *            TARGET_API_VERSION = 800 if mex option -R2018a is present
 *
 *          So the mex programmer could use the following in their code to
 *          make it backward compatible and make TARGET_API_VERSION available
 *          for all versions:
 *
 *            #if !defined(R2017b)
 *            #  define R2017b 700
 *            #endif
 *            #if !defined(R2018a)
 *            #  define R2018a 800
 *            #endif
 *            #if !defined(TARGET_API_VERSION)
 *            #  define TARGET_API_VERSION R2017b
 *            #endif
 *
 *          This macro construct is included as part of this current file.
 *
 * R2019a - First version where the CrossLink field in the mxArray header
 *          has been moved next to the reverse CrossLink field. I.e.,
 *          R2018b-:
 *                  void *RevCrossLink;
 *                  mxClassID ClassID;
 *                  int VariableType;
 *                  mxArray *CrossLink;
 *                     :
 *                    etc.
 *          R2019a+:
 *                  void *RevCrossLink;
 *                  mxArray *CrossLink;
 *                  mxClassID ClassID;
 *                  int VariableType;
 *                     :
 *                    etc.
 *
 * R2019b - First version where the ndim field in the mxArray header
 *          has been moved on top of the N field in a union and one
 *          of the flags bits is now used to indicate 2D matrix. I.e.,
 *          R2019a:
 *                  void *RevCrossLink;
 *                  mxArray *CrossLink;
 *                  mxClassID ClassID;
 *                  int VariableType;
 *                  size_t ndim;
 *                  unsigned int RefCount;
 *                  unsigned int flags;
 *                  union {
 *                    size_t M;
 *                    size_t *dims;
 *                  } Mdims;
 *                  size_t N;
 *                  void *pr;
 *                  union {
 *                    mwIndex *ir;
 *                    mxClassID ClassID;
 *                    char *ClassName;
 *                  } irClassNameID;
 *                  union {
 *                    mwIndex *jc;
 *                    mxClassID ClassID;
 *                  } jcClassID;
 *                  size_t nzmax;
 *          R2019b+:
 *                  void *RevCrossLink;
 *                  mxArray *CrossLink;
 *                  mxClassID ClassID;
 *                  int VariableType;
 *                  unsigned int RefCount;
 *                  unsigned int flags;
 *                  union {
 *                    size_t M;
 *                    size_t *dims;
 *                  } Mdims;
 *                  union {
 *                    size_t N;
 *                    size_t ndim;
 *                  } Nndim;
 *                  void *pr;
 *                  union {
 *                    mwIndex *ir;
 *                    mxClassID ClassID;
 *                    char *ClassName;
 *                  } irClassNameID;
 *                  union {
 *                    mwIndex *jc;
 *                    mxClassID ClassID;
 *                  } jcClassID;
 *                  size_t nzmax;
 *          ndim is only populated if the appropriate flags bit is
 *          reset to 0. If that flags bit is set to 1, then the union
 *          spot is populated with N.
 *
 *          First version where the half( ) data type is provided. However,
 *          it has been implemented as an opaque class so there is no way
 *          to get at the data pointers inside a mex routine. The mxGetData
 *          function will not work for this.
 *
 ****************************************************************************/

#if !defined(MATLAB_VERSION)

#if !defined(matrix_h)
#    include "matrix.h"
#endif

#if !defined(MWSIZE_MAX)
#    define MATLAB_VERSION 0x2006a /* R2006a or earlier */

#elif defined(MX_API_VER) && (MX_API_VER < 0x07040000)
#    define MATLAB_VERSION 0x2006b /* R2006b */

#elif !defined(FMT_PTRDIFF_T)
#    define MATLAB_VERSION 0x2007a /* R2007a */

#elif !defined(CUINT64_T)
#    define MATLAB_VERSION 0x2007b /* R2007b */

#elif defined(mxSetLogical)
#    define MATLAB_VERSION 0x2008a /* R2008a */

#else

#    if !defined(blas_h)
#        include "blas.h"
#    endif
#    if !defined(lapack_h)
#        include "lapack.h"
#    endif

#    if !defined(MATHWORKS_MATRIX_MATRIX_PUB_FWD_H)
#        if defined(CHAR16_T)
#            if !defined(COMPLEX_TYPES)
#                define MATLAB_VERSION 0x2008b /* R2008b */
#            elif !defined(cgeqr2p)
#                define MATLAB_VERSION 0x2010b /* R2010b */
#            else
#                define MATLAB_VERSION 0x2011a /* R2011a */
#            endif
#        else
#            include "emlrt.h"
#            define MATLAB_VERSION EMLRT_VERSION_INFO /* R2011b or later */
#        endif
#    else
#        if !defined(COMPLEX_TYPES)
#            define MATLAB_VERSION 0x2009a /* R2009a */
#        elif !defined(cgbequb)
#            define MATLAB_VERSION 0x2009b /* R2009b */
#        else
#            define MATLAB_VERSION 0x2010a /* R2010a */
#        endif
#    endif

#endif

#endif /* if !defined(MATLAB_VERSION) */

/* For R2018a backward compatibility */
#if !defined(R2017b)
#  define R2017b 700
#endif
#if !defined(R2018a)
#  define R2018a 800
#endif
#if !defined(TARGET_API_VERSION)
#  define TARGET_API_VERSION R2017b
#endif

/* Function prototype for runtime version in matlab_version.c */
int matlab_version(void);
