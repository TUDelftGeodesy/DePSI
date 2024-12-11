/* matlab_version.c */
/*************************************************************************************
 *
 * MATLAB (R) is a trademark of The Mathworks (R) Corporation
 *
 * Filename:    matlab_version.c
 * Programmer:  James Tursa
 * Version:     1.00
 * Date:        May 17, 2018
 * Copyright:   (c) 2018 by James Tursa, All Rights Reserved
 *
 * Change Log:
 * 2018/May/17 --> 1.00, Initial Release (stripped down to function only)
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
 * matlab_version.c contains a function for determining MATLAB runtime version:
 * (version numbering uses hex values to support the "a" and "b" suffix)
 *
 *   0x00000 = UNKNOWN
 *   0x2006a = R2006a or earlier
 *   0x2006b = R2006b
 *   0x2007a = R2007a
 *           :
 *          etc
 *
 *****************************************************************************/

#include <stdio.h>  /* sscanf */

int matlab_version(void)
{
    mxArray *mx_version = NULL;
    char cversion[100+1] = {'\0'};
    char *c = cversion;
    static int version = 0x0000;

    if( !version ) {
        if( !mexCallMATLAB(1,&mx_version,0,NULL,"version") ) {
            mxGetString(mx_version,cversion,100);
            while( *c ) {
                if( *c++ == 'R' ) {
                    if( 100 - (c-cversion) >= 5 ) { 
                        sscanf(c,"%5x",&version);
                    }
                    break;
                }
            }
            mxDestroyArray(mx_version);
        }
    }
    return version;
}
