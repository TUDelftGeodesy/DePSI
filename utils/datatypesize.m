function b = datatypesize(c)
% B=DATATYPESIZE(C)
% return number of bytes B based on string C.
% extended with cpx types for complex, see also freadbk.
% -1 is returned for unknown type
%
% help fread
%        MATLAB    C or Fortran     Description
%        'uchar'   'unsigned char'  unsigned character,  8 bits.
%        'schar'   'signed char'    signed character,  8 bits.
%        'int8'    'integer*1'      integer, 8 bits.
%        'int16'   'integer*2'      integer, 16 bits.
%        'int32'   'integer*4'      integer, 32 bits.
%        'int64'   'integer*8'      integer, 64 bits.
%        'uint8'   'integer*1'      unsigned integer, 8 bits.
%        'uint16'  'integer*2'      unsigned integer, 16 bits.
%        'uint32'  'integer*4'      unsigned integer, 32 bits.
%        'uint64'  'integer*8'      unsigned integer, 64 bits.
%        'single'  'real*4'         floating point, 32 bits.
%        'float32' 'real*4'         floating point, 32 bits.
%        'double'  'real*8'         floating point, 64 bits.
%        'float64' 'real*8'         floating point, 64 bits.
% 
%    The following platform dependent formats are also supported but
%    they are not guaranteed to be the same size on all platforms.
%        'char'    'char*1'         character,  8 bits (signed or unsigned).
%        'short'   'short'          integer,  16 bits.
%        'int'     'int'            integer,  32 bits.
%        'long'    'long'           integer,  32 or 64 bits.
%        'ushort'  'unsigned short' unsigned integer,  16 bits.
%        'uint'    'unsigned int'   unsigned integer,  32 bits.
%        'ulong'   'unsigned long'  unsigned integer,  32 bits or 64 bits.
%        'float'   'float'          floating point, 32 bits.
%        'float16' 'real*2'         floating point, 16 bits.
%
%
% See also FREAD, FREADBK, FSEEK, FTELL
%


% $Revision: 1.4 $  $Date: 2001/09/28 14:24:31 $
% Bert Kampes, 4/3/00

%%% see fread
if (nargin~=1 | ~ischar(c)) helphelp; error('An error occured'); end;


%%% Check complex type request.
cpx=1; % 2 for complex type
if (length(c)>3)
  if (c(1:3)=='cpx')
    c   = c(4:length(c));
    cpx = 2;
  end
end

switch (lower(c))
  case {'uchar','unsigned char','int8','integer*1'}, 	b=1;
  case {'int16','integer*2'}, 				b=2;
  case {'int32','integer*4'}, 				b=4;
  case {'int64','integer*8'}, 				b=8;
  case {'uint8','integer*1'}, 				b=1;
  case {'uint16','integer*2'}, 				b=2;
  case {'uint32','integer*4'}, 				b=4;
  case {'uint64','integer*8'}, 				b=8;
  case {'float16','real*2'},                               b=2;
  case {'single','real*4'},				b=4;
  case {'float32','real*4'},				b=4;
  case {'double','real*8'},				b=8;
  case {'float64','real*8'},				b=8;
  case {'char','char*1'},				b=1;
  case {'short'},					b=2;
  case {'int'},						b=4;
  case {'long'},					b=4;
  case {'ushort'},					b=2;
  case {'uint'},					b=4;
  case {'ulong'},					b=4;
  case {'float'},					b=4;
  otherwise warning('Unknown size for type ...'), 	b=-1;
end

%%% Correct for complex type
if ( b>0 )
  b=b*cpx;
end;

%%% EOF






