function [data, count] = freadbk(infile, lines, bkformat, r0, rN, c0, cN);
% FREADBK  --  Read binary data file.
%
%   DATA = FREADBK pops up a menu to select a filename, asks for the number
%   of lines and format of the input file, and reads in the data (double).
%
%   DATA = FREADBK(FILENAME) fopens the specified file with read permission and
%   returns a row vector DATA.  Format is assumed to be float32 format.
%
%   [DATA, COUNT] = FREADBK(FILENAME) optionally returns the number of elements
%   successfully read.
%
%   [DATA, COUNT] = FREADBK(FILENAME, NUMLINES) includes an optional NUMLINES
%   argument for the number of lines in the file.  Matrix DATA is returned of
%   size NUMLINES rows and the appropriate number of columns.
%
%   [DATA, COUNT] = FREADBK(FILENAME, NUMLINES, BKFORMAT) reads in a file in
%   the format BKFORMAT.  Which is one of the formats specified in FREAD (see
%   help fread), or one of these with a prepended 'cpx' for complex data types.
%   (DATA must be stored major row order, pixel interleaved):
%	'cpxfloat32'	complex floating point, 32 bits.
%	'cpxint16'	complex signed 16 bit short integers, ...
%	'cpx...'	...
%      ('mph'           format: use 'cpxfloat32')
%      ('hgt'           format: use 'freadhgt' function)
%
%   [DATA, COUNT] = FREADBK(FILENAME, NUMLINES, BKFORMAT, r0, rN, c0, cN) only
%   read in the part between: row r0 to row rN, and between column c0 to cN.
%   Indexes start at 1, stop at numlines/width. r0=0 all rows, c0=0 means all
%   columns.
%
%   Examples:
%     To read in a file with 2 channels (complex), stored pixel interleaved, ie,
%     (row1: RE1 - IM1 - RE2 - IM2 - RE3 - ...)
%     (row2: RE1 - IM1 - RE2 - IM2 - RE3 - ...)
%     (      ...                              )
%     in float32 format, which has a height of 100 lines, and a width of 250
%     complex pixels, use:
%       D = freadbk('filename',100,'cpxfloat32');
%
%     To crop this file while reading, between rows 1:10, columns 101:110, use:
%       D = freadbk('filename',100,'cpxfloat32', 1, 10, 101, 110);
%
%   See also FWRITEBK, FOPEN, FREAD, FWRITE, LOAD, SAVE, FPRINTF,
%   FSEEK, FTELL, FREADHGT, FWRITEHGT
%

% $Revision: 1.12 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 04-Mar-2000
%
% Improved by Freek van Leijen, 16-12-04, with a considerable speed
% up (especially for large crops)

%%% Handle input. (varargin?)
%if (nargin >  3) error('working on this... not yet ok, sorry.'); end;
false=0; true=1; 
complextype=false; readwholefile=false;%		defaults

if (nargin <  7) cN = 0; end;%				fall through (not...)
if (nargin <  6) c0 = 0; end;%				fall through
if (nargin <  5) rN = 0; end;%				fall through 
if (nargin <  4) r0 = 0; readwholefile=true; end;%	fall through

if (nargin <  3) bkformat='float32'; end;%		fall through
if (nargin <  2) lines=1; end;%				fall through
if (nargin <  1)
    [infile, inpath] = uigetfile('*', 'Select binary input float file', 0,0);
    infile   = [inpath,infile];
    lines    = input('Enter number of lines in file: ');
    bkformat = input('Format (enter between single quotes): ');% make a gui...
end;
%%% Check if whole file should be read regardless of r0,rNc0,cN.
if (c0==1 & r0==1 & rN==lines)
  bytesperelem = datatypesize(bkformat);%	in bytes
  filesize     = fsize(infile);%		in bytes
  filewidth    = filesize/(bytesperelem*lines);%ok for complex
  if (cN==filewidth)
    readwholefile=true;
  end;
end
if (readwholefile==true)
  disp(['% Reading whole file: ', infile]);
end;

%%% Check bkformat for complex type: 'cpx*'
if (~ischar(bkformat)) error('FREADBK: bkformat must be string.'); end;
if (~ischar(infile))   error('FREADBK: infile must be string.'); end;
if (length(bkformat)==3)
  if (bkformat=='mph') 
    disp('changing mph format to cpxfloat32');
    bkformat='cpxfloat32';
  end
  if (bkformat=='hgt') 
    error('please use freadhgt for hgt format.');
  end
end
%%% complex types defined as prepended 'cpx'
bkformatOrig = bkformat;
if (length(bkformat)>6)
  if (bkformat(1:3)=='cpx')
    complextype = true;
    bkformat=bkformat(4:length(bkformat));
    bkformatOrig = bkformat;
    if strcmp(bkformat,'float16') % to read half/float16 format
      bkformat = 'uint16';
    end
  end
end

%%% Read from file in column vector.
fid = fopen(infile,'r');
if (fid<0)%						try one more time
  [infile, inpath] = uigetfile('*', 'Select binary input file', 0,0);
  infile = [inpath,infile];
  fid    = fopen(infile,'r');
  if (fid<0) error(ferror(fid)); end;
end;

if (readwholefile==true)
  [data,count]=fread(fid,bkformat);%		read data in column vector
				%		count is number of elements, not bytes...
else
  bytesperelem = datatypesize(bkformat);%	in bytes
  filesize     = fsize(fid);%			in bytes
  filewidth    = filesize/(bytesperelem*lines);%ok for complex
  if (r0 > rN) error('r0 > rN'); end;
  if (c0 > cN) error('c0 > cN'); end;
  if (c0==0) c0=1; cN=filewidth; end;
  if (r0==0) r0=1; rN=lines; end;
  if (r0 < 1)            error('r0 < 1'); end;
  if (rN > lines)        error('rN > lines'); end;
  if (c0 < 1)            error('c0 < 1'); end;
  if (cN > filewidth)    error('cN > width'); end;
  if (~isint(filewidth)) error('numlines file seems to be wrong(?)'); end;

  count   = 0;
  offset  = bytesperelem*(c0-1);%		in bytes
  mywidth = cN-c0+1;%				number of elems 2b read
  lines   = rN-r0+1;%				number of lines in matrix
  if (complextype==true)
    mywidth=mywidth*2;%				correction for pix interleaved complex
    offset =offset*2;%				correction for pix interleaved complex
  end;
  
  data = NaN(mywidth*lines,1); %         added by FvL          
  lnr = 0; %                                    added by FvL  

  for ii=r0-1:rN-1
    lnr = lnr + 1;%                             added by FvL  
    start   = ii*filewidth*bytesperelem+offset;
    status  = fseek(fid,start,'bof');%			beginning of file
    if (status ~= 0)
      disp(['r0 rN c0 cN: ',num2str(r0),' ',num2str(rN),...
			' ',num2str(c0),' ',num2str(cN)]);
      disp(['ii offset start: ',...
	     num2str(ii),' ',num2str(offset),' ',num2str(start)]);
      error(ferror(fid));
    end;
    [d cnt] = fread(fid,mywidth,bkformat);%		column vector
    count   = count+cnt;
    data((lnr-1)*mywidth+1:lnr*mywidth) = d;%   changed by FvL
  end;
end;
fclose(fid);

if (length(data)~=count) %                      added by FvL
  error('Something went wrong while reading the data.')
end


if strcmp(bkformatOrig,'float16')

  data = uint16(data); %FvL was double first, transformation needed for NaN initiation
  if verLessThan('matlab','9.5')
    data = halfprecision(data,'single');
  else
    data = half.typecast(data);
    data = single(data);
  end

else

  data = single(data);

end

if (complextype==true)
  % should be possible by something like, but for speed makes no big diff.
  %data = reshape(data,width,2*lines);%		transpose band interleaved if complex
  %data = complex(data(:,1:lines),data(:,lines+1:2*lines)).';
  data = complex(data(1:2:count),data(2:2:count));
  count=count/2;%					correction for complex types
end;

width = count/lines;%
data  = reshape(data,width,lines).';


%%% Export to type if requested 
% data = int8(data) ...


%%% EOF

