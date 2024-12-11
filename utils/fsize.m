function numbytes = fsize(filename);
% FSIZE  return size of file in bytes.
%   b = FSIZE(FILENAME) return size in bytes of specified file.
%   b = FSIZE(FID)      return size in bytes of specified fid (file).
%   -1 is returned if specified file does not exist. 
%   (I found sometimes problems when file already was opened?)
%

% $Revision: 1.5 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 08-Mar-2000

% check if input is filename or fid
if (ischar(filename)) 
  tmpfid=fopen(filename,'r');
else
  tmpfid=filename;
  if (tmpfid<0)
    numbytes=-1;
  else
    oldpos=ftell(tmpfid);
  end
end

if (tmpfid<0)
  numbytes=-1;
else
  status=fseek(tmpfid,0,'eof');
  if (status==-1) error(ferror(tmpfid)); end;
  numbytes=ftell(tmpfid);
  if (numbytes==-1) error(ferror(tmpfid)); end;
  if (ischar(filename))
    fclose(tmpfid);%			close file if opened
  else
    fseek(tmpfid,oldpos,'bof');%	or reset file pointer
  end;
end;

%%% EOF
