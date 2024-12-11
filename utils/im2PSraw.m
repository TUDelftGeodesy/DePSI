function im2PSraw( inIMG, outBin)
%function im2PSraw(inIMG, outBin)
%
% Convert mask.tif (or matlab compatible image file) 
% to mask.raw (uint8) binary file
%
% Usage: im2PSraw <input.tif> <output.raw>
%        im2PSraw zh_mrm_mask.tif zh_mrm_mask.raw
%
% Mahmut Arikan <m.arikan@tudelft.nl> 20060711
% adapted 20071107, FvL (int32 -> uint8)

if nargin<2 || ~exist(inIMG)
   error('Two Arguments Required. or Check path to mask file')
end

% read in image to matrix named mask
mask=imread(inIMG);

% Set output binary file format to unsigned integer*8
type='uint8';

% Main
if ~exist(outBin)
    binwrite(mask, outBin, type);
    disp(sprintf('\nMask is saved as %s having %ix%i dimensions \n', ...
        outBin,size(mask,1),size(mask,2)));    
else
    reply = input(sprintf('Do you want overwrite existing %s? Y/N : ', outBin), 's');
    switch lower(reply)
      case 'y'
        binwrite(mask, outBin, type);
        disp(sprintf('\nMask is saved as %s having %ix%i dimensions \n', ...
        outBin,size(mask,1),size(mask,2)));
      case 'n'
        outBin=input('New file.raw: ', 's');
        binwrite(mask, outBin, type);
        disp(sprintf('\nMask is saved as %s having %ix%i dimensions \n', ...
        outBin,size(mask,1),size(mask,2)));
      otherwise
        disp('Got you! Next time..');
        exit;
    end
end

% generic subfunc for writing binary file
function binwrite(inMat, outBin, type)
    fid = fopen(outBin, 'w');
    fwrite(fid, inMat', type);
    fclose(fid);