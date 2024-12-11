function [orbitnr,dates,filenames_slc,filenames_ifgs,filenames_h2ph,filenames_output,Btemp,Bdop,nSlc,nIfgs,masterIdx,crop,cropIn,cropFinal,Nlines,Npixels,slc_selection,breakpoint,breakpoint2] = ps_read_input_file(input_file,crop,run_mode,slc_selection_input,breakpoint,breakpoint2)

% Read input file
%
% Input:    - input_file              input file
%           - crop                    crop
%           - run_mode                run mode
%           - slc_selection_input     slc selection input
%           - breakpoint              breakpoint
%           - breakpoint2             breakpoint 2
%
% Output: 
%
% ----------------------------------------------------------------------
% File............: ps_read_input_file.m
% Version & Date..: 2.0.1.0, 4-JAN-2016
% Authors.........: Freek van Leijen
%                   Delft University of Technology
% ----------------------------------------------------------------------
%
% This software is developed by Delft University of Technology and is
% intended for scientific use only. Applications for commercial use are
% prohibited.
%
% Copyright (c) 2004-2009 Delft University of Technology, The Netherlands
%
% Change log
%

[orbitnr,dates,filenames_slc,filenames_ifgs,filenames_h2ph, ...
 filenames_output,Bperp,Btemp,Bdop,l0,lN,p0,pN] = ...
  textread(input_file,'%s %s %s %s %s %s %f %f %f %f %f %f %f');
cropIn = [l0 lN p0 pN];
orbitnr = strvcat(orbitnr);
dates = strvcat(dates);
  
if strcmp(filenames_ifgs(end),'dummy')
  filenames_ifgs(end) = [];
  filenames_h2ph(end) = [];
  filenames_output(end) = [];
  Bperp(end) = [];
  Btemp(end) = [];
  Bdop(end) = [];
end
  
nSlc = size(filenames_slc,1);
nIfgs = size(filenames_ifgs,1);
 
[Btemp,Btemp_index] = sort(Btemp);
orbitnr(1:nIfgs,:) = orbitnr(Btemp_index,:);
dates(1:nIfgs,:) = dates(Btemp_index,:);
filenames_slc(1:nIfgs) = filenames_slc(Btemp_index);
filenames_ifgs = filenames_ifgs(Btemp_index);
filenames_h2ph = filenames_h2ph(Btemp_index);
filenames_output = filenames_output(Btemp_index);
Bperp = Bperp(Btemp_index);
Bdop = Bdop(Btemp_index);
cropIn(1:nIfgs,:) = cropIn(Btemp_index,:);
  
%find master index
masterIdx = find(Btemp<0);
if isempty(masterIdx)
  masterIdx = 1;
else
  masterIdx = masterIdx(end)+1;
end

cropFinal = [max(cropIn(:,1)) min(cropIn(:,2)) max(cropIn(:,3)) min(cropIn(:,4))];

if ~isempty(crop)
  crop = [max(crop(1),cropFinal(1)) min(crop(2),cropFinal(2)) ...
          max(crop(3),cropFinal(3)) min(crop(4),cropFinal(4))];
else
  crop = cropFinal;
end
Nlines = crop(2)-crop(1)+1;
Npixels = crop(4)-crop(3)+1;

switch run_mode
  case 'validation'

    if exist('validation','dir')~=7
      mkdir('validation');
    end
      
    orbitnr_valid = [orbitnr(1:masterIdx-1,:);...
                     orbitnr(end,:);...
                     orbitnr(masterIdx:end-1,:)];
    Btemp_valid = round(365*[Btemp(1:masterIdx-1);...                 
                        0; ...
                        Btemp(masterIdx:end)]);
      
    ifgs_valid = [repmat(masterIdx,nIfgs,1) [1:masterIdx-1 masterIdx+1:nSlc]'];
      
    orbitnr_valid_fid = fopen(['validation/orbitnr_valid.txt'],'w');
    fprintf(orbitnr_valid_fid,'%g\n',str2num(orbitnr_valid));
    fclose(orbitnr_valid_fid);
      
    btemp_valid_fid = fopen(['validation/btemp_valid.txt'],'w');
    fprintf(btemp_valid_fid,'%g\n',Btemp_valid);
    fclose(btemp_valid_fid);

    btemp_ifgs_valid_fid = fopen(['validation/btemp_ifgs_valid.txt'],'w');
    fprintf(btemp_ifgs_valid_fid,'%g\n',round(365*Btemp));
    fclose(btemp_ifgs_valid_fid);
      
    ifgs_valid_fid = fopen(['validation/ifgs_valid.txt'],'w');
    fprintf(ifgs_valid_fid,'%g\t%g\n',ifgs_valid');
    fclose(ifgs_valid_fid);
   
end
  
  
if max(abs(Btemp))>33 % to make sure Btemp is in years
  Btemp = Btemp/365;
end
  
if ischar(slc_selection_input)
  slc_selection = ps_determine_slc_selection(orbitnr,slc_selection_input);
elseif isempty(slc_selection_input)
  slc_selection = 1:nSlc;
else
  slc_selection = slc_selection_input;
end
  
if ~isempty(breakpoint)
  if breakpoint>1000
    breakpoint = find(str2num(orbitnr)==breakpoint);
  end
end

if ~isempty(breakpoint2)
  if breakpoint2>1000
    breakpoint2 = find(str2num(orbitnr)==breakpoint2);
  end
end


