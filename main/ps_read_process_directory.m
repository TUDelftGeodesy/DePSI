function [orbitnr,dates,filenames_slc,filenames_ifgs,filenames_h2ph,filenames_output,Btemp,Bdop,nSlc,nIfgs,masterIdx,crop,cropIn,cropFinal,Nlines,Npixels,slc_selection,breakpoint,breakpoint2] = ps_read_process_directory(processDir,startDate,stopDate,excludeDate,ifgsVersion,altimg,crop,run_mode,slc_selection_input,breakpoint,breakpoint2,sensor,processor,master,swath_burst)

% Read process directory
%
% Input:    - processDir              process directory
%           - startDate               start date
%           - stopDate                stop date
%           - excludeDate             dates to be excluded
%           - ifgsVersion             interferogram version
%           - altimg                  optional: alternative polarization
%           - crop                    crop
%           - run_mode                run mode
%           - slc_selection_input     slc selection input
%           - breakpoint              breakpoint
%           - breakpoint2             breakpoint 2
%
% Output: 
%
% ----------------------------------------------------------------------
% File............: ps_read_process_directory.m
% Version & Date..: 2.0.1.0, 04-JAN-2016
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

%% Get image list
[slcs,startDate,stopDate] = get_image_list(processDir,startDate,stopDate,excludeDate,sensor);
nSlc = size(slcs,1);
slcsNum = datenum(char(slcs),'yyyymmdd');

%% Get stack parameters
[cropIn,cropFinal,crop,Nlines,Npixels,Bdop] = get_stack_parameters(processDir,slcs,crop,altimg,sensor,processor,master,swath_burst);

orbitnr = slcs;
dates = datestr(datenum(slcs,'yyyymmdd'),'dd-mmm-yyyy');

if strcmp(sensor,'s1')
  if strcmp(processor,'doris_rippl')
    if isempty(swath_burst) % merged image
      filenames_slc = fullfile(processDir,slcs,['phase_corrected_' altimg '@radar.raw']);
      filenames_ifgs = [];
      filenames_h2ph = fullfile(processDir,slcs,'height_to_phase@radar.raw');
      filenames_output = fullfile(processDir,slcs);
    else % single burst
      filenames_slc = fullfile(processDir,slcs,['slice_' num2str(swath_burst(2)) '_swath_' num2str(swath_burst(1))], ...
                        ['phase_corrected_' altimg '@radar.raw']);
      filenames_ifgs = [];
      filenames_h2ph = fullfile(processDir,slcs,['slice_' num2str(swath_burst(2)) '_swath_' num2str(swath_burst(1))], ...
                         'height_to_phase@radar.raw');
      filenames_output = fullfile(processDir,slcs,['slice_' num2str(swath_burst(2)) '_swath_' num2str(swath_burst(1))]);
    end
  else
    if isempty(swath_burst) % merged image
      filenames_slc = fullfile(processDir,slcs,['slave_rsmp' altimg '.raw']);
      filenames_ifgs = fullfile(processDir,slcs,['cint' ifgsVersion altimg '.raw']);
      filenames_h2ph = fullfile(processDir,slcs,['h2ph' ifgsVersion altimg '.raw']);
      filenames_output = fullfile(processDir,slcs);
    else % single burst
      filenames_slc = fullfile(processDir,slcs,['swath_' num2str(swath_burst(1))],['burst_' num2str(swath_burst(2))], ...
                        ['slave_rsmp' altimg '.raw']);
      filenames_ifgs = fullfile(processDir,slcs,['swath_' num2str(swath_burst(1))],['burst_' num2str(swath_burst(2))], ...
                        ['cint' ifgsVersion altimg '.raw']);
      filenames_h2ph = fullfile(processDir,slcs,['swath_' num2str(swath_burst(1))],['burst_' num2str(swath_burst(2))], ...
                        ['h2ph' ifgsVersion altimg '.raw']);
      filenames_output = fullfile(processDir,slcs,['swath_' num2str(swath_burst(1))],['burst_' num2str(swath_burst(2))]);
    end
  end
else
  filenames_slc = fullfile(processDir,slcs,['slave_rsmp' altimg '.raw']);
  filenames_ifgs = fullfile(processDir,slcs,['cint' ifgsVersion altimg '.raw']);
  filenames_h2ph = fullfile(processDir,slcs,['h2ph' ifgsVersion altimg '.raw']);
  filenames_output = fullfile(processDir,slcs);
end
  
masterIdx = find(Bdop==0);
if length(masterIdx)>1
  error('The master image is not uniquely determined. We have to find another method.');
end
datenums = datenum(dates);
Btemp = datenums-datenums(masterIdx);

[Btemp,Btemp_index] = sort(Btemp);
orbitnr = orbitnr(Btemp_index,:);
dates = dates(Btemp_index,:);
filenames_slc = filenames_slc(Btemp_index);
if ~isempty(filenames_ifgs)
  filenames_ifgs = filenames_ifgs(Btemp_index);
end
filenames_h2ph = filenames_h2ph(Btemp_index);
filenames_output = filenames_output(Btemp_index);
Bdop = Bdop(Btemp_index);
cropIn = cropIn(Btemp_index,:);
  
%find master index
masterIdx = find(Btemp==0);
filenames_slc = [filenames_slc;filenames_slc(masterIdx)];
orbitnr = [orbitnr;orbitnr(masterIdx,:)];
dates = [dates;dates(masterIdx,:)];
cropIn = [cropIn;cropIn(masterIdx,:)];
  
Btemp(masterIdx) = [];
orbitnr(masterIdx,:) = [];
dates(masterIdx,:) = [];
filenames_slc(masterIdx) = [];
if ~isempty(filenames_ifgs)
  filenames_ifgs(masterIdx) = [];
end
filenames_h2ph(masterIdx) = [];
filenames_output(masterIdx) = [];
Bdop(masterIdx) = [];
cropIn(masterIdx,:) = [];

nSlc = size(filenames_slc,1);
%nIfgs = size(filenames_ifgs,1);
nIfgs = nSlc-1;

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
  

