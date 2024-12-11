function ps_make_mrm(param_file)

%
% ----------------------------------------------------------------------
% File............: ps_make_mrm.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Authors.........: Freek van Leijen
%                   Delft Institute of Earth Observation and Space Systems
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
% v1.7.2.16, Freek van Leijen
% - filenames in cells
%


global max_mem_buffer Nlines Npixels az_spacing r_spacing project_id

ps_readinput_parameters; % script to read input parameters

[orbitnr,dates,filenames_slc,filenames_ifgs,filenames_h2ph, ...
 filenames_output,Bperp,Btemp,Bdop,l0,lN,p0,pN] = ...
  textread(input_file,'%s %s %s %s %s %s %f %f %f %f %f %f %f');
crop_in = [l0 lN p0 pN];
orbitnr = strvcat(orbitnr);
dates = strvcat(dates);

if strcmp(cellstr(filenames_ifgs(end)),'dummy')
  filenames_ifgs(end) = [];
  filenames_h2ph(end) = [];
  filenames_output(end) = [];
  Bperp(end) = [];
  Btemp(end) = [];
  Bdop(end) = [];
end

Nslc = size(filenames_slc,1);
Nifgs = size(filenames_ifgs,1);

[Btemp,Btemp_index] = sort(Btemp);
orbitnr(1:Nifgs,:) = orbitnr(Btemp_index,:);
dates(1:Nifgs,:) = dates(Btemp_index,:);
filenames_slc(1:Nifgs) = filenames_slc(Btemp_index);
filenames_ifgs = filenames_ifgs(Btemp_index);
filenames_h2ph = filenames_h2ph(Btemp_index);
filenames_output = filenames_output(Btemp_index);
Bperp = Bperp(Btemp_index);
Bdop = Bdop(Btemp_index);
crop_in(1:Nifgs,:) = crop_in(Btemp_index,:);

if max(abs(Btemp))>33 % to make sure Btemp is in years
  Btemp = Btemp/365;
end

if ischar(slc_selection_input)
  slc_selection = ps_determine_slc_selection(orbitnr,slc_selection_input);
elseif isempty(slc_selection_input)
  slc_selection = 1:Nslc;
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

crop_final = [max(crop_in(:,1)) min(crop_in(:,2)) max(crop_in(:,3)) min(crop_in(:,4))];

if ~isempty(crop)
  crop = [max(crop(1),crop_final(1)) min(crop(2),crop_final(2)) ...
          max(crop(3),crop_final(3)) min(crop(4),crop_final(4))];
else
  crop = crop_final;
end
Nlines = crop(2)-crop(1)+1;
Npixels = crop(4)-crop(3)+1;

[calfactors] = ps_calibration(filenames_slc,...
                              filenames_ifgs,...
                              filename_water_mask,...
                              psc_selection_gridsize,...
                              slc_selection,...
                              crop,...
                              crop_in);

ps_mrm(filenames_slc,...
       slc_selection,...
       psc_selection_gridsize,...
       calfactors,...
       crop,...
       crop_in);

save([project_id '_project.mat']);
  
