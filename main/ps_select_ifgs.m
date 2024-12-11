function [Nifgs] = ps_select_ifgs(ifg_selection,Npsc,Npsc_selections,Npsp,Nifgs);

% Function to make a subset of the available ifgs for
% further processing.
%
% Input:  - ifg_selection     vector with ifgs indices
%         - Npsc              number of ps candidates
%         - Npsc_selections   number of psc selections
%         - Npsp              number of potential ps
%         - Nifgs             number of interferograms (old)
%
% Output: - Nifgs             number of interferograms (new)
%
% ----------------------------------------------------------------------
% File............: ps_select_ifgs.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Author..........: Freek van Leijen
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
%


% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global project_id max_mem_buffer Npar_max



% ----------------------------------------------------------------------
% Loop
% ----------------------------------------------------------------------

for z = 1:Npsc_selections
 
  if ~exist([project_id '_psc_orig_sel' num2str(z) '.raw']) 
    %step not run before
    
    copyfile([project_id '_psc_sel' num2str(z) '.raw'],...
	     [project_id '_psc_orig_sel' num2str(z) '.raw']);
    copyfile([project_id '_psc_results_sel' num2str(z) '.raw'],...
	     [project_id '_psc_results_orig_sel' num2str(z) '.raw']);
    copyfile([project_id '_psc_atmo_sel' num2str(z) '.raw'],...
	     [project_id '_psc_atmo_orig_sel' num2str(z) '.raw']);
  end
  
  psc_fid_in1 = fopen([project_id '_psc_orig_sel' num2str(z) '.raw'],'r');
  psc_fid_out1 = fopen([project_id '_psc_sel' num2str(z) '.raw'],'w');
  psc_data = fread(psc_fid_in1,[2*Nifgs+7 Npsc(z)],'double')';
  
  psc_data = psc_data(:,[1:6 6+ifg_selection (6+Nifgs)+ifg_selection ...
                      2*Nifgs+7]);
  fwrite(psc_fid_out1,psc_data','double');
  fclose(psc_fid_out1);
  clear psc_data
  
  psc_fid_in2 = fopen([project_id '_psc_results_orig_sel' num2str(z) '.raw'],'r');
  psc_fid_out2 = fopen([project_id '_psc_results_sel' num2str(z) '.raw'],'w');
  psc_data = fread(psc_fid_in2,[3*Nifgs+Npar_max+2 Npsc(z)],'double')';
  
  psc_data = psc_data(:,[1:Npar_max+2 Npar_max+2+ifg_selection Npar_max+2+Nifgs+ifg_selection Npar_max+2+2*Nifgs+ifg_selection]);
  fwrite(psc_fid_out2,psc_data','double');
  fclose(psc_fid_out2);
  clear psc_data

  psc_fid_in3 = fopen([project_id '_psc_atmo_orig_sel' num2str(z) '.raw'],'r');
  psc_fid_out3 = fopen([project_id '_psc_atmo_sel' num2str(z) '.raw'],'w');
  psc_data = fread(psc_fid_in3,[Nifgs Npsc(z)],'double')';
  
  psc_data = psc_data(:,ifg_selection);
  fwrite(psc_fid_out3,psc_data','double');
  fclose(psc_fid_out3);
  clear psc_data

end

Npsp_buffer = floor(max_mem_buffer/(6*Nifgs*8));
if (Npsp_buffer>=Npsp(1))
  Npsp_buffer = Npsp(1);
  Npsp_rem_buffer = 0;
  Nbuffers = 1;
else
  Nbuffers = floor(Npsp(1)/Npsp_buffer);
  Npsp_rem_buffer = rem(Npsp(1),Npsp_buffer);
  if (Npsp_rem_buffer > 0)
    Nbuffers = Nbuffers + 1;
  end
end

for z = 1:Npsc_selections

  if ~exist([project_id '_psp_orig_sel' num2str(z) '.raw'])
    %step not run before

    copyfile([project_id '_psp_sel' num2str(z) '.raw'],...
	     [project_id '_psp_orig_sel' num2str(z) '.raw']);
    copyfile([project_id '_psp_atmo_sel' num2str(z) '.raw'],...
	     [project_id '_psp_atmo_orig_sel' num2str(z) '.raw']);
  end
  
  psp_fid_in1(z) = fopen([project_id '_psp_orig_sel' num2str(z) '.raw'],'r');
  psp_fid_out1(z) = fopen([project_id '_psp_sel' num2str(z) '.raw'],'w');
  psp_fid_in2(z) = fopen([project_id '_psp_atmo_orig_sel' num2str(z) '.raw'],'r');
  psp_fid_out2(z) = fopen([project_id '_psp_atmo_sel' num2str(z) '.raw'],'w');

end

for v = 1:Nbuffers
  if (v==Nbuffers)&(Npsp_rem_buffer~=0)
    Npsp_buffer = Npsp_rem_buffer;
  end
  
  for z = 1:Npsc_selections
    psp_data = fread(psp_fid_in1(z),[2*Nifgs+4 Npsp_buffer],'double')';
    psp_data = psp_data(:,[1:4 4+ifg_selection (Nifgs+4)+ifg_selection]);
    fwrite(psp_fid_out1(z),psp_data','double');
    psp_data = fread(psp_fid_in2(z),[Nifgs Npsp_buffer],'double')';
    psp_data = psp_data(:,ifg_selection);
    fwrite(psp_fid_out2(z),psp_data','double');
  end
end
fclose('all');

Nifgs = length(ifg_selection);
