function [Narcs] = ps_form_network(Npsc,Npsc_selections,Nifgs,max_arc_length,network_method,Ncon,Nparts,z)

% Function to determine differential phases between PSCs, based on
% Delaunay triangulation
%
% Input:  - Npsc              number of psc
%         - Npsc_selections   number of psc selections
%         - Nifgs             number of interferograms
%         - max_arc_length    atmospheric correlation length
%         - network_method    network identifier, 'delaunay' or
%                             'spider'
%         - Ncon              minimum number of connections (arcs)
%                             to a psc (network 2 only)
%         - Nparts            number of partitions of a full cycle 
%                             to which the arcs are divided
%                             (network 2 only)
%         - z                 current psc selection
%
% Output: - Narcs            number of arcs
% 
% ----------------------------------------------------------------------
% File............: ps_form_network.m
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
% v1.7.2.9, Freek van Leijen
% - added unwrap2_istep
%



% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global az_spacing r_spacing project_id results_id unwrap_istep unwrap2_istep

Narcs = NaN(Npsc_selections,1);

fid_res = fopen([project_id '_resfile.txt'],'a');

  
% ----------------------------------------------------------------------
% Read data
% ----------------------------------------------------------------------

psc_fid = fopen([project_id '_psc_sel' num2str(z) '.raw'],'r'); % file with temporary psp info
psc_data = fread(psc_fid,[2*Nifgs+7 Npsc(z)],'double')';
fclose(psc_fid);

psc_array = psc_data(:,1:2);
psc_az_orig = psc_data(:,5);
psc_r_orig = psc_data(:,6);
psc_index = find(psc_array(:,2)~=0);
psc_data = psc_data(psc_index,:);

psc_grid_az = psc_data(:,3);
psc_grid_r = psc_data(:,4);
psc_az = psc_data(:,5);
psc_r = psc_data(:,6);
psc_phase = psc_data(:,7:Nifgs+6);
psc_h2ph = psc_data(:,Nifgs+7:2*Nifgs+6);
clear psc_data

switch network_method
  case 'delaunay'
    
    dpsc_phase = NaN(3*Npsc(z)-3,Nifgs);
    dpsc_h2ph = NaN(3*Npsc(z)-3,Nifgs);
    dpsc_arcs = NaN(3*Npsc(z)-3,2);
    % based on theoretical maximum number of arcs for Delaunay
    
    
    
    % ----------------------------------------------------------------------
    % Triangulation by Delaunay
    % ----------------------------------------------------------------------
    
    trian = delaunay(psc_az,(r_spacing/az_spacing)*psc_r);
    % to get a 'nice' triangulation
    trian = sort(trian,2); 
    % operations to extract the unique arcs from the Delaunay triangulation
    A = trian(:,1:2);
    B = trian(:,2:3);
    C = trian(:,1:2:3);
    a = union(A,B,'rows');
    dpsc_arcs_orig = union(a,C,'rows');
    Narcs_orig = size(dpsc_arcs_orig,1);
    
    
    %%% The number of sides can be checked theoretically by: 
    %%% Nsides = 3*Nps - 3 - k, where k is the number of points
    %%% on the so-called convex hull boundary (outside boundary of ps set). So
    %%%
    %%% K = convhull(psc_az,psc_r)
    %%% k = length(K)-1
    %%% theo_Narcs = 3*length(psc_az)-3-k
    %%%
    %%% if (Narcs ~= theo_Narcs)
    %%%   error('Something went wrong in the triangulation process')
    %%% end
    %%% 
    %%% However, because integer pixel cn are used, sometimes a ps lies 
    %%% exactly between two other ps, and 'convhull' makes numerical 
    %%% mistakes. 'delaunay' seems to work correct.
    
  case 'spider'
    dpsc_arcs_orig = ps_construct_redundant_network(psc_az,psc_r,Ncon,Nparts,max_arc_length);
    Narcs_orig = size(dpsc_arcs_orig,1);
    dpsc_phase = NaN(Narcs_orig,Nifgs);
    dpsc_h2ph = NaN(Narcs_orig,Nifgs);
    
  otherwise
    error('You specified a wrong network_method.');
end

% ----------------------------------------------------------------------
% Select phase differences for arcs shorter than atmospheric correlation length
% ----------------------------------------------------------------------

count = 0;
for v = 1:Narcs_orig
  dist = sqrt((az_spacing*(psc_az(dpsc_arcs_orig(v,1)) - ...
                           psc_az(dpsc_arcs_orig(v,2))))^2 + ...
              (r_spacing*(psc_r(dpsc_arcs_orig(v,1))-psc_r(dpsc_arcs_orig(v,2))))^2);
  
  if (dist < max_arc_length)
    count = count + 1;
    dpsc_phase(count,:) = psc_phase(dpsc_arcs_orig(v,2),:)- ...
        psc_phase(dpsc_arcs_orig(v,1),:); 
    dpsc_h2ph(count,:) = (psc_h2ph(dpsc_arcs_orig(v,2),:)+ ...
                          psc_h2ph(dpsc_arcs_orig(v,1),:))/2; % take the mean value of the arc
    dpsc_arcs(count,:) = psc_index(dpsc_arcs_orig(v,:))'; %back to
                                                          %original index
  end
end

Narcs(z) = count;



% ----------------------------------------------------------------------
% Remove NaN lines
% ----------------------------------------------------------------------

dpsc_phase = dpsc_phase(1:Narcs(z),:);
dpsc_h2ph = dpsc_h2ph(1:Narcs(z),:);
dpsc_arcs = dpsc_arcs(1:Narcs(z),:);



% ----------------------------------------------------------------------
% Output to screen
% ----------------------------------------------------------------------

fprintf(1,'In total, %4.0f arcs have been formed for selection %2.0f\n',Narcs(z),z);

fprintf(fid_res,[datestr(now) ', number of arcs in selection %g: %g\n'],z,Narcs(z));



% ----------------------------------------------------------------------
% Write results to file
% ----------------------------------------------------------------------

dpsc_total = [dpsc_arcs dpsc_phase dpsc_h2ph];
dpsc_fid = fopen([project_id '_dpsc_sel' num2str(z) '.raw'],'w');
fwrite(dpsc_fid,dpsc_total','double');
fclose(dpsc_fid);

save([project_id '_network_psc_sel' num2str(z) '_' results_id ...
      '_unwrap_istep' num2str(unwrap_istep) ...
      '_unwrap2_istep' num2str(unwrap2_istep) '.mat'],'psc_r_orig',...
     'psc_az_orig','dpsc_arcs')


fclose(fid_res);

