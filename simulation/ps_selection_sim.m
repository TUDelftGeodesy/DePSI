function [grid_array_az,grid_array_r,Npsc,Npsp] = ps_selection_sim(Nifgs,psc_selection_gridsize,Npsc_selections,psc_selection_method,ifgs,ambi,topo,defo1,defo2,defo3,defo_quad,defo_cubic,defo_nl1,defo_nl2,subpixel,atmo,Bperp,perc)

% Function to read a stack of co-registered, calibrated slc's 
% in buffers and to select persistent scatterer candidates (psc)
% using the amplitude dispersion (ferretti et al., 2001). Also
% potential ps are selected for latter use.
%
% Input:    - Nifgs                   number of interferograms
%           - psc_selection_gridsize  psc selection gridsize
%           - Npsc_selections         number of psc selections
%           - psc_selection_method    psc selection method,
%                                     'threshold' or 'eachgrid'
%           - ifgs                    simulated interferograms (wrapped)
%           - ambi                    simulated ambiguities
%           - topo                    simulated topography
%           - defo1                   simulated linear defo 1
%           - defo2                   simulated linear defo 2
%           - defo_nl1                simulated non-linear defo amplitude
%           - defo_nl2                simulated non-linear defo t0
%           - subpixel                subpixel position (azimuth)
%           - atmo                    atmosphere
%           - Bperp                   perpendicular baselines
%           - perc                    percentage of pixels as psp
%
% Output:   - grid_array_az           array with borders gridcells in
%                                     azimuth direction
%           - grid_array_r            array with borders gridcells in
%                                     range direction
%           - Npsc                    number of psc's
%           - Npsp                    number of psp's
%
% ----------------------------------------------------------------------
% File............: ps_selection_sim.m
% Version & Date..: 1.7.2.8, 19-OCT-2009
% Authors.........: Freek van Leijen
%                   Gini Ketelaar
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

global max_mem_buffer Nlines Npixels az_spacing r_spacing fig
global ps_eval_method detail_plots project_id

Npsc = zeros(Npsc_selections,1);
Npsp = zeros(Npsc_selections,1);

theta = 23*pi/180; % incident angle [deg]
H = 780000;            % satellite vertical height [m]
R = H/cos(theta);      % antenna-target distance [m]
h2ph_orig = Bperp'/(R*sin(theta)); % height to phase factors

mrm = defo1;
mrm_fid = fopen([project_id '_mrm.raw'],'w');
fwrite(mrm_fid,mrm','single'); % transpose is important!
fclose(mrm_fid);

mrm_min = min(mrm(:));
mrm_max = max(mrm(:));
mrm_new = uint8(255*(mrm-mrm_min)./(mrm_max-mrm_min));
mrm_new = imadjust(mrm_new); %stretch

mrm_fid = fopen([project_id '_mrm_uint8.raw'],'w');
fwrite(mrm_fid,mrm_new','uint8');
fclose(mrm_fid);

amp_disp_fid = fopen([project_id '_amp_disp.raw'],'w');
fwrite(amp_disp_fid,rand(Nlines,Npixels)','single'); % transpose is important!
fclose(amp_disp_fid);

sl_mask = rand(Nlines,Npixels);
sl_mask(sl_mask<0.99) = 0;
sl_mask = logical(sl_mask);
side_lobe_fid = fopen([project_id '_side_lobe_mask.raw'],'w');
fwrite(side_lobe_fid,sl_mask','uint8'); % transpose is important!
fclose(side_lobe_fid);


% ----------------------------------------------------------------------
% Open the data files
% ----------------------------------------------------------------------

psc_fid = NaN(Npsc_selections,1);
psc_valid_fid = NaN(Npsc_selections,1);
for z = 1:Npsc_selections
  psc_fid(z) = fopen([project_id '_psc_2orig_sel' num2str(z) '.raw'],'w'); % file with psc info
  psc_valid_fid(z) = fopen([project_id '_psc_valid_sel' num2str(z) '.raw'],'w');
end

switch ps_eval_method
  case 'psp'
    psp_fid = NaN(Npsc_selections,1);
    for z = 1:Npsc_selections
      psp_fid(z) = fopen([project_id '_psp_2orig_sel' num2str(z) '.raw'],'w'); % file with psp info
    end
    psp_valid_fid = fopen([project_id '_psp_valid.raw'],'w');
  case 'whole'
    Npsp = [];
  otherwise
    error('You specified a wrong ps_eval_method.');
end



% ----------------------------------------------------------------------
% Determine buffersize
% ----------------------------------------------------------------------

Nlines_g = round(psc_selection_gridsize/az_spacing);
Npixels_g = round(psc_selection_gridsize/r_spacing);

Ng_az = floor(Nlines/Nlines_g);
Ng_r = floor(Npixels/Npixels_g);
offset_az = floor(rem(Nlines,Nlines_g)/2);
offset_r = floor(rem(Npixels,Npixels_g)/2);

buffer_size_az = floor(max_mem_buffer/(8*8*Nifgs*Npixels*Nlines_g));
% buffer_size_az in number of gridcells in azimuth direction
% 6 because of 2 float array's, 2 complex array's, 1 multi image
% and 1 for the rest, 8 because of double

if (buffer_size_az == 0)
  buffer_size_az = 1;
  buffer_size_r = floor(max_mem_buffer/(8*8*Nifgs*Npixels_g*Nlines_g));
  if (buffer_size_r == 0)
    error('Not enough memory to proceed. Please increase the allocated memory or reduce the grid size');
  end
  Nbuffers_r = floor(Ng_r/buffer_size_r);
  rem_buffer_r = rem(Ng_r,buffer_size_r);
else
  buffer_size_r = Ng_r;
  Nbuffers_r = 1;
  rem_buffer_r = 0;
end

if (buffer_size_az*Nlines_g > Nlines)
  buffer_size_az = Ng_az;
end
   
Nbuffers_az = floor(Ng_az/buffer_size_az);
rem_buffer_az = rem(Ng_az,buffer_size_az);



% ----------------------------------------------------------------------
% Determine borders gridcells
% ----------------------------------------------------------------------

grid_array_az = NaN(Ng_az,3);
grid_array_az(:,1) = [1 offset_az+Nlines_g+1:Nlines_g:offset_az+(Ng_az-1)*Nlines_g+1]';
grid_array_az(:,2) = [offset_az+Nlines_g:Nlines_g:offset_az+ ...
                    (Ng_az-1)*Nlines_g Nlines]';
for v = 1:Nbuffers_az
  buffer = (v-1)*buffer_size_az+1:v*buffer_size_az;
  grid_array_az(buffer,3) = repmat(v,length(buffer),1);
end

if (rem_buffer_az ~= 0)
  buffer = v*buffer_size_az+1:Ng_az;
  grid_array_az(buffer,3) = repmat(v+1,length(buffer),1);
  Nbuffers_az = Nbuffers_az + 1; % add remaining buffer
end

grid_array_r = NaN(Ng_r,3);
grid_array_r(:,1) = [1 offset_r+Npixels_g+1:Npixels_g:offset_r+(Ng_r-1)*Npixels_g+1]';
grid_array_r(:,2) = [offset_r+Npixels_g:Npixels_g:offset_r+ ...
                    (Ng_r-1)*Npixels_g Npixels]';

for v = 1:Nbuffers_r
  buffer = (v-1)*buffer_size_r+1:v*buffer_size_r;
  grid_array_r(buffer,3) = repmat(v,length(buffer),1);
end

if (rem_buffer_r ~= 0)
  buffer = v*buffer_size_r+1:Ng_r;
  grid_array_r(buffer,3) = repmat(v+1,length(buffer),1);
  Nbuffers_r = Nbuffers_r + 1; % add remaining buffer
end

Npsp_grid = floor(max(2,(perc/100)*(size(topo,1)^2)/(Nbuffers_az*Nbuffers_r)));
psc_h2ph_grid = repmat(h2ph_orig,Npsc_selections,1);
psp_h2ph_grid = repmat(h2ph_orig,Npsp_grid,1);
psc_amp_disp_grid = NaN(Npsc_selections,1);

% ----------------------------------------------------------------------
% Main loop per buffer
% ----------------------------------------------------------------------

for v = 1:Nbuffers_az

  
  % ----------------------------------------------------------------------
  % Read data
  % ----------------------------------------------------------------------

  index_az = find(grid_array_az(:,3)==v);
  Ngrid_az = length(index_az);
  begin_buffer_az = grid_array_az(index_az(1),1);
  end_buffer_az = grid_array_az(index_az(end),2);
  buffer_size_az = end_buffer_az-begin_buffer_az+1;

  for w = 1:Nbuffers_r

    index_r = find(grid_array_r(:,3)==w);
    Ngrid_r = length(index_r);
    begin_buffer_r = grid_array_r(index_r(1),1);
    end_buffer_r = grid_array_r(index_r(end),2);
    buffer_size_r = end_buffer_r-begin_buffer_r+1;

    
    % ----------------------------------------------------------------------
    % Select the psc and psp
    % ----------------------------------------------------------------------
    
    for g = 1:Ngrid_az
      for h = 1:Ngrid_r
        
        Nl = length([grid_array_az(index_az(g),1)-begin_buffer_az+1:grid_array_az(index_az(g),2)-begin_buffer_az+1]);
	Np = length([grid_array_r(index_r(h),1)-begin_buffer_r+1:grid_array_r(index_r(h),2)-begin_buffer_r+1]);
        
        Ntot = Np*Nl;
        
        index = unique(ceil(rand(Npsc_selections,1)*Ntot));
        index = [index; ones(Npsc_selections-size(index,1),1)];
        [az,r] = ind2sub([Nl Np],index);
				
        index2 = round(randn(1));

        psc_phase_grid = NaN(Npsc_selections,Nifgs);
        psc_ambi_grid = NaN(Npsc_selections,Nifgs);
        psc_az_grid = NaN(Npsc_selections,1);
        psc_r_grid = NaN(Npsc_selections,1);
        psc_topo_grid = NaN(Npsc_selections,1);
        psc_defo1_grid = NaN(Npsc_selections,1);
        psc_defo2_grid = NaN(Npsc_selections,1);
        psc_defo3_grid = NaN(Npsc_selections,1);
        psc_defo_quad = NaN(Npsc_selections,1);
	psc_defo_cubic = NaN(Npsc_selections,1);
        psc_defo_nl1_grid = NaN(Npsc_selections,1);
        psc_defo_nl2_grid = NaN(Npsc_selections,1);
        psc_subpixel_grid = NaN(Npsc_selections,1);
        psc_m_atmo_grid = NaN(Npsc_selections,1);
        psc_s_atmo_grid = NaN(Npsc_selections,Nifgs);
        Npsc_grid = zeros(Npsc_selections,1);
    
        switch psc_selection_method
          case 'threshold'
            if index2==1
                for z = 1:Npsc_selections
                    psc_az_grid(z) = grid_array_az(index_az(g),1)-begin_buffer_az+az(z);
                    psc_r_grid(z) = grid_array_r(index_r(h),1)-begin_buffer_r+r(z);
                    psc_phase_grid(z,:) = ifgs(psc_az_grid(z),psc_r_grid(z),:);
                    psc_ambi_grid(z,:) = ambi(psc_az_grid(z),psc_r_grid(z),:);
                    psc_topo_grid(z) = topo(psc_az_grid(z),psc_r_grid(z));
                    psc_defo1_grid(z) = defo1(psc_az_grid(z),psc_r_grid(z));
                    psc_defo2_grid(z) = defo2(psc_az_grid(z),psc_r_grid(z));
                    psc_defo3_grid(z) = defo3(psc_az_grid(z),psc_r_grid(z));
                    psc_defo_quad(z) = defo_quad(psc_az_grid(z),psc_r_grid(z));
                    psc_defo_cubic(z) = defo_cubic(psc_az_grid(z),psc_r_grid(z));
                    psc_defo_nl1_grid(z) = defo_nl1(psc_az_grid(z),psc_r_grid(z));
                    psc_defo_nl2_grid(z) = defo_nl2(psc_az_grid(z),psc_r_grid(z));
                    psc_subpixel_grid(z) = subpixel(psc_az_grid(z),psc_r_grid(z));
                    psc_m_atmo_grid(z) = atmo(psc_az_grid(z),psc_r_grid(z),1);
                    psc_s_atmo_grid(z,:) = reshape(atmo(psc_az_grid(z),psc_r_grid(z),2:end),1,Nifgs);
                    Npsc_grid(z) = Npsc_grid(z)+1;
                end
                % in this loop we can built something to avoid sidelobes
            end
          
          case 'eachgrid'
            for z = 1:Npsc_selections
                psc_az_grid(z) = grid_array_az(index_az(g),1)-begin_buffer_az+az(z);
                psc_r_grid(z) = grid_array_r(index_r(h),1)-begin_buffer_r+r(z);
                psc_phase_grid(z,:) = ifgs(psc_az_grid(z),psc_r_grid(z),:);
                psc_ambi_grid(z,:) = ambi(psc_az_grid(z),psc_r_grid(z),:);
                psc_topo_grid(z) = topo(psc_az_grid(z),psc_r_grid(z));
                psc_defo1_grid(z) = defo1(psc_az_grid(z),psc_r_grid(z));
                psc_defo2_grid(z) = defo2(psc_az_grid(z),psc_r_grid(z));
		psc_defo3_grid(z) = defo3(psc_az_grid(z),psc_r_grid(z));
		psc_defo_quad(z) = defo_quad(psc_az_grid(z),psc_r_grid(z));
		psc_defo_cubic(z) = defo_cubic(psc_az_grid(z),psc_r_grid(z));
                psc_defo_nl1_grid(z) = defo_nl1(psc_az_grid(z),psc_r_grid(z));
                psc_defo_nl2_grid(z) = defo_nl2(psc_az_grid(z),psc_r_grid(z));
                psc_subpixel_grid(z) = subpixel(psc_az_grid(z),psc_r_grid(z));
                psc_m_atmo_grid(z) = atmo(psc_az_grid(z),psc_r_grid(z),1);
                psc_s_atmo_grid(z,:) = reshape(atmo(psc_az_grid(z),psc_r_grid(z),2:end),1,Nifgs);
                Npsc_grid(z) = Npsc_grid(z)+1;
            end
            % in this loop we can built something to avoid sidelobes
        
        otherwise
          error('Error: the psc selection method was not specified correctly');
        end
        
        if max(Npsc_grid) ~= 0
          for z = 1:Npsc_selections
            psc_array = [(Npsc(z)+1:Npsc(z)+Npsc_grid(z))' ones(Npsc_grid(z),1)];
            psc_total = [psc_array index_az(g) index_r(h) psc_az_grid(z) ...
                         psc_r_grid(z) psc_phase_grid(z,:) psc_h2ph_grid(z,:) ...
                         psc_amp_disp_grid(z)];
            fwrite(psc_fid(z),psc_total','double');
            psc_validate = [psc_az_grid(z) psc_r_grid(z) psc_phase_grid(z,:) ...
			    psc_ambi_grid(z,:) psc_s_atmo_grid(z,:) ...
			    psc_topo_grid(z) psc_m_atmo_grid(z) ...
			    psc_subpixel_grid(z) psc_defo1_grid(z) ...
			    psc_defo2_grid(z) psc_defo3_grid(z) psc_defo_quad(z) ...
			    psc_defo_cubic(z) psc_defo_nl1_grid(z) ...
                            psc_defo_nl2_grid(z)];
            fwrite(psc_valid_fid(z),psc_validate','double');
          end
        end


        switch ps_eval_method
          case 'psp'
            index = unique(ceil(rand(Npsp_grid,1)*Ntot));
            index = [index; ones(Npsp_grid-size(index,1),1)];
            [az,r] = ind2sub([Nl Np],index);
            
            psp_phase_grid = NaN(Npsp_grid,Nifgs);
            psp_ambi_grid = NaN(Npsp_grid,Nifgs);
            psp_az_grid = NaN(Npsp_grid,1);
            psp_r_grid = NaN(Npsp_grid,1);
            psp_topo_grid = NaN(Npsp_grid,1);
            psp_defo1_grid = NaN(Npsp_grid,1);
            psp_defo2_grid = NaN(Npsp_grid,1);
            psp_defo3_grid = NaN(Npsp_grid,1);
            psp_defo_quad_grid = NaN(Npsp_grid,1);
            psp_defo_cubic_grid = NaN(Npsp_grid,1);
            psp_defo_nl1_grid = NaN(Npsp_grid,1);
            psp_defo_nl2_grid = NaN(Npsp_grid,1);
            psp_subpixel_grid = NaN(Npsp_grid,1);
            psp_m_atmo_grid = NaN(Npsp_grid,1);
            psp_s_atmo_grid = NaN(Npsp_grid,Nifgs);
            
            for z = 1:Npsp_grid
                psp_az_grid(z) = grid_array_az(index_az(g),1)-begin_buffer_az+az(z);
                psp_r_grid(z) = grid_array_r(index_r(h),1)-begin_buffer_r+r(z);
                psp_phase_grid(z,:) = ifgs(psp_az_grid(z),psp_r_grid(z),:);
                psp_ambi_grid(z,:) = ambi(psp_az_grid(z),psp_r_grid(z),:);
                psp_topo_grid(z) = topo(psp_az_grid(z),psp_r_grid(z));
                psp_defo1_grid(z) = defo1(psp_az_grid(z),psp_r_grid(z));
                psp_defo2_grid(z) = defo2(psp_az_grid(z),psp_r_grid(z));
                psp_defo3_grid(z) = defo3(psp_az_grid(z),psp_r_grid(z));
                psp_defo_quad_grid(z) = defo_quad(psp_az_grid(z),psp_r_grid(z));
                psp_defo_cubic_grid(z) = defo_cubic(psp_az_grid(z),psp_r_grid(z));
                psp_defo_nl1_grid(z) = defo_nl1(psp_az_grid(z),psp_r_grid(z));
                psp_defo_nl2_grid(z) = defo_nl2(psp_az_grid(z),psp_r_grid(z));
                psp_subpixel_grid(z) = subpixel(psp_az_grid(z),psp_r_grid(z));
                psp_m_atmo_grid(z) = atmo(psp_az_grid(z),psp_r_grid(z),1);
                psp_s_atmo_grid(z,:) = reshape(atmo(psp_az_grid(z),psp_r_grid(z),2:end),1,Nifgs);
            end
            
            if Npsp_grid ~= 0
                psp_total = [repmat(index_az(g),Npsp_grid,1) repmat(index_r(h),Npsp_grid,1) psp_az_grid psp_r_grid psp_phase_grid psp_h2ph_grid];
                psp_validate = [psp_az_grid psp_r_grid psp_phase_grid ...
                                psp_ambi_grid psp_s_atmo_grid ...
				psp_topo_grid psp_m_atmo_grid ...
				psp_subpixel_grid psp_defo1_grid ...
                                psp_defo2_grid psp_defo3_grid psp_defo_quad_grid ...
				psp_defo_cubic_grid psp_defo_nl1_grid ...
                                psp_defo_nl2_grid];
                
                for z = 1:Npsc_selections
                    fwrite(psp_fid(z),psp_total','double');
                end
                fwrite(psp_valid_fid,psp_validate','double');
            end

            Npsp = Npsp + repmat(Npsp_grid,Npsc_selections,1);
        end
        Npsc = Npsc + Npsc_grid;
      end %for
    end %for
  end %for
end %for



% ----------------------------------------------------------------------
% Output
% ----------------------------------------------------------------------

for z = 1:Npsc_selections
  psc_fid = fopen([project_id '_psc_2orig_sel' num2str(z) '.raw'],'r'); % file with temporary psp info
  psc_data = fread(psc_fid,[2*Nifgs+7 Npsc(z)],'double')';
  fclose(psc_fid);
  
  psc_az_selection = psc_data(:,5);
  psc_r_selection = psc_data(:,6);
  
  save([project_id '_selected_psc_sel' num2str(z) '.mat'],'psc_r_selection','psc_az_selection','grid_array_az','grid_array_r');
end


fid_res = fopen([project_id '_resfile.txt'],'a');
fprintf(fid_res,'\n*******************************************\n');
fprintf(fid_res,'PSC selection\n');
fprintf(fid_res,'*******************************************\n\n');

fprintf(1,'\n');
for z = 1:Npsc_selections
  fprintf(1,'In total, %4.0f PSC have been selected in selection %2.0f\n',Npsc(z),z);

  fprintf(fid_res,'Number of PSCs in selection %g: %g\n',z,Npsc(z));
end
fclose('all');
