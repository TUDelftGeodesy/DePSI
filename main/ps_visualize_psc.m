function ps_visualize_psc(Nifgs,Npsc,Npsc_selections,final_model)

% Function to vizualize the estimated psc's
%
% Input: - Nifgs            number of interferograms
%        - Npsc             number of persistent scatterer candidates
%        - Npsc_selections  number of psc selections
%        - final_model      model used for unwrapped data
%
% ----------------------------------------------------------------------
% File............: ps_visualize_psc.m
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
% v1.7.2.8, Freek van Leijen
% - fixed color axis via ps_model_definitions
% 


% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global fig Nlines Npixels visible_plots Npar_max az_spacing r_spacing
global project_id orbit max_mem_buffer results_id

load cmap_topo
load cmap_defo


% ----------------------------------------------------------------------
% Read data
% ----------------------------------------------------------------------

if r_spacing>=az_spacing
  delta_r = max(1,round(sqrt(Nlines*Npixels)/1000));
  delta_az = max(1,delta_r*round(r_spacing/az_spacing));
else
  delta_az = max(1,round(sqrt(Nlines*Npixels)/1000));
  delta_r = max(1,delta_az*round(az_spacing/r_spacing));
end
Nlines_new = ceil(Nlines/delta_az);
Npixels_new = ceil(Npixels/delta_r);

mrm = zeros(Nlines_new,Npixels_new,'uint8');

Nlines_buffer = min(Nlines,floor(max_mem_buffer/(8*Npixels)));
Nbuffer = ceil(Nlines/Nlines_buffer);
Nlines_rem_buffer = rem(Nlines,Nlines_buffer);

mrm_fid = fopen([project_id '_mrm_uint8.raw'],'r');
sl = 1;
el = 0;
for v = 1:Nbuffer
  if (v==Nbuffer)&&(Nlines_rem_buffer~=0)
    Nlines_buffer = Nlines_rem_buffer;
  end
  temp = fread(mrm_fid,[Npixels Nlines_buffer],'*uint8')';
  temp = temp(1:delta_az:end,1:delta_r:end);
  Nlines_buffer_new = size(temp,1);
  el = el+Nlines_buffer_new;
  mrm(sl:el,:) = temp;
  sl = sl+Nlines_buffer_new;
end
fclose(mrm_fid);

switch orbit
 case 'desc'
  mrm = flipud(fliplr(mrm));
end

%fullscreen = [1 1 1280 1024];
fullscreen=get(0,'Screensize');
[minratio,ratio_index] = min([fullscreen(3)/Npixels_new fullscreen(4)/Nlines_new]);
if ratio_index(1)==1 %(1) to prevent multiple maxima
  figpos = [1 1 fullscreen(3) fullscreen(3)*Nlines_new/Npixels_new];
elseif ratio_index(1)==2
  figpos = [1 1 fullscreen(4)*Npixels_new/Nlines_new fullscreen(4)];
end

[Npar,par_index,covar_index] = ps_model_definitions('Npar',final_model);
[psc_annotation,par_index,covar_index] = ps_model_definitions('annotation',final_model);
[plot_caxis,par_index,covar_index] = ps_model_definitions('caxis',final_model);

% ----------------------------------------------------------------------
% Loop
% ----------------------------------------------------------------------

for z = 1:Npsc_selections
  
  
  % ----------------------------------------------------------------------
  % Read data
  % ----------------------------------------------------------------------
  
  psc_fid = fopen([project_id '_psc_sel' num2str(z) '.raw'],'r'); 
  psc_data = fread(psc_fid,[2*Nifgs+7 Npsc(z)],'double')';
  fclose(psc_fid);
  
  psc_array = psc_data(:,1:2);
  psc_az = psc_data(:,5);
  psc_r = psc_data(:,6);
  clear psc_data
  
  psc_results_fid = fopen([project_id '_psc_results_sel' num2str(z) ...
                      '.raw'],'r'); 
  psc_data = fread(psc_results_fid,[3*Nifgs+Npar_max+2 Npsc(z)],'double')';
  fclose(psc_results_fid);
  
  psc_param = psc_data(:,1:Npar_max);
  psc_ens_coh = psc_data(:,Npar_max+1);
  clear psc_data

  
  psc = find(psc_array(:,2)~=0); % remove isolated PSCs
  psc_az_new = psc_az(psc);
  psc_r_new = psc_r(psc);
  psc_param_new = psc_param(psc,:);
  psc_ens_coh_new = psc_ens_coh(psc);
  
  switch orbit
   case 'desc'
    psc_r_new = Npixels-psc_r_new+1;
    psc_az_new = Nlines-psc_az_new+1;
  end
  
  
  % ----------------------------------------------------------------------
  % plot parameters
  % ----------------------------------------------------------------------
  
  for v = 1:Npar
      %clim = [min(psc_param_new(:,par_index(v))) max(psc_param_new(:,par_index(v)))];
      clim = plot_caxis(par_index(v),:);
      
      if strcmp(char(psc_annotation(par_index(v))),'topo')
        cmap = cmap_topo;
      else
        cmap = cmap_defo;
      end
      
      fig = fig + 1;
      figure(fig);hold on;
      if strcmp(visible_plots,'n')
          set(gcf,'visible','off');
      end
      set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');

      subimage(mrm,gray(256));
      colormap(cmap);
      scatter(psc_r_new/delta_r,psc_az_new/delta_az,20,psc_param_new(:,par_index(v)),'filled');
      caxis(clim);
      title(char(psc_annotation(par_index(v))),'interpreter','none');
      xlabel('range');
      ylabel('azimuth');
      cbar_label = char(psc_annotation(par_index(v)+Npar_max));
      chandle = colorbar;
      set(get(chandle,'title'),'string',cbar_label);
      set(gca,'YDir','normal');
      axis off
      hold off

      print('-dpng',['plots/' project_id '_psc_parameter_' char(psc_annotation(par_index(v))) '_model' num2str(final_model) '_sel' num2str(z) '_' results_id '.png']);
      close(fig);

  end
  

  % ----------------------------------------------------------------------
  % plot ensemble coherence
  % ----------------------------------------------------------------------
  
  clim = [0 1];
  cmap = flipud(cmap_defo);
  
  fig = fig + 1;
  figure(fig);hold on;
  if strcmp(visible_plots,'n')
    set(gcf,'visible','off');
  end
  set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');
  
  subimage(mrm,gray(256));
  colormap(cmap);
  scatter(psc_r_new/delta_r,psc_az_new/delta_az,20,psc_ens_coh_new,'filled');
  caxis(clim);
  title('Ensemble coherence');
  xlabel('range');
  ylabel('azimuth');
  cbar_label = '[-]';
  chandle = colorbar;
  set(get(chandle,'title'),'string',cbar_label);
  set(gca,'YDir','normal');
  axis off
  hold off

  print('-dpng',['plots/' project_id '_psc_parameter_ens_coh_model' num2str(final_model) '_sel' num2str(z) '_' results_id '.png']);
  close(fig);

end

close all;
