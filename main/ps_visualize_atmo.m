function ps_visualize_atmo(Nifgs,Npsp,Npsc,Npsc_selections,dates)

% Function to vizualize the estimated psc's
%
% Input: - Nifgs            number of interferograms
%        - Npsp             number of persistent scatterers
%        - Npsc             number of persistent scatterer candidates
%        - Npsc_selections  number of psc selections
%        - final_model      model used for unwrapped data
%
% ----------------------------------------------------------------------
% File............: ps_visualize_atmo.m
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
% v1.7.2.12, Freek van Leijen
% - bug fixed in filename
% 


% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global fig Nlines Npixels visible_plots Npar_max r_spacing az_spacing
global project_id orbit max_mem_buffer results_id m2ph

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

  psc = find(psc_array(:,2)~=0); % remove isolated PSCs
  psc_az_new = psc_az(psc);
  psc_r_new = psc_r(psc);
  
  %atmo
  psc_results_fid = fopen([project_id '_psc_results_sel' num2str(z) '.raw'],'r'); 
  psc_data = fread(psc_results_fid,[3*Nifgs+Npar_max+2 Npsc(z)],'double')';
  fclose(psc_results_fid);
  
  psc_atmo_master = m2ph*psc_data(psc,2); %second column of psc_param, dangerous...
  clear psc_data
  
  atmo_fid = fopen([project_id '_psc_atmo_filt_sel' num2str(z) '.raw'],'r'); 
  psc_atmo_estimates_orig = fread(atmo_fid,[Nifgs Npsc(z)],'double')';
  fclose(atmo_fid);
  psc_atmo_estimates_orig = psc_atmo_estimates_orig(psc,:);
  psc_atmo_estimates_orig = [psc_atmo_estimates_orig psc_atmo_master psc_atmo_master]; %master atmo

  
  
  delta2 = max(1,round(Npsp(z)/10000));
  Npsp_new = ceil(Npsp(z)/delta2);

  psp_data = NaN(Npsp_new,2*Nifgs+4);
  psp_atmo = NaN(Npsp_new,Nifgs);
  psp_atmo_master1 = NaN(Npsp_new,1);% master atmo (direct)
  psp_atmo_master2 = NaN(Npsp_new,1);% master atmo (mean of ifgs)
  
  Npsp_buffer = min(Npsp(z),floor(max_mem_buffer/(8*(Npar_max+8*Nifgs+6)*delta2))*delta2);
  Npsp_buffer_new = ceil(Npsp_buffer/delta2);
  Nbuffer = ceil(Npsp(z)/Npsp_buffer);
  Npsp_rem_buffer = rem(Npsp(z),Npsp_buffer);
  
  psp_fid = fopen([project_id '_psp_sel' num2str(z) '.raw'],'r'); 
  psp_atmo_fid = fopen([project_id '_psp_atmo_sel' num2str(z) '.raw'],'r');
  psp_atmo_master1_fid = fopen([project_id '_psp_atmo_master1_sel' num2str(z) '.raw'],'r');
  psp_atmo_master2_fid = fopen([project_id '_psp_atmo_master2_sel' num2str(z) '.raw'],'r');
  for v = 1:Nbuffer
    sl = (v-1)*Npsp_buffer_new+1;
    if (v==Nbuffer)&&(Npsp_rem_buffer~=0)
      Npsp_buffer = Npsp_rem_buffer;
      el = Npsp_new;
    else
      el = v*Npsp_buffer_new;
    end
    temp = fread(psp_fid,[2*Nifgs+4 Npsp_buffer],'double')';
    psp_data(sl:el,:) = temp(1:delta2:end,:);
    temp = fread(psp_atmo_fid,[Nifgs Npsp_buffer],'double')';
    psp_atmo(sl:el,:) = temp(1:delta2:end,:);
    temp = fread(psp_atmo_master1_fid,[1 Npsp_buffer],'double')';
    psp_atmo_master1(sl:el) = temp(1:delta2:end,:);% master atmo (direct)
    temp = fread(psp_atmo_master2_fid,[1 Npsp_buffer],'double')';
    psp_atmo_master2(sl:el) = temp(1:delta2:end,:);% master atmo (mean of ifgs)
  end
  fclose(psp_fid);
  fclose(psp_atmo_fid);
  fclose(psp_atmo_master1_fid);
  fclose(psp_atmo_master2_fid);

  psp_az = psp_data(:,3);
  psp_r = psp_data(:,4);
  clear psp_data

  psp_atmo_tot = [psp_atmo psp_atmo_master1 psp_atmo_master2];
  
  %ps_az = [psc_az_new;psp_az];
  %ps_r = [psc_r_new;psp_r];
  %ps_atmo_tot = [psc_atmo_tot_new;psp_atmo_tot];
  
  switch orbit
   case 'desc'
    psc_r_new = Npixels-psc_r_new+1;
    psc_az_new = Nlines-psc_az_new+1;
    psp_r = Npixels-psp_r+1;
    psp_az = Nlines-psp_az+1;
  end
  
  
  % ----------------------------------------------------------------------
  % plot parameters
  % ----------------------------------------------------------------------
  
  clim = [-10 10];
  cmap = cmap_defo;
  
  for v = 1:Nifgs+2 %ifgs+2xmaster
    
    fig = fig + 1;
    figure(fig);hold on;
    if strcmp(visible_plots,'n')
      set(gcf,'visible','off');
    end
    set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');
    
    subimage(mrm,gray(256));
    colormap(cmap);
    scatter(psp_r/delta_r,psp_az/delta_az,5,1000*psp_atmo_tot(:,v)/m2ph,'filled');
    scatter(psc_r_new/delta_r,psc_az_new/delta_az,30,1000*psc_atmo_estimates_orig(:,v)/m2ph,'filled');
    caxis(clim);
    if v<=Nifgs
      title(['atmo ' dates(Nifgs+1,:) '--' dates(v,:)],'interpreter','none');
    elseif v==Nifgs+1
      title(['master atmo (direct) ' dates(Nifgs+1,:)],'interpreter','none');
    else
      title(['master atmo (mean) ' dates(Nifgs+1,:)],'interpreter','none');
    end
    xlabel('range');
    ylabel('azimuth');
    cbar_label = 'mm';
    chandle = colorbar;
    set(get(chandle,'title'),'string',cbar_label);
    set(gca,'YDir','normal');
    axis off
    hold off
    
    if v<=Nifgs
      print('-dpng',['plots/' project_id '_psp_atmo' num2str(v,'%3.0d') '_' ...
		     dates(Nifgs+1,:) '_' dates(v,:) '_sel' num2str(z) '.png']);
    elseif v==Nifgs+1
      print('-dpng',['plots/' project_id '_psp_master_atmo_direct_' ...
		     dates(Nifgs+1,:) '_sel' num2str(z) '.png']);
    else
      print('-dpng',['plots/' project_id '_psp_master_atmo_mean_' ...
		     dates(Nifgs+1,:) '_sel' num2str(z) '.png']);
    end
    
    close(fig);

  end
end

close all;
