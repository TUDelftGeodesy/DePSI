function ps_selection_plots

% Function to plot selected ps
%
% ----------------------------------------------------------------------
% File............: ps_selection_plots.m
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



global project_id

fid = fopen([project_id '_project.mat'],'r');
if fid>0
  fclose(fid);
  load([project_id '_project.mat'])

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
  
  for z = 1:Npsc_selections

    %%%%%%%%%%%% selected psc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fid = fopen([project_id '_selected_psc_sel' num2str(z) '.mat'],'r');
    if fid>0
      fclose(fid);
      load([project_id '_selected_psc_sel' num2str(z) '.mat']) 
      
      switch orbit
       case 'desc'
        psc_az_selection = Nlines-psc_az_selection+1;
        psc_r_selection = Npixels-psc_r_selection+1;
        grid_array_az(:,1) = Nlines-grid_array_az(:,1)+1;
        grid_array_r(:,1) = Npixels-grid_array_r(:,1)+1;
      end

      
      % ----------------------------------------------------------------------
      % Plot results
      % ----------------------------------------------------------------------
      
      fig = fig+1;
      figure(fig);hold on;
      set(gcf,'visible','off');
      set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');
      subimage(mrm,gray(256));
      plot(psc_r_selection/delta_r,psc_az_selection/delta_az,'or','linewidth',2);
      for v = 2:size(grid_array_az,1)
        plot([1 Npixels_new],repmat(grid_array_az(v,1)/delta_az,1,2));
      end
      for v = 2:size(grid_array_r,1)
        plot(repmat(grid_array_r(v,1)/delta_r,1,2),[1 Nlines_new]);
      end
      axis off
      hold off

      print('-dpng',['plots/' project_id '_selected_psc_sel' num2str(z) ...
                     '.png'])
      
    end
  end
end
