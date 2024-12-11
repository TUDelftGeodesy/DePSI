function ps_spatial_unwrapping_plots

% Function to plot spatial unwrapping networks
%
% ----------------------------------------------------------------------
% File............: ps_spatial_unwrapping_plots.m
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
% - improvement of plots
% v1.7.2.8, Freek van Leijen
% - unwrap2_istep
%


global project_id results_id fig unwrap_istep unwrap2_istep

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
    for w = 1:unwrap2_istep

      %%%%%%%%%%%% network psc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fid = fopen([project_id '_network_psc_sel' num2str(z) '_' results_id ...
		   '_unwrap_istep' num2str(unwrap_istep) ...
		   '_unwrap2_istep' num2str(w) '.mat'],'r');
      if fid>0
	fclose(fid);
	load([project_id '_network_psc_sel' num2str(z) '_' results_id ...
	      '_unwrap_istep' num2str(unwrap_istep) ...
	      '_unwrap2_istep' num2str(w) '.mat'])
	
	psc_az = psc_az_orig;
	psc_r = psc_r_orig;
	
	switch orbit
	 case 'desc'
          psc_az = Nlines-psc_az+1;
          psc_r = Npixels-psc_r+1;
	end
	
	Narcs = size(dpsc_arcs,1);
	
	% ----------------------------------------------------------------------
	% Plot results
	% ----------------------------------------------------------------------
	
	fig = fig+1;
	figure(fig);hold on;
	set(gcf,'visible','off');
	set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');
	subimage(mrm,gray(256));
	
	plot(psc_r/delta_r,psc_az/delta_az,'ob','linewidth',2);
	
	for v = 1:Narcs
	  plot(psc_r(dpsc_arcs(v,:))/delta_r,psc_az(dpsc_arcs(v,:))/delta_az,'b')
	end
	axis off
	hold off
	
	print('-dpng',['plots/' project_id '_1_network_psc_sel' num2str(z) ...
		       '_' results_id '_unwrap_istep' num2str(unwrap_istep) ...
		       '_unwrap2_istep' num2str(w) '.png']);
	
	close(fig);
      end
      
      
      %%%%%%%%%%%% apriori removed arcs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fid = fopen([project_id '_apriori_removed_arcs_sel' num2str(z) '_' ...
		   results_id '_unwrap_istep' num2str(unwrap_istep) ...
		   '_unwrap2_istep' num2str(w) '.mat'],'r');
      if fid>0
	fclose(fid);
	load([project_id '_apriori_removed_arcs_sel' num2str(z) '_' ...
	      results_id '_unwrap_istep' num2str(unwrap_istep) ...
		   '_unwrap2_istep' num2str(w) '.mat'])
	
	switch orbit
	 case 'desc'
          psc_az = Nlines-psc_az+1;
          psc_r = Npixels-psc_r+1;
	end
	
	% ----------------------------------------------------------------------
	% Plot results
	% ----------------------------------------------------------------------
	
	fig = fig+1;
	figure(fig);hold on;
	set(gcf,'visible','off');
	set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');
	subimage(mrm,gray(256));
	
	plot(psc_r/delta_r,psc_az/delta_az,'ob','linewidth',2);
	
	Narcs_orig = size(dpsc_arcs_orig,1);
	for v = 1:Narcs_orig
	  plot(psc_r(dpsc_arcs_orig(v,:))/delta_r,...
	       psc_az(dpsc_arcs_orig(v,:))/delta_az,'b')
	end
	
	Narcs_removed1 = size(dpsc_arcs_removed1,1);
	if Narcs_removed1~=0
	  for v = 1:Narcs_removed1
	    h1 = plot(psc_r(dpsc_arcs_removed1(v,:))/delta_r,...
		      psc_az(dpsc_arcs_removed1(v,:))/delta_az,...
		      'color',[255 140 0]./255,'linewidth',2);
        end
	else
	  h1 = [];
	end
	
	Narcs_removed2 = size(dpsc_arcs_removed2,1);
	if Narcs_removed2~=0
	  for v = 1:Narcs_removed2
	    h2 = plot(psc_r(dpsc_arcs_removed2(v,:))/delta_r,...
		      psc_az(dpsc_arcs_removed2(v,:))/delta_az,...
		      'r','linewidth',2);
	  end
	else
	  h2 = [];
	end
	
	if ~isempty(h1)&~isempty(h2)
	  legend([h1 h2],[{'Apriori removed'},{'Untestable'}],...
		 'Location','SouthOutside','Orientation','horizontal');
	elseif ~isempty(h1)
	  legend([h1],[{'Apriori removed'}],...
		 'Location','SouthOutside','Orientation','horizontal');
	elseif ~isempty(h2)
	  legend([h2],[{'Untestable'}],...
		 'Location','SouthOutside','Orientation','horizontal');
	end
	
        axis off
	
	print('-dpng',['plots/' project_id '_2_apriori_removed_arcs_sel' ...
		       num2str(z) '_' results_id '_unwrap_istep' ...
		       num2str(unwrap_istep) ...
		       '_unwrap2_istep' num2str(w) '.png'])
	close(fig);
      end
      
      
      %%%%%%%%%%%% removed arcs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fid = fopen([project_id '_removed_arcs_sel' num2str(z) '_' results_id ...
		   '_unwrap_istep' num2str(unwrap_istep) ...
		   '_unwrap2_istep' num2str(w) '.mat'],'r');
      if fid>0
	fclose(fid);
	load([project_id '_removed_arcs_sel' num2str(z) '_' results_id ...
	      '_unwrap_istep' num2str(unwrap_istep) ...
	      '_unwrap2_istep' num2str(w) '.mat'])
	
	switch orbit
	 case 'desc'
          psc_az = Nlines-psc_az+1;
          psc_r = Npixels-psc_r+1;
          psc_az_new = Nlines-psc_az_new+1;
          psc_r_new = Npixels-psc_r_new+1;
          ref_array(:,2) = Nlines-ref_array(:,2)+1;
          ref_array(:,3) = Npixels-ref_array(:,3)+1;
	end
	
	% ----------------------------------------------------------------------
	% Plot results
	% ----------------------------------------------------------------------
	
	fig = fig+1;
	figure(fig);hold on;
	set(gcf,'visible','off');
	set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');
	subimage(mrm,gray(256));
	
	plot(psc_r_new/delta_r,psc_az_new/delta_az,'ob','linewidth',2);
	
	Narcs_new = size(dpsc_arcs,1);
	for v = 1:Narcs_new
	  plot(psc_r(dpsc_arcs(v,:))/delta_r,psc_az(dpsc_arcs(v,:))/delta_az,'b')
	end
	
	Nref = size(ref_array,1);
	for v = 1:Nref
	  h1 = plot(ref_array(v,3)/delta_r,ref_array(v,2)/delta_az,'g^', ...
		    'markersize',10,'linewidth',2,'MarkerFaceColor','none');
	end
	
	Nadapt = size(adapt_vec,1);
	if Nadapt~=0
	  for v = 1:Nadapt
	    h2 = plot(psc_r(dpsc_arcs(adapt_vec(v,1),:))/delta_r, ...
		      psc_az(dpsc_arcs(adapt_vec(v,1),:))/delta_az,'r', ...
		      'linewidth',2);
	  end
	else
	  h2 = [];
	end
        
	if ~isempty(h2)
	  legend([h1 h2],[{'Reference PS'},{'Removed'}],...
		 'Location','SouthOutside','Orientation','horizontal');
	else
	  legend([h1],[{'Reference PS'}],...
		 'Location','SouthOutside','Orientation','horizontal');
	end
	
	axis off
	
	print('-dpng',['plots/' project_id '_3_removed_arcs_sel' num2str(z) ...
		       '_' results_id '_unwrap_istep' num2str(unwrap_istep) ...
		       '_unwrap2_istep' num2str(w) '.png'])
	close(fig);
	
      end
      
      
      %%%%%%%%%%%% adapted arcs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      fid = fopen([project_id '_adapted_arcs_sel' num2str(z) '_' results_id ...
		   '_unwrap_istep' num2str(unwrap_istep) ...
		   '_unwrap2_istep' num2str(w) '.mat'],'r');
      if fid>0
	fclose(fid);
	load([project_id '_adapted_arcs_sel' num2str(z) '_' results_id ...
	      '_unwrap_istep' num2str(unwrap_istep) ...
		   '_unwrap2_istep' num2str(w) '.mat'])
	
	switch orbit
	 case 'desc'
          psc_az = Nlines-psc_az+1;
          psc_r = Npixels-psc_r+1;
          ref_array(:,2) = Nlines-ref_array(:,2)+1;
          ref_array(:,3) = Npixels-ref_array(:,3)+1;
	end
	
	psc = unique(dpsc_arcs_new);
	psc_az_new = psc_az(psc);
	psc_r_new = psc_r(psc);
	
	% ----------------------------------------------------------------------
	% Plot results
	% ----------------------------------------------------------------------
	
	fig = fig+1;
	figure(fig);hold on;
	set(gcf,'visible','off');
	set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');
	subimage(mrm,gray(256));
	
	plot(psc_r_new/delta_r,psc_az_new/delta_az,'ob','linewidth',2);
	
	Narcs_new = size(dpsc_arcs_new,1);
	for v = 1:Narcs_new
	  plot(psc_r(dpsc_arcs_new(v,:))/delta_r,psc_az(dpsc_arcs_new(v,:))/delta_az,'b')
	end
	
	Nref = size(ref_array,1);
	for v = 1:Nref
	  h1 = plot(ref_array(v,3)/delta_r,ref_array(v,2)/delta_az,'g^', ...
		    'markersize',10,'linewidth',2,'MarkerFaceColor','none');
	end
	
	if ~isempty(adapt_vec)
	  Nel = histc(adapt_vec(:,1),0.5:1:Narcs_new-0.5);
	  cmap = flipud(hot(max(Nel)+1));
	  for v = 1:Narcs_new
	    plot(psc_r(dpsc_arcs_new(v,:))/delta_r,...
		 psc_az(dpsc_arcs_new(v,:))/delta_az,'color',cmap(Nel(v)+1,:),...
		 'linewidth',2);
	  end
	  colormap(cmap);
	  k = colorbar;
	  set(get(k,'title'),'string','Nadapt');
	  set(k,'Ytick',[1:ceil(max(Nel)/5):max(Nel)+1]);
	  set(k,'Yticklabel',cellstr(num2str([1:ceil(max(Nel)/5):max(Nel)+1]'-1)));
	end
	
	legend([h1],[{'Reference PS'}],...
	       'Location','SouthOutside','Orientation','horizontal');

	axis off
	
	print('-dpng',['plots/' project_id '_4_adapted_arcs_sel' num2str(z) ...
		       '_' results_id '_unwrap_istep' num2str(unwrap_istep) ...
		       '_unwrap2_istep' num2str(w) '.png'])
	
	close(fig);
	
	%%%%%%%%%%%% final network %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	fig = fig+1;
	figure(fig);hold on;
	set(gcf,'visible','off');
	set(gcf,'position',figpos,'PaperPositionMode','auto','InvertHardcopy','off');
	subimage(mrm,gray(256));
	
	plot(psc_r_new/delta_r,psc_az_new/delta_az,'ob','linewidth',2);
	
	Narcs_new = size(dpsc_arcs_new,1);
	for v = 1:Narcs_new
	  plot(psc_r(dpsc_arcs_new(v,:))/delta_r,psc_az(dpsc_arcs_new(v,:))/delta_az,'b')
	end
	
	Nref = size(ref_array,1);
	for v = 1:Nref
	  h1 = plot(ref_array(v,3)/delta_r,ref_array(v,2)/delta_az,'g^', ...
		    'markersize',10,'linewidth',2,'MarkerFaceColor','none');
	end
	
	legend([h1],[{'Reference PS'}],...
	       'Location','SouthOutside','Orientation','horizontal');
	
	axis off
	
	print('-dpng',['plots/' project_id '_5_final_network_sel' num2str(z) ...
		       '_' results_id '_unwrap_istep' num2str(unwrap_istep) ...
		       '_unwrap2_istep' num2str(w) '.png'])
	
	close(fig);
      
      end
      
    end
  end 
end
fclose('all');
close all
