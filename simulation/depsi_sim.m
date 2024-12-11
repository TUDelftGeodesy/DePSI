function depsi_sim(param_file,Nifgs,perc,sim_model,sim_method,std_atmo_sim,az_spacing_sim)

%e.g,
%param_file = 'param_file_v1.7.0unwrap_sim1.txt';
%Nifgs = 30;
%perc = 0.001;
%sim_model = 5;
%sim_method = 'shape';
%std_atmo_sim = 0.005; %[m]
%az_spacing_sim = 2;
%Nlines = 256

% depsi_sim
% 
% Usage: depsi_sim('filename',Nifgs,perc,sim_model,sim_method,std_atmo_sim,az_spacing_sim)
%
% This is a simulation script for a PS analysis based on a single
% master stack of interferograms. The ambiguity function,
% bootstrapping as well as integer least-squares can
% be used for the estimations. 
%
% The program consists of a number of modules, which enables
% easy addition of steps and the formulation of alternatives.
%
% The modules are:
% - ps_selection:            reads the data in buffers and
%                            selects the psc, based on amplitude
%                            dispersion or high amplitude in a
%                            minimal number of scenes.
% - form_network:            forms a network between the psc
%                            by Delauney triangulation or a more
%                            redundant network.
% - psc_est:                 estimation of the ambiguities and
%                            parameters of psc by ambiguity
%                            function, bootstrapping or ils.
% - spatial_unwrapping:      spatial unwrapping of the network.
% - filtering:               estimation of atmospheric phase
%                            screen (APS) and non-linear defor-
%                            mation by filtering.
% - atmo_kriging:            atmosphere estimation by Kriging.
% - ps_est:                  estimation of the ambiguities and
%                            parameters of ps by ambiguity
%                            function, bootstrapping or ils.
%
% The program needs an input file (see manual
% for more information). Besides, you might need to change a few
% parameters at the beginning of this script.
%
% Using 'save [filename]' and 'load [filename]' you can easily save
% and load the workspace after a certain module. This prevents the
% need to re-run the whole program for a certain step.
%
% ----------------------------------------------------------------------
% File............: depsi_sim.m
% Version & Date..: 1.7.2.8, 19-OCT-2009
% Authors.........: Freek van Leijen
%                   Gini Ketelaar
%                   Petar Marinkovic
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
% - filenames in cells
% v1.7.2.17, Freek van Leijen
% - added global weighted_unwrap
%



fprintf(1,'\n');
fprintf(1,'PS analysis has started....\n');

if exist('plots','dir')~=7
  mkdir('plots');
end
if exist('plots/atmospheric_filter','dir')~=7
  mkdir('plots/atmospheric_filter');
end
if exist('plots/noise_filter','dir')~=7
  mkdir('plots/noise_filter');
end

global fig ps_method ps_eval_method project_id results_id
global m2ph Npar_max unwrap_istep unwrap2_istep
global max_mem_buffer orbit Nlines Npixels az_spacing r_spacing
global Npsc_selections weighting varfac_threshold
global visible_plots detail_plots ens_coh_threshold
global atmo_range defo_range sat_vel weighted_unwrap
global detrend_flag defo_model_flag

if exist(param_file)
  project_id = readinput('project_id',[],param_file);
else
  error('The parameter file you specified does not exist.');
end

processing_groups = readinput('processing_groups',[],param_file);
if isempty(processing_groups)
  processing_groups = 1:8;
end


if find(processing_groups == 1)

  if nargin==1
    Nifgs = 30;
    perc = 0.1;
    sim_model = 1;
    sim_method = 'shape';
    std_atmo_sim = 0.003; %[m]
  end
  
  ps_readinput_parameters; % script to read input parameters
  if isempty(processing_groups)
    processing_groups = 1:8;
  end

  % ----------------------------------------------------------------------
  % 1
  % Initialization
  % ----------------------------------------------------------------------

  
  % Counters
  % ----------------------------------------------------------------------
  
  fig = 0; %counter of figures
  detrend_flag = 'no';
  defo_model_flag = 'no';
  Nref = NaN(Npsc_selections,1);
  
  
  % General parameters
  % ----------------------------------------------------------------------
  
  if strcmp(sensor,'ers')
    lambda = 0.0565646;    % wavelength [m]
    sat_vel = 7550;        % satellite velocity [m/s]
  elseif strcmp(sensor,'asar')
    lambda = 0.0562356;    % wavelength [m]
    sat_vel = 7550;        % satellite velocity [m/s]
  elseif strcmp(sensor,'rsat')
    lambda = 0.0564105;    % wavelength [m] carrier_f=5300432000 Hz, c=2.99e8
    sat_vel = 7550;        % satellite velocity [m/s]
  elseif strcmp(sensor,'tsx')
    lambda = 0.031000;     % wavelength [m]
    sat_vel = 7600; %???   % satellite velocity [m/s]
  else
    error('The sensor you specified is not supported yet');
  end
  
  m2ph = -4*pi/lambda; % meters to phase constant
  toolbox_version = '1.7.2.8, 19-OCT-2009';
  


  % Model parameters
  % ----------------------------------------------------------------------
  
  Npar_max = ps_model_definitions('Npar_max');
  
  if length(final_model)>1
    warning('You specified multiple final models. Using the first.');
    final_model = final_model(1);
  end
  
  if isempty(std_param)
    std_param = [40 0.1 2 0.02 0.01 0.01 0.02];
    % topo, master_atmo, subpixel, linear, quadratic, cubic, periodic
  end
  
  if ~isempty(breakpoint)
    if breakpoint>1000
      breakpoint = find(str2num(orbitnr)==breakpoint);
    end
  end

  
  fid_res = fopen([project_id '_resfile.txt'],'w');
  fprintf(fid_res,'\n*******************************************\n\n');
  fprintf(fid_res,['Persistent Scatterer simulation ' project_id ' \n']);
  fprintf(fid_res,['Version: ' toolbox_version '\n']);
  fprintf(fid_res,'Result file\n');
  fprintf(fid_res,'%s\n\n',date);
  fprintf(fid_res,'*******************************************\n\n');
  fprintf(fid_res,'Input file:               %s\n',input_file);
  fprintf(fid_res,'Number of interferograms: %3.0f\n',Nifgs);
  fclose(fid_res);


  
  % ----------------------------------------------------------------------
  % 1b
  % PS input simulation
  % ----------------------------------------------------------------------
  
  [topo,...
   defo1,...
   defo2,...
   defo3,...
   defo_quad,...
   defo_cubic,...
   defo_nl1,...
   defo_nl2,...
   subpixel,...
   atmo,...
   ifgs,...
   ifgs_unw,...
   ambi,...
   breakpoint,...
   breakpoint2,...
   noise_level,...
   Btemp,...
   Bperp,...
   Bdop] = ifgs_simulation(Nlines,...
                           Nifgs,...
                           az_spacing_sim,...
                           sat_vel,...
                           sim_model,...
                           sim_method,...
                           std_param(1),...
                           std_param(4),...
                           std_atmo_sim);

  theta = 23*pi/180; %incident angle [deg]
  H = 780000; %satellite vertical height [m]
  R = H/cos(theta); %antenna-target distance [m]
  h2ph = Bperp'/(R*sin(theta)); %height to phase factors
  
  Npixels = Nlines;
  crop = [1 Nlines 1 Npixels];
  crop_in = crop;
  
  datenums = datenum(date)+(rand(Nifgs+1,1)-1)*10*365;
  dates = datestr([sort(datenums(1:Nifgs));datenums(Nifgs+1)],'dd-mmm-yyyy');
  
  filenames_ifgs = cell(Nifgs,1);
  filenames_ifgs_unw = cell(Nifgs,1);
  filenames_h2ph = cell(Nifgs,1);
  filenames_output = cell(Nifgs,1);

  for v = 1:Nifgs
    filenames_ifgs(v) = {[project_id '_ifgs' num2str(v,'%02d') '.raw']};
    filenames_ifgs_unw(v) = {[project_id '_ifgs_unw' num2str(v,'%02d') '.raw']};
    filenames_h2ph(v) = {[project_id '_h2ph' num2str(v,'%02d') '.raw']};
    filenames_output(v) = {[project_id num2str(v,'%02d')]};
    
    ifgs_fid = fopen(char(filenames_ifgs(v)),'w');
    fwrite(ifgs_fid,ifgs(:,:,v)','single');
    fclose(ifgs_fid);
    
    h2ph_fid = fopen(char(filenames_h2ph(v)),'w');
    fwrite(h2ph_fid,repmat(h2ph(v),Npixels,Nlines),'single');
    fclose(h2ph_fid);
    
    ifgs_fid = fopen(char(filenames_ifgs_unw(v)),'w');
    fwrite(ifgs_fid,ifgs_unw(:,:,v)','single');
    fclose(ifgs_fid);
    
  end
  
  save([project_id '_project.mat']);
  clear [^project_id ^processing_groups];

end


if find(processing_groups == 2)

  load([project_id '_project.mat']);
  ps_readinput_parameters; % script to read input parameters
  if isempty(processing_groups)
    processing_groups = 1:8;
  end
  crop = [1 Nlines 1 Npixels];

  % ----------------------------------------------------------------------
  % 2
  % PSC selection
  % ----------------------------------------------------------------------
  
  fprintf(1,'\n');
  fprintf(1,'Initial PSC selection....\n');

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,'\n');
  fprintf(fid_res,[datestr(now) ', group2, start selection ...\n']);
  fclose(fid_res);

  [grid_array_az,...
   grid_array_r,...
   Npsc,...
   Npsp] = ps_selection_sim(Nifgs,...
                            psc_selection_gridsize,...
                            Npsc_selections,...
                            psc_selection_method,...
                            ifgs,...
                            ambi,...
                            topo,...
                            defo1,...
                            defo2,...
			    defo3,...
			    defo_quad,...
			    defo_cubic,...
                            defo_nl1,...
                            defo_nl2,...
                            subpixel,...
                            atmo,...
                            Bperp,...
                            perc);
  
  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group2, end selection.\n']);
  fclose(fid_res);

  save([project_id '_project.mat']);
  
  ps_selection_plots;
  
  clear [^project_id ^processing_groups];
  
end


if find(processing_groups == 3)
  
  load([project_id '_project.mat']);
  ps_readinput_parameters; % script to read input parameters
  if isempty(processing_groups)
    processing_groups = 1:8;
  end
  crop = [1 Nlines 1 Npixels];

  for z = 1:Npsc_selections
    copyfile([project_id '_psc_2orig_sel' num2str(z) '.raw'],[project_id '_psc_sel' num2str(z) '.raw']);
    
    switch ps_eval_method
      case 'psp'
        copyfile([project_id '_psp_2orig_sel' num2str(z) '.raw'],[project_id '_psp_sel' num2str(z) '.raw']);
    end
  end
  
  switch detrend_method
    case 'yes'
      
      fprintf(1,'\n');
      fprintf(1,'Trend removal....\n');
      
      % ----------------------------------------------------------------------
      % 3a
      % PS network
      % ----------------------------------------------------------------------
      
      fprintf(1,'\n');
      fprintf(1,'Network construction....\n');
    
      fid_res = fopen([project_id '_resfile.txt'],'a');
      fprintf(fid_res,'\n');
      fprintf(fid_res,[datestr(now) ', group3, start network construction ...\n']);
      fclose(fid_res);

      results_id = 'orig';
      
      [sig2_est,...
       Nref,...
       final_althyp_index] = ps_network_sim(Nifgs,...
                                            Npsc,...
                                            Npsc_selections,...
                                            Nref,...
                                            max_arc_length,...
                                            psc_model,...
                                            final_model,...
                                            Btemp,...
                                            Bdop,...
                                            std_param,...
                                            breakpoint,...
                                            breakpoint2,...
                                            network_method,...
                                            Ncon,...
                                            Nparts,...
                                            ens_coh_threshold,...
                                            varfac_threshold,...
                                            ref_cn);
      
      fid_res = fopen([project_id '_resfile.txt'],'a');
      fprintf(fid_res,[datestr(now) ', group3, end network construction.\n']);
      fclose(fid_res);

      save([project_id '_project.mat']);
      
      ps_spatial_unwrapping_plots; % make plots of networks

      ps_visualize_psc(Nifgs,Npsc,Npsc_selections,final_model); % make plots of psc
      
      % ----------------------------------------------------------------------
      % 4
      % PS detrend (based on psc)
      % ----------------------------------------------------------------------
      
      fprintf(1,'\n');
      fprintf(1,'Detrend (based on psc) ....\n');

      fid_res = fopen([project_id '_resfile.txt'],'a');
      fprintf(fid_res,[datestr(now) ', group3, start detrend ...\n']);
      fclose(fid_res);
      
      ps_detrend(Nifgs,...
                 Npsc,...
                 Npsp,...
                 Npsc_selections,...
                 filenames_output);

      fid_res = fopen([project_id '_resfile.txt'],'a');
      fprintf(fid_res,[datestr(now) ', group3, end detrend.\n']);
      fclose(fid_res);

      detrend_flag = 'yes';
      
    otherwise
      
      fprintf(1,'\n');
      fprintf(1,'No trend is removed....\n');
      
      detrend_flag = 'no';
      for z = 1:Npsc_selections
        copyfile([project_id '_psc_2orig_sel' num2str(z) '.raw'],[project_id '_psc_3detr_sel' num2str(z) '.raw']);
        
        switch ps_eval_method
          case 'psp'
            copyfile([project_id '_psp_2orig_sel' num2str(z) '.raw'],[project_id '_psp_3detr_sel' num2str(z) '.raw']);
        end
      end

  end

  save([project_id '_project.mat']);
  clear [^project_id ^processing_groups];

end


if find(processing_groups == 4)

  load([project_id '_project.mat']);
  ps_readinput_parameters; % script to read input parameters
  if isempty(processing_groups)
    processing_groups = 1:8;
  end
  crop = [1 Nlines 1 Npixels];

  for z = 1:Npsc_selections
    copyfile([project_id '_psc_3detr_sel' num2str(z) '.raw'],[project_id '_psc_sel' num2str(z) '.raw']);
    
    switch ps_eval_method
      case 'psp'
        copyfile([project_id '_psp_3detr_sel' num2str(z) '.raw'],[project_id '_psp_sel' num2str(z) '.raw']);
    end
  end

  % ----------------------------------------------------------------------
  % 3b
  % PS network
  % ----------------------------------------------------------------------
  
  fprintf(1,'\n');
  fprintf(1,'Network construction....\n');

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,'\n');
  fprintf(fid_res,[datestr(now) ', group4, start network construction ...\n']);
  fclose(fid_res);

  if strcmp(detrend_flag,'yes')
    results_id = 'detr';
  elseif strcmp(detrend_flag,'no')
    results_id = 'orig';
  end
  
  [sig2_est,...
   Nref,...
   final_althyp_index] = ps_network_sim(Nifgs,...
					Npsc,...
					Npsc_selections,...
					Nref,...
					max_arc_length,...
					psc_model,...
					final_model,...
					Btemp,...
					Bdop,...
					std_param,...
					breakpoint,...
					breakpoint2,...
					network_method,...
					Ncon,...
					Nparts,...
					ens_coh_threshold,...
					varfac_threshold,...
					ref_cn);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group4, end network construction.\n']);
  fclose(fid_res);
  
  save([project_id '_project.mat']);
  
  ps_spatial_unwrapping_plots; % make plots of networks
  ps_visualize_psc(Nifgs,Npsc,Npsc_selections,final_model); % make plots of psc    
  
  
  % ----------------------------------------------------------------------
  % 5a
  % PS filtering
  % ----------------------------------------------------------------------
  
  fprintf(1,'\n');
  fprintf(1,'Separation atmosphere and unmodeled deformation ....\n');

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group4, start filtering ...\n']);
  fclose(fid_res);
  
  ps_filtering(Btemp,...
	       ts_atmo_filter,...
	       ts_atmo_filter_length,...
	       Nifgs,...
	       Npsc,...
	       Npsc_selections);
  
  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group4, end filtering.\n']);
  fclose(fid_res);
  
  
  
  % ----------------------------------------------------------------------
  % 5b
  % Atmosphere estimation by Kriging
  % ----------------------------------------------------------------------
  
  fprintf(1,'\n');
  fprintf(1,'Atmosphere estimation by Kriging....\n');
  
  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group4, start atmosphere estimation ...\n']);
  fclose(fid_res);

  ps_atmo_kriging(Nifgs,...
		  Npsc,...
		  Npsc_selections,...
		  Npsp,...
		  filenames_output);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group4, end atmosphere estimation.\n']);
  fclose(fid_res);
  
  ps_visualize_atmo(Nifgs,Npsp,Npsc,Npsc_selections,dates);
  
  save([project_id '_project.mat']);
  clear [^project_id ^processing_groups];
  
end


if find(processing_groups == 5)

  load([project_id '_project.mat']);
  ps_readinput_parameters; % script to read input parameters
  if isempty(processing_groups)
    processing_groups = 1:8;
  end
  crop = [1 Nlines 1 Npixels];

  for z = 1:Npsc_selections
    copyfile([project_id '_psc_4atmo_sel' num2str(z) '.raw'],[project_id '_psc_sel' num2str(z) '.raw']);
    
    switch ps_eval_method
      case 'psp'
        copyfile([project_id '_psp_4atmo_sel' num2str(z) '.raw'],[project_id '_psp_sel' num2str(z) '.raw']);
    end
  end

  % ----------------------------------------------------------------------
  % 6
  % Make subset of available ifgs (optional)
  % ----------------------------------------------------------------------

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,'\n');
  fclose(fid_res);

  if isempty(ifg_selection_input)
    % do nothing
  else

    fid_res = fopen([project_id '_resfile.txt'],'a');
    fprintf(fid_res,[datestr(now) ', group5, start interferogram selection ...\n']);
    fclose(fid_res);

    if exist('Nifgs_orig')
      Nifgs = Nifgs_orig;
      Btemp = Btemp_orig;
      Bdop = Bdop_orig;
      Bperp = Bperp_orig;
    else
      Nifgs_orig = Nifgs;
      Btemp_orig = Btemp;
      Bdop_orig = Bdop;
      Bperp_orig = Bperp;
    end
    
    if ischar(ifg_selection_input) % 
      ifg_selection = ps_determine_ifg_selection(orbitnr(1:Nifgs,:),ifg_selection_input);
    else
      ifg_selection = ifg_selection_input;
    end
    Nifgs = ps_select_ifgs(ifg_selection,Npsc,Npsc_selections,Npsp,Nifgs);
    Btemp = Btemp(ifg_selection);
    Bdop = Bdop(ifg_selection);
    Bperp = Bperp(ifg_selection);

    fid_res = fopen([project_id '_resfile.txt'],'a');
    fprintf(fid_res,[datestr(now) ', group5, end interferogram selection.\n']);
    fclose(fid_res);
    
  end


  
  % ----------------------------------------------------------------------
  % 3c
  % PS network
  % ----------------------------------------------------------------------
  
  fprintf(1,'\n');
  fprintf(1,'Network construction....\n');

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group5, start network construction ...\n']);
  fclose(fid_res);

  results_id = 'atmo';

  [sig2_est,...
   Nref,...
   final_althyp_index] = ps_network_sim(Nifgs,...
                                        Npsc,...
                                        Npsc_selections,...
                                        Nref,...
                                        max_arc_length,...
                                        psc_model,...
                                        final_model,...
                                        Btemp,...
                                        Bdop,...
                                        std_param,...
                                        breakpoint,...
					breakpoint2,...
                                        network_method,...
                                        Ncon,...
                                        Nparts,...
                                        ens_coh_threshold,...
					varfac_threshold,...
                                        ref_cn);
  

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group5, end network construction.\n']);
  fclose(fid_res);

  save([project_id '_project.mat']);

  ps_spatial_unwrapping_plots; % make plots of networks
  ps_visualize_psc(Nifgs,Npsc,Npsc_selections,final_model); % make plots of psc
  

  % ----------------------------------------------------------------------
  % 7a
  % PS densification after atmosphere correction
  % ----------------------------------------------------------------------
  
  fprintf(1,'\n');
  fprintf(1,['PS densification by ' ps_method ' ....\n']);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group5, start densification ...\n']);
  fclose(fid_res);
  
  Nps_atmo = ps_densification_sim(Btemp,...
                                  Bdop,...
                                  grid_array_az,...
                                  grid_array_r,...
                                  Nifgs,...
                                  Npsc,...
                                  Npsc_selections,...
                                  Npsp,...
                                  Nref,...
                                  ps_model,...
                                  final_model,...
                                  final_althyp_index,...
                                  std_param,...
                                  sig2_est,...
                                  breakpoint,...
				  breakpoint2,...
                                  Namp_disp_bins,...
                                  Ndens_iterations,...
                                  densification_flag,...
                                  filenames_output,...
                                  filenames_ifgs,...
                                  filenames_h2ph,...
                                  filenames_ifgs_unw,...
                                  detrend_flag,...
                                  defo_model_flag,...
                                  ps_area_of_interest,...
                                  dens_method,...
                                  dens_check,...
                                  filename_water_mask,...
				  crop,...
				  crop_in,...
                                  Nest);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group5, end densification.\n']);
  fclose(fid_res);
  
  save([project_id '_project.mat']);
  
  ps_visualize_ps(Nifgs,Nps_atmo,Npsc,Npsc_selections,final_model); % make plots of ps

  save([project_id '_project.mat']);
  clear [^project_id ^processing_groups];

end


if find(processing_groups == 6)

  load([project_id '_project.mat']);
  ps_readinput_parameters; % script to read input parameters
  if isempty(processing_groups)
    processing_groups = 1:8;
  end
  crop = [1 Nlines 1 Npixels];

  if ~isempty(defo_method)
    
    % -------------------------------------------------------------------
    % 8
    % Estimation of the deformation model
    % -------------------------------------------------------------------
    
    fprintf(1,'\n');
    fprintf(1,['Estimation of deformation model....\n']);

    fid_res = fopen([project_id '_resfile.txt'],'a');
    fprintf(fid_res,'\n');
    fprintf(fid_res,[datestr(now) ', group6, start deformation modeling ...\n']);
    fclose(fid_res);
    
    ps_defo_model(Btemp,...
                  Nifgs,...
                  Npsc,...
                  Npsp,...
                  Nps_atmo,...
                  Npsc_selections,...
                  std_param,...
                  xc0,...
                  yc0,...
                  zc0,...
                  r0,...
                  r10,...
                  epoch,...
                  defo_method,...
                  filenames_output);

    fid_res = fopen([project_id '_resfile.txt'],'a');
    fprintf(fid_res,[datestr(now) ', group6, end deformation modeling.\n']);
    fclose(fid_res);
     
    defo_model_flag = 'yes';
    
    for z = 1:Npsc_selections
      copyfile([project_id '_psc_5defo_sel' num2str(z) '.raw'],[project_id '_psc_sel' num2str(z) '.raw']);
      
      switch ps_eval_method
        case 'psp'
          copyfile([project_id '_psp_5defo_sel' num2str(z) '.raw'],[project_id '_psp_sel' num2str(z) '.raw']);
      end
    end
    
    % ----------------------------------------------------------------------
    % 3d
    % PS network
    % ----------------------------------------------------------------------
    
    fprintf(1,'\n');
    fprintf(1,'Network construction....\n');

    fid_res = fopen([project_id '_resfile.txt'],'a');
    fprintf(fid_res,[datestr(now) ', group6, start network construction ...\n']);
    fclose(fid_res);
    
    results_id = 'defo';

    [sig2_est,...
     Nref,...
     final_althyp_index] = ps_network_sim(Nifgs,...
                                          Npsc,...
                                          Npsc_selections,...
                                          Nref,...
                                          max_arc_length,...
                                          psc_model,...
                                          final_model,...
                                          Btemp,...
                                          Bdop,...
                                          std_param,...
                                          breakpoint,...
                                          breakpoint2,...
                                          network_method,...
                                          Ncon,...
                                          Nparts,...
                                          ens_coh_threshold,...
                                          varfac_threshold,...
                                          ref_cn);

    fid_res = fopen([project_id '_resfile.txt'],'a');
    fprintf(fid_res,[datestr(now) ', group6, end network construction.\n']);
    fclose(fid_res);

    save([project_id '_project.mat']);

    ps_spatial_unwrapping_plots; % make plots of networks
    ps_visualize_psc(Nifgs,Npsc,Npsc_selections,final_model); % make plots of psc
  
  
    % -----------------------------------------------------------------
    % 7b
    % PS densification after defo correction
    % -----------------------------------------------------------------
    
    fprintf(1,'\n');
    fprintf(1,['PS densification by ' ps_method ' ....\n']);

    fid_res = fopen([project_id '_resfile.txt'],'a');
    fprintf(fid_res,[datestr(now) ', group6, start densification ...\n']);
    fclose(fid_res);
    
    Nps_defo = ps_densification_sim(Btemp,...
                                    Bdop,...
                                    grid_array_az,...
                                    grid_array_r,...
                                    Nifgs,...
                                    Npsc,...
                                    Npsc_selections,...
                                    Npsp,...
                                    Nref,...
                                    ps_model,...
                                    final_model,...
                                    final_althyp_index,...
                                    std_param,...
                                    sig2_est,...
                                    breakpoint,...
                                    breakpoint2,...
                                    Namp_disp_bins,...
                                    Ndens_iterations,...
                                    densification_flag,...
                                    filenames_output,...
                                    filenames_ifgs,...
                                    filenames_h2ph,...
                                    filenames_ifgs_unw,...
                                    detrend_flag,...
                                    defo_model_flag,...
                                    ps_area_of_interest,...
                                    dens_method,...
                                    dens_check,...
                                    filename_water_mask,...
                                    crop,...
                                    crop_in,...
                                    Nest);

    fid_res = fopen([project_id '_resfile.txt'],'a');
    fprintf(fid_res,[datestr(now) ', group6, end densification.\n']);
    fclose(fid_res);
    
    save([project_id '_project.mat'])
    
    ps_visualize_ps(Nifgs,Nps_defo,Npsc,Npsc_selections,final_model);

  else
    
    defo_model_flag = 'no';
    Nps_defo = [];

  end

  save([project_id '_project.mat']);
  clear [^project_id ^processing_groups];

end


if find(processing_groups == 7)

  load([project_id '_project.mat']);
  ps_readinput_parameters; % script to read input parameters
  if isempty(processing_groups)
    processing_groups = 1:8;
  end
  crop = [1 Nlines 1 Npixels];

  % ----------------------------------------------------------------------
  % 9
  % Output generation
  % ----------------------------------------------------------------------
  
  fprintf(1,'\n');
  fprintf(1,['Calculate spatio-temporal consistency ....\n']);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,'\n');
  fprintf(fid_res,[datestr(now) ', group7, start spatio-temporal consistency ...\n']);
  fclose(fid_res);
  
  ps_spatio_temporal_consistency(Nifgs,...
				 Nps_atmo,...
				 Nps_defo,...
				 Npsc_selections,...
				 defo_model_flag,...
				 stc_min_max,...
				 'orig',...
				 ts_noise_filter,...
				 ts_noise_filter_length);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group7, end spatio-temporal consistency.\n']);
  fclose(fid_res);
  
  %fprintf(1,'\n');
  %fprintf(1,['Create output ....\n']);
  
  %ps_output_sim(Nifgs,Nps_atmo,Nps_defo,Npsc_selections,final_model,defo_model_flag);
  
  %switch ps_eval_method
  %  case 'whole'  
  %   switch densification_flag
  %     case 'yes'
  %       [Nps_atmo_nodens,Nps_defo_nodens] = ps_output_nodens_sim(Nifgs,Nps_atmo,Nps_defo,Npsc_selections,final_model,defo_model_flag,master_res,linelo,pixlo,ref_height,output_format,dates,Btemp);
  %    
  %       ps_visualize_ps_nodens(Nifgs,Nps_atmo_nodens,Npsc,Npsc_selections,final_model);
  %       if ~isempty(defo_method)
  %         ps_visualize_ps_nodens(Nifgs,Nps_defo_nodens,Npsc,Npsc_selections,final_model);
  %       end
  %   end
  %end

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group7, start check results ...\n']);
  fclose(fid_res);

  ps_check_results_sim(Nifgs,Nps_atmo,Nps_defo,Npsc,Npsp,Nref,Npsc_selections,final_model,defo_model_flag);
  
  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group7, end check results.\n']);
  fclose(fid_res);

  fclose('all');

  save([project_id '_project.mat']);
  clear [^project_id ^processing_groups];

end

if find(processing_groups == 8)
  
  load([project_id '_project.mat']);
  ps_readinput_parameters; % script to read input parameters
  if isempty(processing_groups)
    processing_groups = 1:8;
  end
  crop = [1 Nlines 1 Npixels];

  % ----------------------------------------------------------------------
  % 10
  % Time series filtering
  % ----------------------------------------------------------------------
  
  fprintf(1,'\n');
  fprintf(1,['Filtering time series ....\n']);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,'\n');
  fprintf(fid_res,[datestr(now) ', group8, start noise filtering ...\n']);
  fclose(fid_res);

  ps_noise_filter(Nifgs,...
		  Nps_atmo,...
		  Nps_defo,...
		  Npsc_selections,...
		  Btemp,...
		  Bdop,...
		  defo_model_flag,...
		  final_model,...
		  final_althyp_index,...
		  std_param,...
		  sig2_est,...
		  ts_noise_filter,...
		  ts_noise_filter_length);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group8, end noise filtering.\n']);
  fclose(fid_res);
 
  fprintf(1,'\n');
  fprintf(1,['Calculate spatio-temporal consistency filtered....\n']);
  
  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group8, start spatio-temporal consistency ...\n']);
  fclose(fid_res);

  ps_spatio_temporal_consistency(Nifgs,...
				 Nps_atmo,...
				 Nps_defo,...
				 Npsc_selections,...
				 defo_model_flag,...
				 stc_min_max,...
				 'filt',...
				 ts_noise_filter,...
				 ts_noise_filter_length);

  fid_res = fopen([project_id '_resfile.txt'],'a');
  fprintf(fid_res,[datestr(now) ', group8, end spatio-temporal consistency.\n']);
  fclose(fid_res);

  fclose('all');

  save([project_id '_project.mat']);
  clear all

end

fprintf(1,'\n');
fprintf(1,'Done!\n');


% ----------------------------------------------------------------------
% The end
% ----------------------------------------------------------------------



