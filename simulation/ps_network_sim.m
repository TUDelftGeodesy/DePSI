function [sig2_est,Nref,final_althyp_index] = ps_network_sim(Nifgs,Npsc,Npsc_selections,Nref,max_arc_length,psc_model,final_model,Btemp,Bdop,std_param,breakpoint,breakpoint2,network_method,Ncon,Nparts,ens_coh_threshold,varfac_threshold,ref_cn)

% Network estimation
%
% Input:    - Nifgs              number of interferograms
%           - Npsc               number of psc's
%           - Npsc_selections    number of psc_selections
%           - Nref               number of reference points
%           - max_arc_length     maximum arc length
%           - psc_model          vector with model indices to be
%                                evaluated (see help model_definitions
%                                for more information)
%           - final_model        model used for unwrapped data
%           - Btemp              temporal baselines [year]
%           - Bdop               Doppler baselines [Hz]
%           - std_param          standard deviations for pseudo-
%                                observations
%           - breakpoint         breakpoint in case of double
%                                linear model
%           - breakpoint2        second breakpoint, should be larger
%                                than first breakpoint 
%           - network_method     network identifier, 'delaunay' or
%                                'spider'
%           - Ncon               minimum number of connections (arcs)
%                                to a psc (network 2 only)
%           - Nparts             number of partitions of a full cycle 
%                                to which the arcs are divided
%                                (network 2 only)
%           - ens_coh_threshold  ensemble coherence threshold
%           - varfac_threshold   variance factor threshold
%           - ref_cn             coordinates reference point
%  
% Output:   - sig2_est           variance components
%           - Nref               number of reference points
%           - final_althyp_index althyp_index used for unwrapped data
%
% ----------------------------------------------------------------------
% File............: ps_network.m
% Version & Date..: 1.7.2.8, 19-OCT-2009
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
%



% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global ps_method unwrap_istep results_id visible_plots
  
  
% ----------------------------------------------------------------------
% 4.1
% Variance Component Estimation (VCE)
% ----------------------------------------------------------------------
  
fprintf(1,'\n');
fprintf(1,'Variance component estimation ....\n');

sig2_est = ps_vce(Npsc,...
                  Nifgs,...
                  max_arc_length,...
                  psc_model,...
                  Btemp,...
                  Bdop,...
                  std_param,...
                  breakpoint,...
                  breakpoint2);



for z = 1:Npsc_selections

  network_flag = 0;
  unwrap_istep = 0;
  
  while network_flag==0

    unwrap_istep = unwrap_istep+1;

    % ----------------------------------------------------------------------
    % 4.2
    % Form network between initial ps
    % ----------------------------------------------------------------------
    
    fprintf(1,'\n');
    fprintf(1,'Forming network ....\n');
    
    Narcs = ps_form_network_sim(Npsc,...
                                 Npsc_selections,...
                                 Nifgs,...
                                 max_arc_length,...
                                 network_method,...
                                 Ncon,...
                                 Nparts,...
                                 z);
    
    
  
    
    % ----------------------------------------------------------------------
    % 4.3
    % Ps estimation
    % ----------------------------------------------------------------------
    
    fprintf(1,'\n');
    fprintf(1,['Ps estimation by ' ps_method ' ....\n']);

    ps_estimation(Btemp,...
                  Bdop,...
                  psc_model,...
                  Nifgs,...
                  Narcs,...
                  Npsc_selections,...
                  std_param,...
                  sig2_est,...
                  breakpoint,...
                  breakpoint2,...
                  z);
    
    
    
    
    % ----------------------------------------------------------------------
    % 4.4
    % Spatial unwrapping
    % ----------------------------------------------------------------------
    
    fprintf(1,'\n');
    fprintf(1,'Spatial unwrapping ....\n');
    
    [Nref,...
     final_althyp_index,...
     network_flag] = ps_spatial_unwrapping_sim(ens_coh_threshold,...
                                               varfac_threshold,...
                                               ref_cn,...
                                               Nifgs,...
                                               Npsc,...
                                               Npsc_selections,...
                                               Nref,...
                                               Narcs,...
                                               Btemp,...
                                               Bdop,...
                                               final_model,...
                                               std_param,...
                                               sig2_est,...
                                               breakpoint,...
                                               breakpoint2,...
                                               network_flag,...
                                               z);
    
  end
end
