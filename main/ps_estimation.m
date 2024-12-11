function ps_estimation(Btemp,Bdop,model,Nifgs,Narcs,Npsc_selections,std_param,sig2_est,breakpoint,breakpoint2,z)

% Function to estimate the ambiguities using LAMBDA and to
% calculate the residual phases
%
% Input:   - Btemp               temporal baselines [year]
%          - Bdop                Doppler baselines [Hz]
%          - model               vector with model indices to be
%                                evaluated (see help model_definitions
%                                for more information)
%          - Nifgs               number of interferograms
%          - Narcs               number of arcs
%          - Npsc_selections     number of psc selections
%          - std_param           standard deviations for pseudo-
%                                observations
%          - sig2_est            estimated variance components
%          - breakpoint          breakpoint in case of double
%                                linear model
%          - breakpoint2         second breakpoint, should be larger
%                                than first breakpoint 
%          - z                   current psc selection
%
% ----------------------------------------------------------------------
% File............: ps_estimation.m
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

global project_id ps_method

step1_orig = 1; % [m]
step2_orig = 0.0001; % [m], = 0.1 mm




% ----------------------------------------------------------------------
% Read data
% ----------------------------------------------------------------------

dpsc_fid = fopen([project_id '_dpsc_sel' num2str(z) '.raw'],'r'); 
dpsc_data = fread(dpsc_fid,[2*Nifgs+2 Narcs(z)],'double')';
fclose(dpsc_fid);

dpsc_phase = dpsc_data(:,3:Nifgs+2);
dpsc_h2ph = dpsc_data(:,Nifgs+3:2*Nifgs+2);
clear dpsc_data


% ----------------------------------------------------------------------
% Psc estimation
% ----------------------------------------------------------------------

switch ps_method
  case 'perio'
    
    [dpsc_param,dpsc_acheck,dpsc_model_info,dpsc_ens_coh] = ps_periodogram_estimation(dpsc_h2ph,Btemp,dpsc_phase,model,std_param,step1_orig,step2_orig);

    dpsc_varfac_best = NaN(Narcs(z),1);
    
  case {'ils','boot'}
    
    [dpsc_param,dpsc_acheck,dpsc_model_info,dpsc_varfac_best,dpsc_ens_coh] = ps_ils_and_bootstrap_estimation(dpsc_h2ph,Btemp,Bdop,dpsc_phase,model,std_param,sig2_est,breakpoint,breakpoint2);
    
    if isempty(find(~isnan(dpsc_model_info(:,1))));
      error('Something went wrong in the estimation of the network. Most likely the stochastic model used is incorrect.');
    end
    
  otherwise
    
    error('Please select a valid ps method: ''period'', ''ils'' or ''boot''');
    
end

% ----------------------------------------------------------------------
% Write results to file
% ----------------------------------------------------------------------

dpsc_total = [dpsc_param dpsc_ens_coh dpsc_acheck dpsc_model_info dpsc_varfac_best];
dpsc_fid = fopen([project_id '_dpsc_results_sel' num2str(z) '.raw'],'w');
fwrite(dpsc_fid,dpsc_total','double');
fclose(dpsc_fid);
