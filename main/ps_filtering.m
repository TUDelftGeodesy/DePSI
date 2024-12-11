function ps_filtering(Btemp,ts_atmo_filter,ts_atmo_filter_length,Nifgs,Npsc,Npsc_selections)

% Function to estimate the master APS by averaging, substract the
% master APS, apply a low-pass filter to the time series to
% remove the non-linear deformation. What remains is the APS per
% slave. The master APS and slave APS's are combined again.
%
% Input:  - Btemp                  temporal baselines [year]
%         - ts_atmo_filter         atmo filter: block, triangle,
%                                  gaussian or prediction
%         - ts_atmo_filter_length  atmo filter length [year]
%         - filter_length          length filter for low-pass filtering
%                                  non-linear deformation [year]
%         - Nifgs                  number of interferograms
%         - Npsc                   number of psc's
%         - Npsc_selections        number of psc selections
%
% Output: - atmo_estimates         atmosphere estimates per psc
% (to file)
%
%
% This software is developed by Delft University of Technology and is
% intended for scientific use only. Applications for commercial use are
% prohibited.
%  
% ----------------------------------------------------------------------
% File............: ps_filtering.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Author..........: Freek van Leijen
%                   Astrid Humme
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
% - bug fixed regarding master atmosphere
% v1.7.2.8, Astrid Humme
% - switch for various filters



% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global project_id Npar_max ps_method m2ph fig


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
  clear psc_data
  
  psc_results_fid = fopen([project_id '_psc_results_sel' num2str(z) '.raw'],'r'); 
  psc_data = fread(psc_results_fid,[3*Nifgs+Npar_max+2 Npsc(z)],'double')';
  % the transpose at the end is obviously very important...
  fclose(psc_results_fid);
  
  psc_atmo_master = psc_data(:,2); %second column of psc_param, dangerous...
  psc_phase_res = psc_data(:,Npar_max+3:Npar_max+Nifgs+2);
  clear psc_data
  
  % ----------------------------------------------------------------------
  % Remove the isolated points, initialize
  % ----------------------------------------------------------------------

  psc = find(psc_array(:,2)~=0); 
  psc_phase_res = psc_phase_res(psc,:);
  psc_atmo_master = psc_atmo_master(psc);
  Npsc_new = size(psc_phase_res,1);
  Npsc_tot = Npsc(z);
   


  % ----------------------------------------------------------------------
  % Estimate master ATMO and subtract
  % ----------------------------------------------------------------------
  
  %THIS IS WRONG: (FvL, 04-10-09)
  %switch ps_method 
  %  case 'perio'
  %    atmo_master = mean(psc_phase_res,2);
  %    psc_phase_res = psc_phase_res-repmat(atmo_master,1,Nifgs);
  %  case {'ils','boot'}
  %    index = find(~isnan(psc_atmo_master));
  %    if ~isempty(index)
  %        atmo_master = m2ph*psc_atmo_master;
  %    else
  %        error('Something went wrong with the estimated master atmosphere');
  %    end
  %    %psc_phase_res is already without master atmosphere
  %end

  index = find(~isnan(psc_atmo_master));
  if ~isempty(index)
    atmo_master = m2ph*psc_atmo_master;
  else
    error('Something went wrong with the estimated master atmosphere');
  end
  %psc_phase_res is already without master atmosphere

  
  % ----------------------------------------------------------------------
  % Construct filter
  % ----------------------------------------------------------------------            
  
  time_span = max(Btemp)-min(Btemp);
  filter_nl = zeros(1,ceil(time_span*1000)+1);                
  
  switch ts_atmo_filter
    
   case 'block'
    
    filter_nl(1:round(0.5*ts_atmo_filter_length*1000)+1) = 1;
    
   case 'triangle'
    
    filter_nl(1:round(0.5*ts_atmo_filter_length*1000)+1) = fliplr(0:1/round(0.5*ts_atmo_filter_length*1000):1);
    % factor 0.5 because of mirrored filter
    
   case 'gaussian'
    
    filter_length_ts = ts_atmo_filter_length / 3; % 3*std
    filter_nl = 1/(sqrt(2*pi)*0.5*filter_length_ts*1000) * exp (-0.5*((0:(ceil(time_span*1000)+1))/(0.5*filter_length_ts*1000)).^2);
    
   case 'prediction'
    
    error(['Atmosphere filtering by least-squares prediction is not implemented' ...
	   ' yet']);
    
  end
  
  % ----------------------------------------------------------------------
  % Estimate non-linear deformation
  % ----------------------------------------------------------------------
  
  non_linear = zeros(Npsc_new,Nifgs);
  for v = 1:Nifgs
    weights = filter_nl(floor(abs(Btemp'-Btemp(v))*1000)+1);
    non_linear(:,v) = (psc_phase_res*weights')/sum(weights);
  end
  
  
  
  % ----------------------------------------------------------------------
  % Estimate slave atmosphere and construct atmosphere per interferogram
  % ----------------------------------------------------------------------
  
  atmo_slave = psc_phase_res - non_linear; % = actually minus atmo_slave
  atmo_estimates = NaN(Npsc_tot,Nifgs);
  atmo_estimates(psc,:) = atmo_slave + repmat(atmo_master,1,Nifgs);
  % atmo_master+atmo_slave because you actually estimate -atmo_slave

  % ----------------------------------------------------------------------
  % Write results to file
  % ----------------------------------------------------------------------

  atmo_fid = fopen([project_id '_psc_atmo_filt_sel' num2str(z) '.raw'],'w');
  fwrite(atmo_fid,atmo_estimates','double');
  fclose(atmo_fid);
  
  % ----------------------------------------------------------------------
  % Create random plots to visualize filtering of atmosphere
  % ----------------------------------------------------------------------

  psc_random = unique(ceil(rand(10,1)*length(psc)));

  for v = 1:length(psc_random)
    fig = fig+1;
    figure(fig);hold on;
    set(gcf,'visible','off');
    h1=plot(Btemp,psc_phase_res(psc_random(v),:),'r','linewidth',2);
    plot(Btemp,psc_phase_res(psc_random(v),:),'ro','linewidth',2);
    h2=plot(Btemp,atmo_slave(psc_random(v),:),'b','linewidth',2);
    plot(Btemp,atmo_slave(psc_random(v),:),'bo','linewidth',2);
    h3=plot(Btemp,non_linear(psc_random(v),:),'g','linewidth',2);
    plot(Btemp,non_linear(psc_random(v),:),'go','linewidth',2);
    grid on;
    xlabel('Temporal baseline [year]');
    ylabel('Phase [rad]');
    title(['Atmospheric filter, PSC ' num2str(psc(psc_random(v))) ' (randomly selected)']);
    legend([h1 h2 h3],[{'Residual phase'},{'Slave atmosphere'}, ...
		       {'Non-linear deformation'}],'Location',...
	   'SouthOutside','Orientation','horizontal');

    print('-dpng',['plots/atmospheric_filter/' project_id ...
                   '_atmospheric_filter_sel' num2str(z) ...
		   '_psc' num2str(psc(psc_random(v))) '.png']);

  end
  close all
end

