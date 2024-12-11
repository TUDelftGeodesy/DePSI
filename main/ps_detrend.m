function detrend_param = ps_detrend(Nifgs,Npsc,Npsp,Npsc_selections,filenames_output)

% Function to remove a trend from the interferograms based on the psc
%
% Input:    - Nifgs              number of interferograms
%           - Npsc               number of psc's
%           - Npsp               number of psp's
%           - Npsc_selections    number of psc's selections
%           - filenames_output   filenames of output
%
% Output:   - detrend_param      parameters of the estimated plane
%                                for each interferogram
%
% ----------------------------------------------------------------------
% File............: ps_detrend.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Author..........: Gini Ketelaar
%                   Freek van Leijen
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
% - filenames_output in cells
% v1.7.2.16, Freek van Leijen
% - increased psp buffer
% v1.7.7.3, Freek van Leijen
% - saving and loading of previously estimated trends
%


% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global project_id Npar_max detail_plots visible_plots
global fig ps_eval_method max_mem_buffer Nlines Npixels

maxit=500; % maximum number of iterations
pcrej=0.5; % fraction of ehat above critical value to be rejected per iteration (to speed up)
ehatmax=100/180*pi; % maximum residual (absolute) 
% datasnooping is to exclude topography and defo fringes
% of course only works if majority of area is not subsiding and flat


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
  
  psc_az = psc_data(:,5);
  psc_r = psc_data(:,6);
  psc_phase = psc_data(:,7:6+Nifgs);
  
  psc_results_fid = fopen([project_id '_psc_results_sel' num2str(z) '.raw'],'r'); 
  psc_results_data = fread(psc_results_fid,[3*Nifgs+Npar_max+2 Npsc(z)],'double')';
  fclose(psc_results_fid);
  
  psc_acheck = psc_results_data(:,Nifgs+Npar_max+3:2*Nifgs+Npar_max+2);
  clear psc_results_data
  
  psc_phase_unw = 2*pi*psc_acheck + psc_phase; % unwrap
  
  
  % ----------------------------------------------------------------------
  % Estimate a plane for each interferogram
  % ----------------------------------------------------------------------

  if strcmp(detail_plots,'y')
    fig = fig + 1;  
    figure(fig);hold on
    if strcmp(visible_plots,'n')
      set(gcf,'visible','off');
    end
  end
  
  detrend_param = NaN(3,Nifgs);
  
  for v=1:Nifgs
    
    if strcmp(detail_plots,'y')
      figure(fig);
      subplot(2,3,1);scatter(psc_r,psc_az,5,psc_phase(:,v));set(gca,'Clim',[-pi pi]);colorbar; title('phase, wrapped');
      subplot(2,3,2);scatter(psc_r,psc_az,5,psc_phase_unw(:,v));set(gca,'Clim',[-4*pi 4*pi]);colorbar; title('phase, unwrapped');
    end
    
    ehat=inf; 
    ind=find(~isnan(psc_phase_unw(:,v)));
    A=[psc_az(ind) psc_r(ind) ones(size(ind))];
    y=psc_phase_unw(ind,v);
    it=1;
    while and(max(ehat)>ehatmax,it<maxit)
      xhat=A\y;
      ehat=y-A*xhat;
      indrej=find(abs(ehat)>ehatmax);
      [temp indtemp]=sort(abs(ehat(indrej)));
      nrej=ceil(pcrej*numel(temp));
      indrej=indrej(indtemp(end-nrej+1:end));
      A(indrej,:)=[];
      y(indrej)=[];
      it=it+1;
    end
    
    detrend_param(:,v) = xhat;
    
    if strcmp(detail_plots,'y')
      
      % ----------------------------------------------------------------
      % Plot results
      % ----------------------------------------------------------------
      
      % calculate updates phase values
      temp_new = psc_phase_unw(:,v)-[psc_az psc_r  ones(size(psc_az))]*xhat;
      % wrap
      temp_new = mod((temp_new+pi),2*pi)-pi;
      
      figure(fig);
      subplot(2,3,3);scatter(psc_r,psc_az,5,mod([psc_az psc_r  ones(size(psc_az))]*xhat+pi,2*pi)-pi);set(gca,'Clim',[-pi pi]);colorbar;title('correction, wrapped');
      subplot(2,3,4);scatter(A(:,2),A(:,1),5,y-A*xhat);set(gca,'Clim',[-4*pi 4*pi]);colorbar;title('residuals accepted observations');
      subplot(2,3,5);scatter(psc_r,psc_az,5,psc_phase_unw(:,v)-[psc_az psc_r  ones(size(psc_az))]*xhat);title('corrected, unwrapped');set(gca,'Clim',[-4*pi 4*pi]);colorbar;
      subplot(2,3,6);scatter(psc_r,psc_az,5,temp_new);title('corrected, wrapped');set(gca,'Clim',[-pi pi]);colorbar;
      
      clear temp_new;
    end
  end
  

  if exist([project_id '_trend_param_sel' num2str(z) '.mat']);
    load([project_id '_trend_param_sel' num2str(z) '.mat']);
  else
    save([project_id '_trend_param_sel' num2str(z) '.mat'],'detrend_param');
  end

  % ----------------------------------------------------------------
  % Adapt phase psc
  % ----------------------------------------------------------------
  
  psc_trend = [psc_az psc_r  ones(size(psc_az))]*detrend_param;
  psc_data(:,7:6+Nifgs) = psc_data(:,7:6+Nifgs) - psc_trend;
  psc_data(:,7:6+Nifgs) = mod(psc_data(:,7:6+Nifgs)+pi,2*pi)-pi;
  
  psc_fid = fopen([project_id '_psc_3detr_sel' num2str(z) '.raw'],'w'); 
  fwrite(psc_fid,psc_data','double');
  fclose(psc_fid);

  psc_trend_fid = fopen([project_id '_psc_trend_sel' num2str(z) '.raw'],'w');
  fwrite(psc_trend_fid,psc_trend','double');
  fclose(psc_trend_fid);

  
  % ----------------------------------------------------------------
  % Adapt phase psp/original interferograms
  % ----------------------------------------------------------------
  
  switch ps_eval_method
    case 'psp'
      Npsp_buffer = floor(max_mem_buffer/(10*Nifgs*8));
      if (Npsp_buffer>=Npsp(1))
        Npsp_buffer = Npsp(1);
        Npsp_rem_buffer = 0;
        Nbuffers = 1;
      else
        Nbuffers = floor(Npsp(1)/Npsp_buffer);
        Npsp_rem_buffer = rem(Npsp(1),Npsp_buffer);
        if (Npsp_rem_buffer > 0)
          Nbuffers = Nbuffers + 1;
        end
      end
      
      Npsp_sofar = 0;
      psp_fid = fopen([project_id '_psp_sel' num2str(z) '.raw'],'r'); 
      psp_out_fid = fopen([project_id '_psp_3detr_sel' num2str(z) '.raw'],'w'); 
      psp_trend_fid = fopen([project_id '_psp_trend_sel' num2str(z) '.raw'],'w');

      for v = 1:Nbuffers
        if (v==Nbuffers)&&(Npsp_rem_buffer~=0)
          Npsp_buffer = Npsp_rem_buffer;
        end
        
        psp_data = fread(psp_fid,[2*Nifgs+4 Npsp_buffer],'double')';
        psp_az = psp_data(:,3);
        psp_r = psp_data(:,4);
        
        psp_trend = [psp_az psp_r ones(size(psp_az))]*detrend_param;
        psp_data(:,5:4+Nifgs) = psp_data(:,5:4+Nifgs) - psp_trend;
        psp_data(:,5:4+Nifgs) = mod(psp_data(:,5:4+Nifgs)+pi,2*pi)-pi;
        
        fwrite(psp_out_fid,psp_data','double');
        fwrite(psp_trend_fid,psp_trend','double');
        Npsp_sofar = Npsp_sofar+Npsp_buffer;
      end
      fclose(psp_fid);
      fclose(psp_out_fid);
      fclose(psp_trend_fid);
      
    case 'whole'
      
      % -----------------------------------------------------------------
      % Write estimated trend to file
      % -----------------------------------------------------------------
      
      for v = 1:Nifgs
        [trend_r,trend_az] = meshgrid((1:Npixels)*detrend_param(2,v),(1:Nlines)*detrend_param(1,v));
        
        trend = trend_az+trend_r+repmat(detrend_param(3,v),Nlines,Npixels);
        
        trend_fid = fopen([char(filenames_output(v)) '_trend_sel' ...
                          num2str(z) '.raw'],'w');
        fwrite(trend_fid,trend','single');
        fclose(trend_fid);
      end
      
    end
end





