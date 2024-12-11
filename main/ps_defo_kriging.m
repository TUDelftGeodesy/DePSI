function ps_defo_kriging(ps_data,Nifgs,Npsc,Npsp,std_param,filenames_output,z)

% Deformation model estimation by Kriging
%
% Input:    - ps_data           ps data (after data reduction)
%           - Nifgs             number of interferograms
%           - Npsc              number of psc
%           - Npsp              number of psp
%           - std_param         a priori standard deviations of
%                               parameters
%           - filenames_output  filenames for output
%           - z                 psc selection
%
% Output:   - (to file only)
%
% ----------------------------------------------------------------------
% File............: ps_defo_kriging.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
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
% v1.7.2.12, Freek van Leijen
% - filenames_output in cells
%


% -----------------------------------------------------
% Initialize
% -----------------------------------------------------

global Nlines Npixels az_spacing r_spacing Npar_max m2ph
global defo_range max_mem_buffer project_id
global ps_eval_method detail_plots fig

Nps = size(ps_data,1);
Nmax = 50;
dx_max = 100000; %[m], FvL changed to high number (100000 ipv 5000)
kriging_method = 4; %Universal Kriging
covariance_model = 3; % gaussian
a0 = 500;
b0 = 0;

ps_az = ps_data(:,2)*az_spacing;
ps_r = ps_data(:,3)*r_spacing;


% --------------------------------------------------------------------
% Read data
% --------------------------------------------------------------------

psc_fid = fopen([project_id '_psc_sel' num2str(z) '.raw'],'r'); 
psc_data = fread(psc_fid,[2*Nifgs+7 Npsc(z)],'double')';
fclose(psc_fid);

psc_array = psc_data(:,1:2);
psc_az = psc_data(:,5);
psc_r = psc_data(:,6);
clear psc_data

psc_defo = NaN(Npsc(z),Nifgs);

psc = find(psc_array(:,2)~=0); % remove isolated PSCs
psc_az_new = psc_az(psc);
psc_r_new = psc_r(psc);
psc_azx = psc_az*az_spacing; % coordinates in meters
psc_rx = psc_r*r_spacing;
psc_azx_new = psc_azx(psc);
psc_rx_new = psc_rx(psc);
clear psc_az psc_r


switch ps_eval_method
  case 'psp'
    
    % --------------------------------------------------------------------
    % Read data
    % --------------------------------------------------------------------

    Npsp_buffer = floor(max_mem_buffer/(4*Nifgs*8));
    if Npsp_buffer>=Npsp(z)
      Npsp_buffer = Npsp(z);
      Npsp_rem_buffer = 0;
      Nbuffers = 1;
    else
      Nbuffers = floor(Npsp(z)/Npsp_buffer);
      Npsp_rem_buffer = rem(Npsp(z),Npsp_buffer);
      if (Npsp_rem_buffer > 0)
        Nbuffers = Nbuffers + 1;
      end
    end
    
    psp_fid = fopen([project_id '_psp_sel' num2str(z) '.raw'],'r'); 
    psp_az = NaN(Npsp(z),1);
    psp_r = NaN(Npsp(z),1);
    psp_grid_az = NaN(Npsp(z),1);
    psp_grid_r = NaN(Npsp(z),1);
    
    for v = 1:Nbuffers
      sl = (v-1)*Npsp_buffer+1;
      el = v*Npsp_buffer;
      if (v==Nbuffers)&&(Npsp_rem_buffer~=0)
        Npsp_buffer = Npsp_rem_buffer;
        el = Npsp(z);
      end
      
      psp_data = fread(psp_fid,[2*Nifgs+4 Npsp_buffer],'double')';
      
      psp_grid_az(sl:el) = psp_data(:,1);
      psp_grid_r(sl:el) = psp_data(:,2);
      psp_az(sl:el) = psp_data(:,3);
      psp_r(sl:el) = psp_data(:,4);
      clear psp_data
    end
    fclose(psp_fid);
    
    psp_azx = psp_az*az_spacing; % coordinates in meters
    psp_rx = psp_r*r_spacing;

    psp_defo_fid = fopen([project_id '_psp_defo_temp_sel' num2str(z) '.raw'],'w');
    psp = (1:Npsp(z))';
    
    
    % ------------------------------------------------------------------
    % Kriging of psc and psp
    % ------------------------------------------------------------------
    
    grid_az = unique(psp_grid_az);
    Ngrid_az = length(grid_az);
    grid_r = unique(psp_grid_r);
    Ngrid_r = length(grid_r);
    
    for v = 1:Nifgs
      
      fprintf(1,'Processing interferogram %g of %g for selection %g\n',v,Nifgs,z);   
      x = [ps_az ps_r m2ph*ps_data(:,Npar_max+6+v)];
      var_defo = var(x(:,3));

      
      % ------------------------------------------------------------------
      % Fit variogram
      % ------------------------------------------------------------------
      
      c10 = 0.8*var_defo;
      c20 = 0.2*var_defo;
      [a,b,c1,c2,emp_vario] = ps_fit_vario(x(:,1),x(:,2),x(:,3),covariance_model,a0,b0,c10,c20);
      c2 = max(c2,c1/100); % to avoid c2~0 (->singular Kriging system!?)
      cv_model = [1 NaN c2;covariance_model a c1];
      
      if strcmp(detail_plots,'y')
        figure(fig);
        print('-dpng',['plots/' project_id '_defo_variogram_ifgs' num2str(v) '.png']);
      end

      for w = 1:length(psc)
        x0 = [psc_azx_new(w) psc_rx_new(w)];
        psc_defo(psc(w),v) = ps_kriging(x,x0,Nmax,dx_max,kriging_method,cv_model);
      end
      
      psp_defo = NaN(Npsp(z),1);
      
      for w = 1:Ngrid_az
        index1 = find(psp_grid_az==w);
        
        for u = 1:Ngrid_r
          index2 = find(psp_grid_r(index1)==u);
          index = index1(index2);
          if ~isempty(index)
            x0 = [psp_azx(index) psp_rx(index)];
            
            psp_defo(index) = ps_kriging(x,x0,Nmax,dx_max,kriging_method,cv_model);
          end
        end
      end
      fwrite(psp_defo_fid,psp_defo,'double');
    end

    
  case 'whole'
    
    % ------------------------------------------------------------------
    % Kriging of whole image
    % ------------------------------------------------------------------
    
    Sgrid_az = round(200/az_spacing);
    Sgrid_r = round(200/r_spacing);
    Ngrid_az = ceil(Nlines/Sgrid_az);
    Ngrid_r = ceil(Npixels/Sgrid_r);
    
    ps_az = ps_data(:,2)*az_spacing; %in meters
    ps_r = ps_data(:,3)*r_spacing;

    for v = 1:Nifgs
      
      fprintf(1,'Processing interferogram %g of %g for selection %g\n',v,Nifgs,z);   
      
      defo_fid = fopen([char(filenames_output(v)) '_defo_sel' num2str(z) '.raw'],'w');
      
      x = [ps_az ps_r m2ph*ps_data(:,Npar_max+6+v)];
      var_defo = var(x(:,3));
      

      % ------------------------------------------------------------------
      % Fit variogram
      % ------------------------------------------------------------------
      
      c10 = 0.8*var_defo;
      c20 = 0.2*var_defo;
      [a,b,c1,c2,emp_vario] = ps_fit_vario(x(:,1),x(:,2),x(:,3),covariance_model,a0,b0,c10,c20);
      c2 = max(c2,c1/100); % to avoid c2~0 (->singular Kriging system!?)
      cv_model = [1 NaN c2;covariance_model a c1];
      
      if strcmp(detail_plots,'y')
        figure(fig);
        print('-dpng',['plots/' project_id '_defo_variogram_ifgs' num2str(v) '.png']);
      end
      
      for w = 1:Ngrid_az
        az_start = (w-1)*Sgrid_az+1;
        az_stop = min(w*Sgrid_az,Nlines);
        Sgrid_az2 = az_stop-az_start+1;
        defo_array = NaN(Sgrid_az2,Npixels);
        
        for u = 1:Ngrid_r
          r_start = (u-1)*Sgrid_r+1;
          r_stop = min(u*Sgrid_r,Npixels);
          Sgrid_r2 = r_stop-r_start+1;
          
          x0 = NaN(Sgrid_az2*Sgrid_r2,2);
          x0(:,1) = kron(ones(Sgrid_r2,1),(az_start:az_stop)'*az_spacing);
          x0(:,2) = kron((r_start:r_stop)'*r_spacing,ones(Sgrid_az2,1));
          
          sol = ps_kriging(x,x0,Nmax,dx_max,kriging_method,cv_model);
          defo_array(:,r_start:r_stop) = reshape(sol,Sgrid_az2,Sgrid_r2);
          %defo in radians
        end
        
        % -------------------------------------------------------
        % Store psc values
        % -------------------------------------------------------
        
        index = find((psc_az_new>=az_start)&(psc_az_new<=az_stop));
        psc_az_buffer = psc_az_new(index)-az_start+1;
        psc_r_buffer = psc_r_new(index);
        psc_defo(psc(index),v) = diag(defo_array(psc_az_buffer,psc_r_buffer));
        
        fwrite(defo_fid,defo_array','single');
      end
      fclose(defo_fid);
    end
  
  otherwise
    error('You specified a wrong ps_eval_method..');
end


% -----------------------------------------------------------------
% Write psc values to file
% -----------------------------------------------------------------

psc_defo_fid = fopen([project_id '_psc_defo_temp_sel' num2str(z) '.raw'],'w');
fwrite(psc_defo_fid,psc_defo','double');
fclose(psc_defo_fid);


