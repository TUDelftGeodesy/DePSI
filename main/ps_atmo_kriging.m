function ps_atmo_kriging(Nifgs,Npsc,Npsc_selections,Npsp,filenames_output)
  
% Atmosphere (APS) estimation by Kriging
%
% Input:    - Nifgs            number of interferograms
%           - Npsc             number of psc's
%           - Npsc_selections  number of psc_selections
%           - Npsp             number of psp's
%           - filenames_output filenames for output
%
% ----------------------------------------------------------------------
% File............: ps_atmo_kriging.m
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
% v1.7.2.16, Freek van Leijen
% - psp radar coordinates in single precision to reduce memory
% vsvn, Freek van Leijen
% - Nmax changed to 100
%



% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global max_mem_buffer Nlines Npixels az_spacing r_spacing Npar_max
global fig ps_eval_method project_id m2ph
global detail_plots visible_plots atmo_range

Nmax = 100; % FvL changed to higher number (100 ipv 20)
dx_max = 100000; %[m] , FvL changed to high number (100000 ipv 10000)
kriging_method = 4; %Universal Kriging
covariance_model = 3; % gaussian
a0 = 1000;
b0 = 0;


for z = 1:Npsc_selections

  % --------------------------------------------------------------------
  % Read data
  % --------------------------------------------------------------------
  
  psc_fid = fopen([project_id '_psc_sel' num2str(z) '.raw'],'r'); 
  psc_data = fread(psc_fid,[2*Nifgs+7 Npsc(z)],'double')';
  fclose(psc_fid);
  
  psc_array = psc_data(:,1:2);
  psc = find(psc_array(:,2)~=0); 
  psc_az = psc_data(psc,5);
  psc_r = psc_data(psc,6);
  clear psc_data
  Npsc_new = length(psc_az);

  psc_results_fid = fopen([project_id '_psc_results_sel' num2str(z) '.raw'],'r'); 
  psc_data = fread(psc_results_fid,[3*Nifgs+Npar_max+2 Npsc(z)],'double')';
  fclose(psc_results_fid);
  
  psc_atmo_master = m2ph*psc_data(psc,2); %second column of psc_param, dangerous...
  clear psc_data
  
  atmo_fid = fopen([project_id '_psc_atmo_filt_sel' num2str(z) '.raw'],'r'); 
  atmo_estimates = fread(atmo_fid,[Nifgs Npsc(z)],'double')';
  fclose(atmo_fid);
  atmo_estimates = atmo_estimates(psc,:);
  atmo_estimates = [atmo_estimates psc_atmo_master]; %master atmo
  
  psc_atmo = NaN(Npsc(z),Nifgs); 
  psc_atmo_master = NaN(Npsc(z),1);
  
  psc_azx = psc_az*az_spacing; % coordinates in meters
  psc_rx = psc_r*r_spacing;

  switch ps_eval_method
    case 'psp'
      Npsp_buffer = floor(max_mem_buffer/(8*Nifgs*8));
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
      psp_az = NaN(Npsp(z),1,'single'); %changed to single 
      psp_r = NaN(Npsp(z),1,'single');
      psp_grid_az = NaN(Npsp(z),1,'single');
      psp_grid_r = NaN(Npsp(z),1,'single');
      
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
      clear psp_az psp_r
      
      psp_atmo_fid = fopen([project_id '_psp_atmo_temp_sel' num2str(z) '.raw'],'w');

      
      % ------------------------------------------------------------------
      % Kriging of psc and psp
      % ------------------------------------------------------------------
      
      grid_az = unique(psp_grid_az);
      Ngrid_az = length(grid_az);
      grid_r = unique(psp_grid_r);
      Ngrid_r = length(grid_r);
      
      for v = 1:Nifgs+1 % +1 for master

	if v<=Nifgs
	  fprintf(1,'Processing atmospheric phase screen %g of %g for selection %g\n',v,Nifgs,z);   
	else
	  fprintf(1,'Processing master atmospheric phase screen\n');   
	end
	
        x = [psc_azx psc_rx atmo_estimates(:,v)];
        var_atmo = var(x(:,3));
        
        
        % ------------------------------------------------------------------
        % Fit variogram
        % ------------------------------------------------------------------
        
        c10 = 0.8*var_atmo;
        c20 = 0.2*var_atmo;

        [a,b,c1,c2,emp_vario] = ps_fit_vario(x(:,1),x(:,2),x(:,3),covariance_model,a0,b0,c10,c20);
        c2 = max(c2,c1/100); % to avoid c2~0 (->singular Kriging system!?)
        cv_model = [1 NaN c2;covariance_model a c1];
        
        if strcmp(detail_plots,'y')
          figure(fig);
          print('-dpng',['plots/' project_id '_atmo_variogram_ifgs' num2str(v) '.png']);
        end
        
        for w = 1:Npsc_new

          x0 = [psc_azx(w) psc_rx(w)];
	  if v<=Nifgs
	    psc_atmo(psc(w),v) = ps_kriging(x,x0,Nmax,dx_max,kriging_method,cv_model);
	  else
	    psc_atmo_master(psc(w)) = ps_kriging(x,x0,Nmax,dx_max,kriging_method,cv_model);
	  end
        end
	
        psp_atmo = NaN(Npsp(z),1);
        
        for w = 1:Ngrid_az
          index1 = find(psp_grid_az==grid_az(w));
          
          for u = 1:Ngrid_r

            index2 = find(psp_grid_r(index1)==grid_r(u));
            index = index1(index2);
            if ~isempty(index)
              x0 = [psp_azx(index) psp_rx(index)];
              
              psp_atmo(index) = ps_kriging(x,x0,Nmax,dx_max,kriging_method,cv_model);
            end
          end
        end
	%%%psp_atmo = zeros(size(psp_atmo)); %%%TEMP FIX

	if v<=Nifgs
	  fwrite(psp_atmo_fid,psp_atmo,'double');
	else
	  psp_atmo_master_fid = fopen([project_id '_psp_atmo_master1_sel' num2str(z) '.raw'],'w');
	  fwrite(psp_atmo_master_fid,psp_atmo,'double');
	  fclose(psp_atmo_master_fid);
	end
      end
          
    case 'whole'
      
      % ------------------------------------------------------------------
      % Kriging of whole image
      % ------------------------------------------------------------------
      
      Sgrid_az = round(200/az_spacing);
      Sgrid_r = round(200/r_spacing);
      Ngrid_az = ceil(Nlines/Sgrid_az);
      Ngrid_r = ceil(Npixels/Sgrid_r);
      
      for v = 1:Nifgs+1 % +1 for master

	if v<=Nifgs
	  fprintf(1,'Processing atmospheric phase screen %g of %g for selection %g\n',v,Nifgs,z);   
	else
	  fprintf(1,'Processing master atmospheric phase screen\n');   
	end
        
	if v<=Nifgs
	  atmo_fid = fopen([char(filenames_output(v)) '_atmo_sel' num2str(z) '.raw'],'w');
	else
	  atmo_fid = fopen([project_id '_atmo_master1_sel' num2str(z) '.raw'],'w');
	end
	
	x = [psc_azx psc_rx atmo_estimates(:,v)];
        var_atmo = var(x(:,3));


        % ------------------------------------------------------------------
        % Fit variogram
        % ------------------------------------------------------------------
        
        c10 = 0.8*var_atmo;
        c20 = 0.2*var_atmo;
        [a,b,c1,c2,emp_vario] = ps_fit_vario(x(:,1),x(:,2),x(:,3),covariance_model,a0,b0,c10,c20);
        c2 = max(c2,c1/100); % to avoid c2~0 (->singular Kriging system!?)
        cv_model = [1 NaN c2;covariance_model a c1];
        
        if strcmp(detail_plots,'y')
          figure(fig);
          print('-dpng',['plots/' project_id '_atmo_variogram_ifgs' num2str(v) '.png']);
        end
        
        for w = 1:Ngrid_az
          az_start = (w-1)*Sgrid_az+1;
          az_stop = min(w*Sgrid_az,Nlines);
          Sgrid_az2 = az_stop-az_start+1;
          atmo_array = NaN(Sgrid_az2,Npixels);
          
          for u = 1:Ngrid_r
            r_start = (u-1)*Sgrid_r+1;
            r_stop = min(u*Sgrid_r,Npixels);
            Sgrid_r2 = r_stop-r_start+1;
            
            x0 = NaN(Sgrid_az2*Sgrid_r2,2);
            x0(:,1) = kron(ones(Sgrid_r2,1),(az_start:az_stop)'*az_spacing);
            x0(:,2) = kron((r_start:r_stop)'*r_spacing,ones(Sgrid_az2,1));
            
            sol = ps_kriging(x,x0,Nmax,dx_max,kriging_method,cv_model);
            atmo_array(:,r_start:r_stop) = reshape(sol,Sgrid_az2,Sgrid_r2);
          end

          
          % -------------------------------------------------------
          % Store psc values
          % -------------------------------------------------------
          
          index = find((psc_az>=az_start)&(psc_az<=az_stop));
          psc_az_buffer = psc_az(index)-az_start+1;
          psc_r_buffer = psc_r(index);
          psc_atmo(psc(index),v) = diag(atmo_array(psc_az_buffer,psc_r_buffer));

          
          % -------------------------------------------------------
          % Output to file
          % -------------------------------------------------------

          fwrite(atmo_fid,atmo_array','single'); 
          
        end
        fclose(atmo_fid);
      end
    
    otherwise
      error('You specified a wrong ps_eval_method..');
  end
  
  
  % -----------------------------------------------------------------
  % Write psc values to file
  % -----------------------------------------------------------------

  %%%psc_atmo = zeros(size(psc_atmo)); %%%TEMP FIX
  
  psc_atmo_fid = fopen([project_id '_psc_atmo_temp_sel' num2str(z) '.raw'],'w');
  fwrite(psc_atmo_fid,psc_atmo','double');
  fclose(psc_atmo_fid);

  psc_atmo_master_fid = fopen([project_id '_psc_atmo_master1_sel' num2str(z) '.raw'],'w');
  fwrite(psc_atmo_master_fid,psc_atmo_master,'double');
  fclose(psc_atmo_master_fid);

  psc_atmo_master2 = mean(psc_atmo,2);
  psc_atmo_master_fid = fopen([project_id '_psc_atmo_master2_sel' num2str(z) '.raw'],'w');
  fwrite(psc_atmo_master_fid,psc_atmo_master2,'double');
  fclose(psc_atmo_master_fid);

end  

clear psp_azx psp_rx


% -------------------------------------------------------------------
% Subtract the atmospheric delays
% -------------------------------------------------------------------

for z = 1:Npsc_selections
  psc_fid = fopen([project_id '_psc_sel' num2str(z) '.raw'],'r'); 
  psc_data = fread(psc_fid,[2*Nifgs+7 Npsc(z)],'double')';
  
  psc_array = psc_data(:,1:2);
  psc_grid_az = psc_data(:,3);
  psc_grid_r = psc_data(:,4);
  psc_az = psc_data(:,5);
  psc_r = psc_data(:,6);
  psc_phase = psc_data(:,7:Nifgs+6);
  psc_h2ph = psc_data(:,Nifgs+7:2*Nifgs+6);
  psc_amp_disp = psc_data(:,2*Nifgs+7);
  clear psc_data
  fclose(psc_fid);
  
  psc_fid = fopen([project_id '_psc_4atmo_sel' num2str(z) '.raw'],'w'); 

  psc_atmo_fid = fopen([project_id '_psc_atmo_sel' num2str(z) '.raw'],'w');

  psc_atmo_temp_fid = fopen([project_id '_psc_atmo_temp_sel' num2str(z) '.raw'],'r'); 
  psc_atmo = fread(psc_atmo_temp_fid,[Nifgs Npsc(z)],'double')';
  
  psc_phase = mod(psc_phase-psc_atmo+pi,2*pi)-pi;
  
  psc_total = [psc_array psc_grid_az psc_grid_r psc_az psc_r psc_phase psc_h2ph psc_amp_disp];
  
  fwrite(psc_fid,psc_total','double');
  fclose(psc_fid);
  
  frewind(psc_atmo_fid); % fixed bug in version 1.5
  fwrite(psc_atmo_fid,psc_atmo','double');
  fclose(psc_atmo_fid);
end


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
    psp_atmo_temp_fid = NaN(Npsc_selections,1);
    psp_atmo_fid = NaN(Npsc_selections,1); %write again in
                                                 %correct order
    psp_fid = NaN(Npsc_selections,1);
    psp_fid_out = NaN(Npsc_selections,1);
    for z = 1:Npsc_selections
        psp_fid(z) = fopen([project_id '_psp_sel' num2str(z) '.raw'],'r'); 
        psp_fid_out(z) = fopen([project_id '_psp_4atmo_sel' num2str(z) '.raw'],'w'); 
        psp_atmo_temp_fid(z) = fopen([project_id '_psp_atmo_temp_sel' num2str(z) '.raw'],'r'); 
        psp_atmo_fid(z) = fopen([project_id '_psp_atmo_sel' num2str(z) '.raw'],'w');
	psp_atmo_master_fid(z) = fopen([project_id '_psp_atmo_master2_sel' num2str(z) '.raw'],'w');
    end
    
    for v = 1:Nbuffers
        if (v==Nbuffers)&&(Npsp_rem_buffer~=0)
            Npsp_buffer = Npsp_rem_buffer;
        end
        
        for z = 1:Npsc_selections
            psp_data = fread(psp_fid(z),[2*Nifgs+4 Npsp_buffer],'double')';
            
            psp_grid_az = psp_data(:,1);
            psp_grid_r = psp_data(:,2);
            psp_az = psp_data(:,3);
            psp_r = psp_data(:,4);
            psp_phase_orig = psp_data(:,5:Nifgs+4);
            psp_h2ph = psp_data(:,Nifgs+5:2*Nifgs+4);
            clear psp_data
            
            psp_atmo = NaN(Npsp_buffer,Nifgs);
            for w = 1:Nifgs
                fseek(psp_atmo_temp_fid(z),(w-1)*8*Npsp(z)+Npsp_sofar*8,-1);
                psp_atmo(:,w) = fread(psp_atmo_temp_fid(z),[Npsp_buffer 1],'double');
            end
            psp_phase = mod(psp_phase_orig-psp_atmo+pi,2*pi)-pi;
            
            psp_total = [psp_grid_az psp_grid_r psp_az psp_r psp_phase psp_h2ph];
            fwrite(psp_fid_out(z),psp_total','double');
            fwrite(psp_atmo_fid(z),psp_atmo','double');
            
	    psp_atmo_master2 = mean(psp_atmo,2);
	    fwrite(psp_atmo_master_fid(z),psp_atmo_master2,'double');
	    
        end
        Npsp_sofar = Npsp_sofar+Npsp_buffer;
    end
end

fclose('all');
