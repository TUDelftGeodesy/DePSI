function ps_check_results_sim(Nifgs,Nps_atmo,Nps_defo,Npsc,Npsp,Nref,Npsc_selections,final_model,defo_model_flag)

%
% Input:    
%    
% Output:   
%    
% ----------------------------------------------------------------------
% File............: ps_check_results_sim.m
% Version & Date..: 1.7.2.8, 19-OCT-2009
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



global project_id orbit
global Npar_max Nlines Npixels az_spacing r_spacing m2ph max_mem_buffer


% ----------------------------------------------------------------------
% Create the correct annotation
% ----------------------------------------------------------------------
 
[Npar,par_index,covar_index] = ps_model_definitions('Npar',final_model);
[ps_annotation,par_index,covar_index] = ps_model_definitions('annotation',final_model);

switch defo_model_flag
  case 'yes'
    results_id = [{'atmo'};{'defo'}];
    Nps = [Nps_atmo Nps_defo];
    Nresults = 2;
  case 'no'
    results_id = {'atmo'};
    Nps = Nps_atmo;
    Nresults = 1;
end


for z = 1:Npsc_selections

  for w = 1:Nresults

    % ----------------------------------------------------------------------
    % Determine buffersize
    % ----------------------------------------------------------------------
    
    Nps_buffer = floor(max_mem_buffer/(5*(2*Nifgs+2)*8));
    
    if (Nps_buffer>=Nps(z,w))
      Nps_buffer = Nps(z,w);
      Nbuffers = 1;
      Nps_rem_buffer = 0;
    else
      Nbuffers = floor(Nps,w(z)/Nps_buffer);
      Nps_rem_buffer = rem(Nps(z,w),Nps_buffer);
      if (Nps_rem_buffer > 0)
        Nbuffers = Nbuffers + 1;
      end
    end
    
    Nps_buffer_orig = Nps_buffer;

    
    
    % ----------------------------------------------------------------------
    % Open files
    % ----------------------------------------------------------------------
    
    ps_fid = fopen([project_id '_ps_results_' char(results_id(w)) ...
                    '_sel'  num2str(z) '.raw'],'r');
    psc_valid_fid = fopen([project_id '_psc_valid_sel' num2str(z) '.raw'],'r');
    psp_valid_fid = fopen([project_id '_psp_valid.raw'],'r');
    ref_fid = fopen([project_id '_ref_sel' num2str(z) '.raw'],'r');
    ref_array = fread(ref_fid,[3 Nref(z)],'double')';
    fclose(ref_fid);

    psc_valid = fread(psc_valid_fid,[3*Nifgs+Npar_max+2 Npsc(z)],'double')';
    psp_valid = fread(psp_valid_fid,[3*Nifgs+Npar_max+2 Npsp(z)],'double')';
    ps3_valid = [psc_valid;psp_valid];
    
    psc_az = psc_valid(:,1);
    psc_r = psc_valid(:,2);
    psc_phase = psc_valid(:,3:Nifgs+2);
    psc_ambi = psc_valid(:,Nifgs+3:2*Nifgs+2);
    psc_s_atmo = psc_valid(:,2*Nifgs+3:3*Nifgs+2);
    psc_param = psc_valid(:,3*Nifgs+3:3*Nifgs+Npar_max+2);
    
    psp_az = psp_valid(:,1);
    psp_r = psp_valid(:,2);
    psp_phase = psp_valid(:,3:Nifgs+2);
    psp_ambi = psp_valid(:,Nifgs+3:2*Nifgs+2);
    psp_s_atmo = psp_valid(:,2*Nifgs+3:3*Nifgs+2);
    psp_param = psp_valid(:,3*Nifgs+3:3*Nifgs+Npar_max+2);
    
    ps3_az = ps3_valid(:,1);
    ps3_r = ps3_valid(:,2);
    ps3_phase = ps3_valid(:,3:Nifgs+2);
    ps3_ambi = ps3_valid(:,Nifgs+3:2*Nifgs+2);
    ps3_s_atmo = ps3_valid(:,2*Nifgs+3:3*Nifgs+2);
    ps3_param = ps3_valid(:,3*Nifgs+3:3*Nifgs+Npar_max+2);
    
    Nps_buffer = Nps_buffer_orig;
    
    for v = 1:Nbuffers
      
      if (v == Nbuffers)&&(Nps_rem_buffer~=0)
        Nps_buffer = Nps_rem_buffer;
      end
      
      % ----------------------------------------------------------------------
      % Read data
      % ----------------------------------------------------------------------
      
      ps_data = fread(ps_fid,[6+Npar_max+3*Nifgs Nps_buffer],'double')';
      
      ps_comment = ps_data(:,1);
      ps_az = ps_data(:,2);
      ps_r = ps_data(:,3);
      ps_param = ps_data(:,4:Npar_max+3);
      ps_ens_coh = ps_data(:,Npar_max+4);
      ps_ens_coh_local = ps_data(:,Npar_max+5);
      ps_sig2hat = ps_data(:,Npar_max+6);
      ps_phase_unw = ps_data(:,Npar_max+7:Nifgs+Npar_max+6);
      ps_defo_tot = ps_data(:,Nifgs+Npar_max+7:2*Nifgs+Npar_max+6);
      ps_atmo = ps_data(:,2*Nifgs+Npar_max+7:3*Nifgs+Npar_max+6);
      clear ps_data
      
      
      % ----------------------------------------------------------------------
      % Compare input and output
      % ----------------------------------------------------------------------
      
      [dummy,index1a,index1b] = intersect([psc_az psc_r],[ps_az ps_r],'rows');
      [dummy,index2a,index2b] = intersect([psp_az psp_r],[ps_az ps_r],'rows');
      [dummy,index3a,index3b] = intersect([ps3_az ps3_r],[ps_az ps_r],'rows');
      
      psc_phase_unw = 2*pi*psc_ambi+psc_phase;
      psp_phase_unw = 2*pi*psp_ambi+psp_phase;
      ps3_phase_unw = 2*pi*ps3_ambi+ps3_phase;
      psc_phase_unw_orig = psc_phase_unw;
      psc_phase_unw = psc_phase_unw - repmat(psc_phase_unw_orig(ref_array(1,1),:),Npsc(z),1);
      psp_phase_unw = psp_phase_unw - repmat(psc_phase_unw_orig(ref_array(1,1),:),Npsp(z),1);
      ps3_phase_unw = ps3_phase_unw - repmat(psc_phase_unw_orig(ref_array(1,1),:),Npsc(z)+Npsp(z),1);
    ps_atmo(find(ps_comment==0),:)
      ps_phase_unw_atmo = ps_phase_unw+ps_atmo-repmat(ps_atmo(find(ps_comment==0),:),Nps_buffer,1);
      
      figure;
      imagesc(psc_phase_unw(index1a,:)-ps_phase_unw(index1b,:));colorbar
      title('psc-ps');
      figure;
      imagesc(psc_phase_unw(index1a,:)-ps_phase_unw_atmo(index1b,:));colorbar
      title('psc-ps_atmo');
      figure;
      imagesc(psp_phase_unw(index2a,:)-ps_phase_unw(index2b,:));colorbar
      title('psp-ps');
      figure;
      imagesc(psp_phase_unw(index2a,:)-ps_phase_unw_atmo(index2b,:));colorbar
      title('psp-ps_atmo');
      figure;
      imagesc(ps3_phase_unw(index3a,:)-ps_phase_unw(index3b,:));colorbar
      title('ps3-ps');
      figure;
      imagesc(ps3_phase_unw(index3a,:)-ps_phase_unw_atmo(index3b,:));colorbar
      title('ps3-ps_atmo');

      figure;
      scatter(ps3_r(index3a),ps3_az(index3a),5,1000*(ps3_param(index3a,2)-psc_param(ref_array(1,1),2)));colorbar
      figure;
      scatter(ps_r(index3b),ps_az(index3b),5,1000*ps_param(index3b,2));colorbar

      figure;
      scatter(ps_r(index3b),ps_az(index3b),5,1000*(ps3_param(index3a,2)-psc_param(ref_array(1,1),2)-ps_param(index3b,2)));colorbar

      figure;
      scatter(ps_r(index3b),ps_az(index3b),5,1000*mean(ps3_s_atmo(index3a,:),2));colorbar

      for w = 1:Nifgs
      figure;
      scatter(ps3_r(index3a),ps3_az(index3a),5,1000*(ps3_param(index3a,2)-psc_param(ref_array(1,1),2)-ps3_s_atmo(index3a,w)-psc_s_atmo(ref_array(1,1),w)));colorbar
      figure;
      scatter(ps_r(index3b),ps_az(index3b),5,1000*ps_atmo(index3b,w)/m2ph);colorbar
      figure;
      scatter(ps_r(index3b),ps_az(index3b),5,1000*(ps3_param(index3a,2)-psc_param(ref_array(1,1),2)-ps3_s_atmo(index3a,w)-psc_s_atmo(ref_array(1,1),w))-1000*ps_atmo(index3b,w)/m2ph);colorbar
      display('pause')
      pause
      end
      
      %figure;
      %imagesc(psc_phase_unw(index1a,:));colorbar
      %title('psc');
      %figure;
      %imagesc(ps_phase_unw(index1b,:));colorbar
      %title('ps');
      %figure;
      %imagesc(ps_phase_unw_atmo(index1b,:));colorbar
      %title('ps_atmo');
      
      %figure;
      %imagesc(psp_phase_unw(index2a,:));colorbar
      %title('psc');
      %figure;
      %imagesc(ps_phase_unw(index2b,:));colorbar
      %title('ps');
      %figure;
      %imagesc(ps_phase_unw_atmo(index2b,:));colorbar
      %title('ps_atmo');
      
      %figure;
      %imagesc(ps3_phase_unw(index3a,:));colorbar
      %title('ps3');
      %figure;
      %imagesc(ps_phase_unw(index3b,:));colorbar
      %title('ps');
      %figure;
      %imagesc(ps_phase_unw_atmo(index3b,:));colorbar
      %title('ps_atmo');
      
      error('dfdf')
    end
  end
end

