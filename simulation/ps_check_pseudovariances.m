clear all
close all

global fig ps_method ps_eval_method project_id
global m2ph Npar_max
global max_mem_buffer orbit Nlines Npixels az_spacing r_spacing
global Npsc_selections alpha0 gamma0 weighting
global visible_plots detail_plots std_noise std_atmo std_defo
global atmo_range defo_range sat_vel

project_id = 'test_v1.7.1';

load([project_id '_project.mat']);
ps_method = 'boot';
psc_model = 8;

array_H = NaN([NstdH NstdA NstdS NstdD]);
array_A = NaN([NstdH NstdA NstdS NstdD]);
array_S = NaN([NstdH NstdA NstdS NstdD]);
array_D = NaN([NstdH NstdA NstdS NstdD]);
array_ambi = NaN([NstdH NstdA NstdS NstdD]);

for z = 1:2
  dpsc_valid_fid = fopen([project_id '_dpsc_valid' num2str(z) '.raw'],'r');
  dpsc_valid = fread(dpsc_valid_fid,[3*Nifgs+Npar_max+2 Narcs(z)],'double')';
  fclose(dpsc_valid_fid);
  dpsc_ambi_valid = dpsc_valid(:,3+Nifgs:2*Nifgs+2);
  dpsc_param_valid = dpsc_valid(:,3*Nifgs+3:3*Nifgs+Npar_max+2);
  for a = 1:NstdH
    for b = 1:NstdA
        for c = 1:NstdS
            for d = 1:NstdD
                load(['psc_results_z_' num2str(z) '_method_' ps_method ...
                '_model_' num2str(psc_model) '_' num2str(a) num2str(b) ...
               num2str(c) num2str(d) '.mat']);
               dpsc_param = dpsc_data(:,1:Npar_max);
               dpsc_ambi = dpsc_data(:,Npar_max+2:Npar_max+Nifgs+1);
               
               array_H(a,b,c,d) = sum(abs(dpsc_param(:,1)-dpsc_param_valid(:,1)));
               array_A(a,b,c,d) = sum(abs(dpsc_param(:,2)-dpsc_param_valid(:,2)));
               array_S(a,b,c,d) = sum(abs(dpsc_param(:,3)-dpsc_param_valid(:,3)));
               array_D(a,b,c,d) = sum(abs(dpsc_param(:,4)-dpsc_param_valid(:,4)));
               array_ambi(a,b,c,d) = sum(abs(dpsc_ambi(:)-dpsc_ambi_valid(:)));
            
            end
        end
    end
end
end