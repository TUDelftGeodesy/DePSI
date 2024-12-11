%Input_file template for ps_analysis.m
%(See the function ps_readinput_parameters.m for more information)

% General parameters
% ----------------------------------------------------------------------

max_mem_buffer = 50e7
visible_plots = 'n'
detail_plots = 'n'
processing_groups = []
run_mode = 'normal'


% Project parameters
% ----------------------------------------------------------------------

project_id = 'marken_rsat2'
input_file = []
processDir = '../../stacks/nl_noordholland_rsat2_asc_t109_xf2/process'
startDate = ''
stopDate = ''
excludeDate = ['20200509';'20200623']
ifgsVersion = '_srd'
altimg = ''
master = '20200923'
swath_burst = [2,504]
sensor = 's1'
orbit = 'asc'
processor = 'doris_rippl'
project = 'nobv_zegveld'
crop = []
az_spacing = 1.24
r_spacing = 1.33
slc_selection_input = []
ifg_selection_input = []
ref_cn = []
Ncv = 25
ps_method = 'perio'
psc_model = [1]
ps_model = [1]
final_model = [2]
breakpoint = []
breakpoint2 = []
ens_coh_threshold = 0.5
varfac_threshold = 3
detrend_method = 'yes'
output_format = 1
stc_min_max = [30,100]
do_apriori_sidelobe_mask = 'yes'
do_aposteriori_sidelobe_mask = 'no'


% Geocoding parameters
%----------------------------------------------------------------------

master_res = 'slave.res'
ref_height = 0
demFile = 'dem_radar.raw'

% Psc parameters
%----------------------------------------------------------------------

amplitude_calibration = 'yes'
psc_selection_method = 'threshold'
psc_selection_gridsize = 300
psc_threshold = 0.2
max_arc_length = 5000
network_method = 'spider'
Ncon = 16
Nparts = 8
Npsc_selections = 1
filename_water_mask = []
gamma_threshold = 0.45
psc_distribution = 'uniform'
weighted_unwrap = 'yes'

% threshold is percentage of slc's that has an amplitude peak
livetime_threshold = 0.2;
% include local near maxima 
peak_tolerance = 0.9;


% Ps parameters
% ----------------------------------------------------------------------

psp_selection_method = 'ampdisp'
psp_threshold1 = 0.4
psp_threshold2 = []
ps_eval_method = 'psp'
Namp_disp_bins = 100
Ndens_iterations = 5
densification_flag = 'yes'
ps_area_of_interest = 'markenRsat2Asc_dike_mask.raw'
dsFileName = []
dens_method = 'orig'
dens_check = 'nocheck'
Nest = 1;


% Stochastic model parameters
% ----------------------------------------------------------------------

std_param = [30,0.005,1,0.02,0.01,0.01,0.005];
defo_range = 5000
weighting = 'vce'
ts_atmo_filter = 'gaussian'
ts_atmo_filter_length = 12/12
ts_noise_filter = 'gaussian'
ts_noise_filter_length = 8/12


% Bowl parameters
%-----------------------------------------------------------------------

defo_method = []
xc0 = []
yc0 = []
zc0 = []
r0 = []
r10 = []
epoch = []
