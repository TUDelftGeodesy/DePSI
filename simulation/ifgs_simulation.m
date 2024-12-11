function [topo,defo1,defo2,defo3,defo_quad,defo_cubic,defo_nl1,defo_nl2,subpixel,atmo,ifgs,ifgs_unw,ambi,breakpoint,breakpoint2,noise_level,t,a,Bdop] = ifgs_simulation(M,Nifgs,az_spacing,sat_vel,model,method,std_H,std_D,std_atmo,power,regime)

%[topo,defo1,defo2,defo_nl1,defo_nl2,subpixel,atmo,ifgs,ambi,breakpoint,noise_level,t,a,Bdop] = ifgs_simulation(M,Nifgs,az_spacing,sat_vel,model,method,std_H,std_D,std_atmo,power,regime)
%
% Simulate Nifgs SAR interferograms sampled with PSs
% 
% INPUT
%
% M           size of the images
% Nifgs       number of interferograms
% az_spacing  azimuth spacing
% sat_vel     velocity of satellite
% model       1 DEM error + lin_defo + nl_defo + atmo + orbit + noise
%             2 DEM error + lin_defo + nl_defo + atmo + noise
%             3 DEM error + lin_defo + atmo + noise
%             4 DEM error + lin_defo + nl_defo
%             5 DEM error + lin_defo + atmo
%             6 DEM error + lin_defo + noise
%             7 DEM error + lin_defo 
%             8 DEM error
%             9 DEM error + breakpoint
%             10 DEM error + breakpoint + nl_defo
%             11 noise
%             12 atmo + noise
%             13 atmo
%             14 DEM error + lin_defo + nl_defo + atmo
%             15 DEM error + lin_defo + atmo + subpixel
%             16 DEM error + lin_defo + subpixel
%             (more options can easily be specified)
% method      simulation method, 'frac' for fractals,
%             'shape' for pre-defined shaped, e.g. peaks, cone,
%             bowl
% std_H       standard deviation of simulated heights
% std_D       standard deviation of simulated deformation
% power       [exp1 exp2 exp2]  
%             power law in the atmo delay simulation as isotropic 2D fractal surface 
%             default value: [5/3 8/3 2/3]
% regime      [p1 p2 p3] 
%             cumulative percentage of spectrum covered by a specific beta
%             defaule value: [95 99 100]
%        
%
% OUTPUT
%
% topo        simulated topography [m]
% defo1        simulated velocity field1
% defo2        simulated velocity field2
% defo_nl1      amplitude of periodic (= non-linear) deformation
% defo_nl2          periodic offset of non-linear deformation
% atmo        simulated atmospheric delay for Nifgs interferograms
% ifgs        simulated interferograms
% ifgs_unw    simulated unwrapped interferograms
% ambi        simulated ambiguities
% breakpoint  breakpoint in deformation
% noise_level level of simulated noise
% t           temporal baselines
% a           perpendicular baselines
%
% Example:
% [topo,vel,defo_nl1,defo_nl2,acq,totatmo,totint,table,rel] = simulation_sim(128,1,'c',5,1,10);

% set default values for optional input parameter
if nargin==11 regime=[95 99 100]; end;
if nargin==10 power=[5/3 8/3 2/3]; regime=[95 99 100]; end;

% constant values
theta = 23*pi/180;     % incident angle [deg]
lambda = 0.0565646;    % wavelength [m]
H = 780000;            % satellite vertical height [m]
R = H/cos(theta);      % antenna-target distance [m]
m2ph = -4*pi/lambda;
rand('state',sum(100*clock));
randn('state',sum(100*clock));

switch model
  
 case 1
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 1;
  breakpoint_flag = 0;
  def_nl_flag = 1;
  atmo_flag = 1;
  orbit_flag = 1;
  noise_flag = 1;
  subpixel_flag = 0;
  
 case 2
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 1;
  breakpoint_flag = 0;
  def_nl_flag = 1;
  atmo_flag = 1;
  orbit_flag = 0;
  noise_flag = 1;
  subpixel_flag = 0;
  
 case 3
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 1;
  breakpoint_flag = 0;
  def_nl_flag = 0;
  atmo_flag = 1;
  orbit_flag = 0;
  noise_flag = 1;
  subpixel_flag = 0;
 
 case 4
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 1;
  breakpoint_flag = 0;
  def_nl_flag = 1;
  atmo_flag = 0;
  orbit_flag = 0;
  noise_flag = 0;
  subpixel_flag = 0;
 
 case 5
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 1;
  breakpoint_flag = 0;
  def_nl_flag = 0;
  atmo_flag = 1;
  orbit_flag = 0;
  noise_flag = 0;
  subpixel_flag = 0;

 case 6
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 1;
  breakpoint_flag = 0;
  def_nl_flag = 0;
  atmo_flag = 0;
  orbit_flag = 0;
  noise_flag = 1;
  subpixel_flag = 0;
 
 case 7
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 1;
  breakpoint_flag = 0;
  def_nl_flag = 0;
  atmo_flag = 0;
  orbit_flag = 0;
  noise_flag = 0;
  subpixel_flag = 0;
 
 case 8
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 0;
  breakpoint_flag = 0;
  def_nl_flag = 0;
  atmo_flag = 0;
  orbit_flag = 0;
  noise_flag = 0;
  subpixel_flag = 0;
 
 case 9
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 0;
  breakpoint_flag = 1;
  def_nl_flag = 0;
  atmo_flag = 0;
  orbit_flag = 0;
  noise_flag = 0;
  subpixel_flag = 0;
 
 case 10
  topo_flag = 0;
  DEM_error_flag = 1;
  def_lin_flag = 0;
  breakpoint_flag = 1;
  def_nl_flag = 1;
  atmo_flag = 0;
  orbit_flag = 0;
  noise_flag = 0;
  subpixel_flag = 0;
 
 case 11
  topo_flag = 0;
  DEM_error_flag = 0;
  def_lin_flag = 0;
  breakpoint_flag = 0;
  def_nl_flag = 0;
  atmo_flag = 0;
  orbit_flag = 0;
  noise_flag = 1;
  subpixel_flag = 0;
 
 case 12
  topo_flag = 0;
  DEM_error_flag = 0;
  def_lin_flag = 0;
  breakpoint_flag = 0;
  def_nl_flag = 0;
  atmo_flag = 1;
  orbit_flag = 0;
  noise_flag = 1;
  subpixel_flag = 0;
 
 case 13
  topo_flag = 0;
  DEM_error_flag = 0;
  def_lin_flag = 0;
  breakpoint_flag = 0;
  def_nl_flag = 0;
  atmo_flag = 1;
  orbit_flag = 0;
  noise_flag = 0;
  subpixel_flag = 0;
  
  case 14
    topo_flag = 0;
    DEM_error_flag = 1;
    def_lin_flag = 1;
    breakpoint_flag = 0;
    def_nl_flag = 1;
    atmo_flag = 1;
    orbit_flag = 0;
    noise_flag = 0;
    subpixel_flag = 0;
    
  case 15
    topo_flag = 0;
    DEM_error_flag = 1;
    def_lin_flag = 1;
    breakpoint_flag = 0;
    def_nl_flag = 0;
    atmo_flag = 1;
    orbit_flag = 0;
    noise_flag = 0;
    subpixel_flag = 1;
  
  case 16
    topo_flag = 0;
    DEM_error_flag = 1;
    def_lin_flag = 1;
    breakpoint_flag = 0;
    def_nl_flag = 0;
    atmo_flag = 0;
    orbit_flag = 0;
    noise_flag = 0;
    subpixel_flag = 1;
    
  otherwise
    error(['Error: please specify a possible ''model''']);
end; % switch

% generate orthogonal baseline values [Bn err errBn] and times
fprintf(1,'Acquisition simulation...\n');
fprintf(1,'\n');
[a b c] = baseline(Nifgs);
acqdates = ERS_time_series(2,Nifgs);
tt =  acqdates(:,1); % [days]
t  = ((tt - min(tt))/365) - 5;    % [year], zero somewhere in the
                                  % middle
Bdop = randn(Nifgs,1)*40;


if (breakpoint_flag == 1)
  breakpoint = ceil(rand(1)*(Nifgs-2))+1; % to make sure 2<=b<=Nifgs-1
  if t(breakpoint)> 0
    t1 = [t(1:breakpoint);repmat(t(breakpoint),Nifgs-breakpoint,1)];
    t2 = [zeros(breakpoint,1);t(breakpoint+1:end)- ...
        repmat(t(breakpoint),Nifgs-breakpoint,1)];
  elseif t(breakpoint)< 0
    t1 = [t(1:breakpoint-1)-repmat(t(breakpoint),breakpoint-1,1);zeros(Nifgs-breakpoint+1,1)];
    t2 = [repmat(t(breakpoint),breakpoint-1,1);t(breakpoint:end)];
  else
    error(['Somehow you managed to make a master-master interferogram????']);
  end
else
  breakpoint = [];
  t1 = t;
  t2 = zeros(Nifgs,1);
end

%TODO
breakpoint2 = []; %create this option later
defo3 = zeros(M);
defo_quad = zeros(M);
defo_cubic = zeros(M);

acq = [a b t];
fprintf(1,'\n');


% big switch for the simulation products

if (topo_flag == 1) & (DEM_error_flag == 1)
  error('Error: you can only simulate topography or DEM error');
elseif (topo_flag == 1)
  % generate topography
  fprintf(1,'Topography simulation...');
  switch method
   case 'frac'
     topo = simtopo(M);      
   case 'shape'
     topo = cone(M)*std_H;
    otherwise
     error('You specified a wrong ''method''.');
  end
  fprintf(1,'DONE\n');
elseif (DEM_error_flag == 1)
  % generate DEM error
  fprintf(1,'DEM error simulation...');
  topo = 0.5*randn(M)*std_H;
  fprintf(1,'DONE\n');
else
  topo = zeros(M);
end

if (def_lin_flag == 1)|(breakpoint_flag == 1)
  % generate velocity field
  fprintf(1,'Linear deformation simulation...');
  switch method
   case 'frac'
     defo1 = fracsurf(M,3,'n');
   case 'shape'
     %defo1 = peaks(M)/1000;
     %defo1 = peaks(M)/1000+std_D*bowl2(M);
     defo1 = std_D*bowl2(M);
   otherwise
     error('You specified a wrong ''method''.');
  end
  
  if (breakpoint_flag == 1)
    switch method
     case 'frac'
       defo2 = fracsurf(M,3,'n');
     case 'shape'
       defo2 = rot90(peaks(M)/1000,-1); %rotate peaks 90deg
     otherwise
       error('You specified a wrong ''method''.');
    end
  else
    defo2 = zeros(M);
  end
  fprintf(1,'DONE\n');
else
  defo1 = zeros(M);
  defo2 = zeros(M);
end

if (def_nl_flag == 1)
  % generate amplitude and defo_nl2 for periodic deformation
  fprintf(1,'Periodic deformation simulation...');
  defo_nl1 = fracsurf(M,3,'n');
  defo_nl2 = mod(100*fracsurf(M,3,'n')+0.5,1)-0.5;
  defo_nl2(defo_nl1<0) = defo_nl2(defo_nl1<0)+0.5;
  defo_nl2 = mod(defo_nl2+0.5,1)-0.5;
  defo_nl1 = abs(defo_nl1); 
  %defo_nl2 = 10*fracsurf(M,3,'n'); 
  % factor 100 causes ~ -0.5<defo_nl2<0.5 
  fprintf(1,'DONE\n');
else
  defo_nl1 = zeros(M);
  defo_nl2 = zeros(M);
end

if (atmo_flag == 1)
  % create master atmosphere and Nifgs slave atmospheres with std of 1.7
  % mm (P0 = 1) (P0 = 5, std = 8 mm)
  fprintf(1,'Atmospheric delay simulation...');
  %std_atmo_vec = randn(Nifgs+1,1)*std_atmo;
  std_atmo_vec = rand(Nifgs+1,1)*std_atmo;
  atmo = NaN([M M Nifgs+1]); % first is master 
  for v = 1:Nifgs+1
    atmo_temp = simatmo(M, [5/3 8/3 2/3], [10  90 100],1,'n');
    std_temp = std(atmo_temp(:));
    atmo(:,:,v) = atmo_temp*std_atmo_vec(v)/std_temp;
  end; 
  fprintf(1,'DONE\n');         
else
  atmo = zeros(M,M,Nifgs+1);
end
%atmo(:,:,1) = 0;


if (subpixel_flag == 1)
  fprintf(1,'Subpixel position simulation...');
    subpixel = rand(M,M)*az_spacing;
else
    subpixel = zeros(M,M);
end

if (noise_flag == 1)
  fprintf(1,'Noise simulation...');
  noise_std = NaN([M M Nifgs+1]);
  noise_level = NaN(Nifgs+1,1);

  noise_level(1) = pi*15/180;
  %noise_level(1) = pi*1/180;
  noise_std(:,:,1) = randn(M)*noise_level(1);
  for v = 1:Nifgs
    noise_level(v+1) = (pi*20/180)+randn(1)*(pi*5/180);
    %noise_level(v+1) = (pi*2/180)+randn(1)*(pi*1/180);
    noise_std(:,:,v+1) = randn(M)*noise_level(v+1);
  end
  fprintf(1,'DONE\n');         
else
  noise_std = zeros(M,M,Nifgs+1);
  noise_level = zeros(Nifgs+1,1);
end


fprintf(1,'Generation of interferograms...\n')

ifgs = NaN([M M Nifgs]);
ifgs_unw = NaN([M M Nifgs]);

for v = 1:Nifgs
    
    topo_ph = m2ph/(sin(theta)*R)*topo*a(v); % topographic phase term

    if (orbit_flag == 1)
      orbit_ph  = m2ph/(sin(theta)*R)*topo*b(v);  % residual fringes (error in orbit determination)
    else
      orbit_ph = zeros(M);
    end
 
    defo_ph = m2ph*(defo1*t1(v)+defo2*t2(v));      % deformation phase term

    defo_nl_ph = m2ph*(defo_nl1.*sin(2*pi*(t(v)-defo_nl2)) + ...
                       defo_nl1.*sin(2*pi*defo_nl2));% non-linear deformation
    
    atmo_ph = m2ph*(atmo(:,:,1) - atmo(:,:,v+1)); % atmospheric
                                                  % phase term
    subpixel_ph = (2*pi/sat_vel)*Bdop(v)*subpixel;
    
    noise_ph = noise_std(:,:,1) + noise_std(:,:,v+1);

    ifgs_unw(:,:,v) = (topo_ph + orbit_ph + defo_ph + defo_nl_ph + ...
                       atmo_ph + subpixel_ph + noise_ph);
    
    ifgs(:,:,v)=(mod(ifgs_unw(:,:,v)+pi,2*pi)-pi);

    dummy = [v Nifgs];
    fprintf(1,'Interferogram %3.0f of %3.0f\n',dummy);

end; %for (v)

ambi = round((ifgs_unw-ifgs)/(2*pi));

