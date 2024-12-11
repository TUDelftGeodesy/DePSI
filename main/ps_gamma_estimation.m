function [estim_gamma phase_residual] = ps_gamma_estimation(cpx_ifgs,dist,h2ph,gammas,plot_label);

% Function estim the gamma ,coherence as defined in A. Hooper
% Thesis, with two differences:
% the distances are taken into account when calculating the 
% mean values of the complex numbers  and subpixel position is calculated. 
% 
% Input:    - cpx_ifgs        interferometric complex numbers of
%                             size nof_pixel x Nof_igfs
%           - dist            distances between the points, 1 
%                             column arranged as the rows in ifgs_ph
%           - h2ph            height to phase factors
%           - gammas          estimated gammas used for weigthing,
%                             if non gammas is assigned to 1
%
% Output:   - estim_gamma     estimated coherence (Gamma)
%           - dem_err         estimated dem error
%				 
% ----------------------------------------------------------------------
% File............: ps_gamma_estimation.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Authors.........: Miguel Caro Cuenca
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

global m2ph sat_vel

%m2ph = -4*pi/lambda;
 
filter_size=800;%atm decorrelation length
%time_size=abs( corr(Btemp,h2ph')*10);
max_dem_error=30;

 
[nof_ps  Nifgs] = size( cpx_ifgs );

%n_trial_wraps=round(Nifgs/7.5);%this is used to set the search spce for the DEM-error estmations

n_trial_wraps=abs(round(max_dem_error*max(abs(h2ph))*m2ph/2/pi));

if n_trial_wraps<4
n_trial_wraps=4;
end

if nargin <5 
  plot_label = [];
%  gammas=repmat(1,nof_ps,1) ;

end



%if nargin <5 
%  gammas=repmat(1,nof_ps,1) ;

%end


ind_remove = find(gammas < 0.2);

ind_center = find(dist<0.0001);

%keyboard
if ~isempty(intersect (ind_remove, ind_center))

	estim_gamma =0;
	phase_residual = zeros(1,Nifgs);

else

	gammas(ind_remove) =0;
	gammas(ind_center)=0;%remove the PS in question from the averaging
%	gammas(ind_center)=1;

	%GAUSSIAN window seems more appropiate since in the spcetral domain
	%there is not side lobes, the fourier transform of a gaussian is another gaussian
	%However, since the domain is finite, in rality what we do is a gaussian times a rect, producing
	%the sidelobes in the spectral domain.  
	%inv_weight= sqrt(dist+min(dist(ind_rest))/2);%Weighting with the inverse of the square root of the distance 
	%weight = ( gammas ./inv_weight)/sum( gammas./inv_weight);
	gaus_win = exp(dist.^2/(-2*filter_size^2));
	weight = (gammas .*gaus_win/sum( gammas.*gaus_win));


	filter_cpx_ph=NaN([Nifgs 1]);

%keyboard
	for n_ifg = 1:Nifgs
  	%Weighted Average
  	filter_cpx_ph(n_ifg,1) = transp( cpx_ifgs(:,n_ifg)) *weight ;
	end

	%figure;plot(Btemp,angle(filter_cpx_ph),'o');grid on


	%subtracting the filter_cpx_ph which should be atmo and defo from the original phases
	res_cpx_ifgs  =transp( cpx_ifgs(ind_center,:)) .* conj(filter_cpx_ph);


	%plot_label='plots/dem'
	[dem_err_phase, C0, estim_gamma,  phase_residual] = ps_topofit(res_cpx_ifgs,h2ph',n_trial_wraps,plot_label);

	phase_residual =  cpx_ifgs(ind_center,:).*exp(-j*dem_err_phase*h2ph);
	%figure;plot(h2ph,angle(smoothed_ts.*exp(-j*dem_err_phase*h2ph)),'o');grid on

end %~isempty(intersect (ind_remove, ind_center))






