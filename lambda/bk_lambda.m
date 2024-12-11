%%% script to test my lambda estimator
%%% 1) simulate data, lin.defo,topo, bias  + sinus!
%%% 2) setup system of equations as y=Aa+Bb, A are ambiguities
%%% 3) add 3 pseudo observations h=0,v=0,b=0
%%% 4) LS exact float solution for a,b
%%% 5) decorrelate Qa?
%%% 6) search/integer solution for ambiguities
%%%
%%% BK 14-Jan-2002
%%% clear all
clear all
close all

more off

N          = 60
noiselevel = 40
noiselevel = deg2rad(noiselevel);
noise      = noiselevel.*randn(N,1);
btemp      = 4*randn(N,1);% wrt some master at btemp=0?
bperp      = 400*randn(N,1);% wrt some master at btemp=0?

%%% functional model
lambda = 0.056;
rsintheta = 850000*sin(deg2rad(23));
KK     = -4*pi/lambda;
h2p    = KK*bperp./rsintheta;
v2p    = KK*btemp*1e-3;
v2p2   = KK*sin(btemp)*1e-3;
b2p    = ones(N,1);

%%% signal
topo   =  20.1;
defo   =  15.6;
defo2  =  5.6;
bias   =  deg2rad(60.0);
%%%
phase  = h2p*topo + v2p*defo + v2p2*defo2 + b2p*bias + noise;
y      = wrap(phase);
%Qy     = deg2rad(40)*eye(N);% simple model for now
A      = 2*pi*eye(N);
B      = [h2p,v2p,v2p2,b2p];
NP     = 4;
NP     = size(B,2);

%%% true number of input wraps...
true_wraps = -round(phase./(2*pi));
disp('signal simulated');

%%% pseudo-observations for float solution and vc-matrix
y2           = [y; zeros(NP,1)];
A2           = [A; zeros(NP,N)];
B2           = [B; eye(NP)];
%Qy2 = Qy
%sigma_phase  = deg2rad(noiselevel);% on diff. see what's the effect
sigma_phase  = deg2rad(50);% on diff. see what's the effect, 50 seems OK, if noise=40
sigma_topo   = 25;% [m]
sigma_defo   = 15;% [mm/y]
sigma_defo2  = 10;% [mm/y]
sigma_bias   = deg2rad(40);% [rad]
var_phase    = sigma_phase^2;
var_topo     = sigma_topo^2;
var_defo     = sigma_defo^2;
var_defo2    = sigma_defo2^2;
var_bias     = sigma_bias^2;
Qy           = var_phase*eye(N);% simple model for now
Qy2          = var_phase*eye(N+NP);% simple model for now
%for i=1:NP
Qy2(N+1,N+1) = var_topo;% [m]
Qy2(N+2,N+2) = var_defo;% [mm/y]
Qy2(N+3,N+3) = var_defo2;% [mm/y]
Qy2(N+NP,N+NP) = var_bias;% [rad]


%%% float solution for ambiguities, partioned model
p_a2   = eye(N+NP) - B2 * inv(B2.'*inv(Qy2)*B2) * B2.'*inv(Qy2);
A2_    = p_a2*A2;
Qahat  = inv(A2_.'*inv(Qy2)*A2_);% is this correct???
ahat   = Qahat*A2_.'*inv(Qy2)*y2;% float solution partioned model
disp('float solution computed');

%%% check if xhat is correct partioned computed
%AA    = [A2,B2];
%ahat2 = inv(AA.'*inv(Qy2)*AA)*AA.'*inv(Qy2)*y2;
%disp('solution for topo etc.');
%disp(ahat2(N+1:N+NP));
%%% it is, but of course, the computed ambiguities are exactly through points
%%% and have nothing to do with the true ambiguities...
%%% one may want to skip the first smallest expected ambiguities?


%%% lambda method
afloat = ahat;
ncands = 1;
n      = size (Qahat,1);
afixed = zeros(n,ncands);
sqnorm = zeros(1,ncands);
% ----------------------------------------------------------------------
incr   = afloat - rem(afloat,1);
afloat = rem(afloat,1);
% ----------------------------------------------------------------------
% Compute the Z-transformation based on L and D of Q, ambiguities
% are transformed according to \hat{z} = Z^T\hat{a}, Q is transformed
% according to Q_\hat{z} = Z^T * Q_\hat{a} * Z
% ----------------------------------------------------------------------
[Qz,Z,L,D,zfloat] = bk_decorrel(Qahat,afloat);
% RH Berekenen successrate:
%   successrate    = prod (2 * normcdf (1./(2.*sqrt(D)),0,1) - 1);
%   fprintf(1,'Successrate : %6.3f\n',successrate);
%PJ: Let op: Deze successrate is de "werkelijke" successrate bij gebruik van
%de bootstrap estimator, het is dus een ONDERGRENS voor de lambda methode
%(omdat die LAMBDA-methode beter werkt (of even goed) als de
%bootstrap-methode. Na decorrelatie is het wel een scherpe ondergrens.
% ----------------------------------------------------------------------
% Compute a suitable Chi^2 such that we have the requested number of 
% candidates at minimum; use an 'eps' to make sure the candidates are 
% inside the ellipsoid and not exactly on the border.
% ----------------------------------------------------------------------

%+++++++++++++++++++
%%% FOR A NUMBER OF TESTS< SEE IF WE CAN RECOVER INPUT
ntries    = 1000
ntries    = 40
if (ntries >= 10) profile on; end;
qfactor = 1.0
n_notfound = 0;
n_skipped  = 0;
ttot      = 0;
TBK      = 0;
TRH      = 0;
qtt = [];
all_topo  = 20.*randn(ntries,1);% more than specified in psuedo obs.
all_defo  = 25.*randn(ntries,1);
all_defo2 = 15.*randn(ntries,1);
all_bias = deg2rad(40.*randn(ntries,1));
invZ = inv(Z);
invQy = inv(Qy);
invQy2 = inv(Qy2);
invL = inv(L);% offer to new bk_lsearch
invD = 1 ./D;% offer to new bk_lsearch
LS_projector = Qahat*A2_.'*invQy2;
for t=1:ntries
  disp(t);
  %%% signal
  topo   = all_topo(t);
  defo   = all_defo(t);
  defo2  = all_defo2(t);
  bias   = all_bias(t);
  noise  = noiselevel.*randn(N,1);
  input  = [topo,defo,defo2,bias]
  %%%
  phase  = h2p*topo + v2p*defo + v2p2*defo2 + b2p*bias + noise;
  y      = wrap(phase);
  y2     = [y; zeros(NP,1)];
t0=cputime;
afloat   = LS_projector*y2;% float solution partioned model
%% Make estimates in 'a' between -1 and +1 by subtracting an
%%% BK: always for insar?
if (max(abs(afloat)) > 1) 
  warning('something is wrong, afloat not in -1,1??');
end
%incr   = afloat - rem(afloat,1);
%afloat = rem(afloat,1);
zfloat = Z.' * afloat;%
%%% does Chi2 change? yes: diff with chistart and new (RH)
Chi2 = bk_chistart(D,L,zfloat,ncands,qfactor);
Chi2 = Chi2
crit_skip = N;%%% when to look for other initial guess...
if (Chi2 > crit_skip) 
%if (Chi2 > N*2.0) 
  warning('chi2>N*1.5, may take rel. long..., reducing it');
  warning('looking for other initial values before skipping this one...');
  ierr = 1;
  zfixed = zeros(n,1);
  %%% skip this...
  n_skipped = n_skipped+1;
  %%% TEST IF OTHER PSEUDO OBS> WOULD HELP
  %%% or connect close by same amb. (is already so via Qy??)
y3 = y2;
%%% try a couple of times (does this do any good???)...
for hmm=1:10
  y3(n+1) = 15*randn;% topo
  y3(n+2) = 10*randn;% defo1
  y3(n+3) = 10*randn;% defo2
  y3(n+4) = deg2rad(20*randn);% bias
  y3(n+1:n+4).'
  ahat2   = LS_projector*y3;% float solution partioned model
  %if (max(abs(ahat2)) > 1) 
    %warning('something is wrong, afloat not in -1,1??');
  incr   = ahat2 - rem(ahat2,1);% seems to happen if not 0 pseudo obs.
  ahat2 = rem(ahat2,1);
  %end
  zfloat2 = Z.' * ahat2;%
  %zfloat.' - zfloat2.'%%%see if it changed...
  % see if Chi3 < Chi2
  Chi3 = bk_chistart(D,L,zfloat2,ncands,qfactor)
  %%% use this one... is faster... but not realistic?
  %%% assume we can always decrease and still find correct max?
  %%% by always doing this..., and take the min. of 10 tries
  %%% requires speedup of chistart
  if (Chi3 < crit_skip)
	  n_skipped=n_skipped -1;
    [zfixed,sqnorm,ierr] = bk_lsearch(zfloat2,invL,invD,Chi3,ncands);
    afixed = (zfixed' * inv(Z))';
    afixed = round(afixed + repmat(incr,1,ncands)); % replace with increments
    break;
  end
end
%%% if Chi3 is much smaller than Chi2, use it...
%keyboard
else

%t0=cputime;
%Chi2 = 80
%t0=cputime;
[zfixed,sqnorm,ierr] = bk_lsearch (zfloat,invL,invD,Chi2,ncands);
%TBK = TBK + cputime-t0;
%t0=cputime;
%[zfixed2,sqnorm2,ierr2] = lsearch (zfloat,L,D,Chi2,ncands);
%TRH = TRH + cputime-t0;
%qqqq= max(abs(zfixed-zfixed2))
% ----------------------------------------------------------------------
% --- Perform the back-transformation and add the increments
% ----------------------------------------------------------------------
%afixed = (afixed' * inv(Z))'; %RH Jul 11 15:43:09 CEST 2001
%afixed = (zfixed' * inv(Z))';
afixed = round(zfixed.' * invZ).';
end%%%skip not likely one
% ----------------------------------------------------------------------
% End of routine: lambda2
% ----------------------------------------------------------------------
if ierr == 1
  n_notfound = n_notfound+1;
  %zfixed
  fprintf(1,'Not enough candidates were found!!\n'); %RH 21 aug
end

ttot = ttot+ (cputime-t0);


%%% compute xhat for fixed solution
%disp('computing final parameter solution');
phase_uw = y-A*afixed;
%xhat     = LS_projector*phase_uw;
xhat     = inv(B.'*invQy*B)*B.'*invQy*phase_uw;

xhat  = xhat.';


%%% check results.
true_wraps = -round(phase./(2*pi));% +-?
truth = true_wraps.';
fixed = afixed.';
truth_min_fixed = max(abs(truth-fixed));
if (truth_min_fixed ~= 0) 
  warning('wrong solution obtaine...');
  xhat=xhat
end
qtt = [qtt,truth_min_fixed];
%%% store results

end%FOR

%+++++++++++++++++++
q=find(qtt==0);
n_notfound = n_notfound
n_skipped  = n_skipped
successrate = 100*(length(q))/(ntries-n_notfound)

time_per = ttot/ntries
if (ntries >= 10) profile plot; profile off; end;

%%% EOF

