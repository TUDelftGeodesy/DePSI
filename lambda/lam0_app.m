function [lam0ap]= lam0_app(lambda,df,gam0,cv)
%LAM0_APP: approximation of the non-centrality parameter lambda
%
% This routines computes the first approximation of the non-centrality
% parameter lambda. It is based on the FORTRAN routine CHILNC.F by F.Kenselaar
% See Abramowitz/Stegun formula 26.4.28. This routine uses the normal
% distribution.
%
% Input:
%  - lambda : approximate value non-centrality parameter lambda
%  - df     : degrees of freedom
%  - cv     : critical value
%  - gam0   : right-hand probability
%
% Output:
%  - lam0ap : approximate value for non-centrality parameter lambda

% -----------------------------------------------------------------------------
% File......: lam0_app.m
% Author....: Marcel Martens
%             Mathematical Geodesy and Positioning
%             Delft University of Technology
% Date......: 20-dec-1998
% Modified..: 16-MAY-2000, P. Joosten MGP/TUD
% -----------------------------------------------------------------------------

a    = df + lambda;
b    = lambda/a;
nn   = a/(1.d0 + b);
q29n = 2.d0/(9.d0*nn);
xn   = ((cv/a)^(1.d0/3.d0) - (1.d0 - q29n))/sqrt(q29n);

if xn == 0.d0;
  qn = 0.5d0;
else
  XX = 0.5d0*xn*xn;
  qn = (1.d0-gamcdf(XX,0.5,1))/2.d0;
end;

if xn <= 0.d0; qn= 1.d0 - qn; end;

lam0ap = qn - gam0;

% --------------------------------
% --- End of function lam0_app ---
% --------------------------------
