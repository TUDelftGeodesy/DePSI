function [lam0ac] = lam0_acc(lambda,df,gam0,cv)
%LAM0_ACC: computation of the accurate value of the non-centrality
%          parameter lambda
%
% This routine computes the accurate value of the non-centrality parameter
% lambda. It is based on a FORTRAN routine CHILNC.F by F.Kenselaar.
%
% Input:
%  - lambda : approximate value non-centrality parameter lambda
%  - df     : degrees of freedom
%  - cv     : critical value
%  - gam0   : right-hand probability
%
% Output:
%  - lam0ac : accurate value for non-centrality parameter lambda

% -----------------------------------------------------------------------------
% File......: lam0_acc.m
% Author....: Marcel Martens
%             Mathematical Geodesy and Positioning
%             Delft University of Technology
% Date......: 20-dec-1998
% Modified..: 16-MAY-2000, P. Joosten MGP/TUD
% -----------------------------------------------------------------------------
    
if lambda <0.d0;  lambda = 0.d0; end;
lam0ac = (1.d0 - ncx2cdf (cv,df,lambda)) - gam0;

% --------------------------------
% --- End of function lam0_acc ---
% --------------------------------
