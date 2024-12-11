function [lam0] = lambda0(gam0, df, cv)
%LAMBDA0: computation of the non-centrality parameter lambda
%
% This routine computes the non-centrality parameter lambda out of the
% non-central chi-square distribution by bisection iteration of the
% probability function. It is based on CHILNC.F by F.Kenselaar.
%
% The routine depends on the availability of:
% - lam0_acc and lam0_app
% - Matlab's statistics toolbox
%
% Input:
%  - gam0   : right-hand probability
%  - df     : degrees of freedom
%  - cv     : critical value
%
% Output:
%  - lam0   : non-centrality parameter lambda

% -----------------------------------------------------------------------------
% File......: lambda0.m
% Author....: Marcel Martens
%             Mathematical Geodesy and Positioning
%             Delft University of Technology
% Date......: 20-dec-1998
% Modified..: 16-MAY-2000, P. Joosten MGP/TUD
% -----------------------------------------------------------------------------

% -----------------------
% --- Initializations ---
% -----------------------

lam0  = 0.d0;
stop1 = 0;
stop2 = 0;
mxit  = 50;
fac   = 1.6D0;

% -----------------------
% --- Checks on input ---
% -----------------------
if df   <  1                  ; error('LAMBDA0: degree of freedom smaller than 1.'); end;
if cv   <  0.e-9              ; error('LAMBDA0: critical value less than zero.'); end;
if gam0 <= 0.d0 | gam0 > 1.d0 ; error('LAMBDA0: gamma0 not between 0 and 1'); end;

cv1 = chi2inv ((1.d0-gam0),df);

if cv <  cv1 ; error('LAMBDA0: this gamma0 and Critical Value lead to negative lambda'); end;
if cv == cv1 ; lam0 = 0.d0; stop1 = 1 ; end;


if stop1 == 0;

  % ---------------------------------------------
  % --- compute startvalue for lambda (lstrt) ---
  % ---------------------------------------------

  x1 = 0;
  x2 = 5;
  ii = 0;
  f1 = lam0_app(x1, df, gam0, cv);
  f2 = lam0_app(x2, df, gam0, cv);
  step = x2 - x1;

  while ii <=  mxit & stop2==0;

    ii=ii+1;

    if f1.*f2 < 0.d0; stop2=1; end;

    if abs(f1) < abs(f2);
      x1 = x1 - step;
      f1 = lam0_app(x1,df,gam0,cv);
    else
      x2 = x2 + step;
      f2 = lam0_app(x2,df,gam0,cv);
    end;
    step = step.*fac;
  end;

  lstrt = (x1+x2)/2.d0;

  % -------------------------------------------------------------------
  % --- compute approximate value for lambda using an approximation ---
  % --- formula for non-central chi-square distribution             ---
  % -------------------------------------------------------------------
  
  if (cv-cv1) < 0.1d0;
    lambda = 0.d0;
  else
    options = optimset('TolX',1.e-4);
    lambda = fzero(@(lambda) lam0_app(lambda,df,gam0,cv),lstrt,options);
  end

  % -------------------------------------------
  % -- compute an accurate value for lambda ---
  % -------------------------------------------
  options = optimset('TolX',1.e-6);
  lam0 = fzero(@(lambda) lam0_acc(lambda,df,gam0,cv),lambda,options);  %fzero('lam0_acc',lambda,1.e-6,[],df,gam0,cv);

end;

% -------------------------------
% --- End of function lambda0 ---
% -------------------------------
