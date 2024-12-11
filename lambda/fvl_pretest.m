function [lam0,k1,kb,alphab] = fvl_pretest(b,alpha0,gamma0)

%PRETEST: perform preparations for statistical hypothesis testing
%
% This routine performs some preparations to perform hypothesis testing
% using 'testa' or 'testb'
%
% Following parameters need to exist:
%   b     : redundancy (number of conditions)
%   alpha0: level of significance, optional, default = 0.001
%   gamma0: power, optional, default = 0.80
%
% After this script, following parameters of interest are available:
%   lam0  : non-centrality parameter
%   k1    : critical value for 1-dimensional test
%   kb    : critical value for b-dimensional test
%   alphab: level of significance for b-dimensional test

% ----------------------------------------------------------------------
% File.....: pretest.m
% Date.....: 18-MAY-2000
% Version..: 1.0
% Author...: Christian Tiberius / Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% --------------------------------------------------------------------
% --- check if input is available and correct, set defaults if not ---
% --------------------------------------------------------------------

if exist('b')      ~= 1; disp ('redundancy parameter b missing'); end;
if exist('alpha0') ~= 1; alpha0=0.001; end; % set default for level of sign.
if exist('gamma0') ~= 1; gamma0=0.800; end; % set default for power

if (size(alpha0,1) ~= 1) | (size(alpha0,2) ~= 1); disp('alpha0 not a scalar'); end
if (size(gamma0,1) ~= 1) | (size(gamma0,2) ~= 1); disp('gamma0 not a scalar'); end
 
if (alpha0 <= 0) | (alpha0 >= 1); disp('alpha0 not in <0,1>'); end;
if (gamma0 <= 0) | (gamma0 >= 1); disp('gamma0 not in <0,1>'); end;

% ------------------------------------------------------
% --- compute non-centrality parameter lambda nought ---
% ------------------------------------------------------

k1   = chi2inv(1.0-alpha0,1);          % critical value for b=1 dof (chi-squared)
lam0 = lambda0(gamma0,1,k1);           % non-centrality parameter

% -----------------------------------------
% --- compute critical values k1 and kb ---
% -----------------------------------------

if (b>0);
  
  kb     = ncx2inv(1.0-gamma0,b,lam0); % critical value for b dof
  alphab = (1.0-chi2cdf(kb,b));        % lev. of sign. for b-dimensional test

else

  disp('message: zero redundancy');
  kb     = NaN;
  alphab = NaN;

end;

k1 = sqrt(k1);                           % critical value for b=1 dof (normal)

% -------------------------------
% --- End of script pretest.m ---
% -------------------------------
