function [Chi2,afixed] = ps_chistart(D,L,a,ncands,qfactor)
%CHISTART: Computes the initial size of the search ellipsoid
%
% This routine computes or approximates the initial size of the search
% ellipsoid. If the requested number of candidates is not more than the
% dimension + 1, this is done by computing the squared distances of partially
% conditionally rounded float vectors to the float vector in the metric of the
% covariance matrix. Otherwise an approximation is used.
%
% Input arguments
%    L,D   : LtDL-decomposition of the variance-covariance matrix of
%            the float ambiguities (preferably decorrelated)
%    a     : float ambiguites (preferably decorrelated)
%    ncands: Requested number of candidates (default = 2)
%    qfactor: Multiplication factor for the volume of the resulting
%            search ellipsoid (default = 1.5)
%
% Output arguments:
%    Chi2  : Size of the search ellipsoid

% ----------------------------------------------------------------------
% File.....: chistart.m
% Date.....: 19-MAY-1999
% Modified.: 05-MAR-2001, by P. Joosten
% Author...: Peter Joosten
%            Mathematical Geodesy and Positioning
%            Delft University of Technology
% ----------------------------------------------------------------------

% 12-03-09, changed by FvL to output the best afixed
% 05-02-12, FvL, bug fix for k=i=n

% ------------------
% --- Initialize ---
% ------------------

if nargin < 4; ncands = 2  ; end;
if nargin < 5; qfactor = 1.5; end;

n = max(size(a));

% ----------------------------------------------------------------------
% --- Computation depends on the number of candidates to be computed ---
% ----------------------------------------------------------------------

if ncands <= n+1;

  % --------------------------------------------------------
  % --- Computation based on the bootstrapping estimator ---
  % --------------------------------------------------------

  %Chi = [];
  %%% SPEED UP BK
  Chi = zeros(n+1,1);
  invLtDL = inv(L.'*diag(D)*L);
  
  afixed_array = NaN(n,n+1);
  
  for k = n:-1:0;

    afloat = a;
    afixed = a;
  
    for i = n:-1:1;
      %      dw = 0;
      %      for j = n:-1:i;
      %	dw = dw + L(j,i) * (afloat(j) - afixed(j));
      %      end;
      % speed up use mat mult
      dw = L(i:n,i).' * (afloat(i:n) - afixed(i:n));% standing a
      afloat(i) = afloat(i) - dw;
      if (i ~= k);
	afixed(i) = round (afloat(i));
      else;
	%if isequal (afloat(i),afixed(i));
	if (afloat(i) == afixed(i));
	  %afixed(i) = afixed(i) + 1;
	  afixed(i) = round (afloat(i) + sign(afloat(i)));%FvL, 05-02-12
	else;
	  afixed(i) = round (afloat(i) + sign (afloat(i) - afixed(i)));
	end;
      end;
    end;
  
    %%% speed up BK
    q = a-afixed;
    Chi(k+1) = q.' * invLtDL * q;
    %Chi = [Chi (a-afixed)' * inv(L'*diag(D)*L) * (a-afixed)];
  
    afixed_array(:,k+1) = afixed; %FvL, from original
    
  end;

  % ---------------------------------------------------------------
  % --- Sort the results, and return the appropriate number     ---
  % --- Add an "eps", to make sure there is no boundary problem ---
  % ---------------------------------------------------------------

  %keyboard speed up sort BK
if (ncands == 1)
  [Chi,min_index] = min(Chi); %FvL
  Chi2 = Chi + 1d-6; %FvL
  afixed = afixed_array(:,min_index); %FvL
else
  [Chi,sort_index]  = sort(Chi); %FvL
  Chi2 = Chi(1:ncands) + 1d-6; %FvL
  afixed = afixed_array(:,sort_index(1:ncands)); %FvL
end
%warning('no qfactor?');   
else

  % -----------------------------------------------------
  % An approximation for the squared norm is computed ---
  % -----------------------------------------------------
  
  Linv = inv(L);
  Dinv = 1./D;
   
  Vn   = (2/n) * (pi ^ (n/2) / gamma(n/2));
  Chi2 = qfactor * (ncands / sqrt((prod(1 ./ Dinv)) * Vn)) ^ (2/n);
  afixed = [];
end;

% ----------------------------------------------------------------------
% End of routine: chistart
% ----------------------------------------------------------------------

