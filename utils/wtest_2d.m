function [w]= wtest_2d(Qecheck,echeck,Nifgs)

% 2D w-test (spatial-temporal)
%
% Input:    - Qecheck             diagonal of variance matrix of the residuals
%           - echeck              residuals
%           - Nifgs               number of interferograms
%
% Output:   - w                   w-tests
%
% ----------------------------------------------------------------------
% File............: wtest_2d.m
% Version & Date..: 1.6, 31-JAN-2006
% Author..........: Freek van Leijen
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
% change log
% v1.7.2.5, only diagonal of Qecheck needed, no Qy involved (Shizhuo Liu)
%


% remove lose edges, which cannot be tested
index = find(abs(Qecheck) < 1e-10);
if ~isempty(index)
  echeck(index,:) = [];
  Qecheck(index) = [];
  %Qecheck(index:) = [];
  %Qecheck(:,index) = [];
  %invQy(index,:) = [];
  %invQy(:,index) = [];
end

% testing
% H4     = invQy*echeck;
%H4 = echeck ; % since invQy is identity matrix
%H5     = repmat(diag(invQy*Qecheck*invQy),1,Nifgs); %somehow, this gave a numerical error sometimes when invQy=I
%H5     = repmat(diag(Qecheck),1,Nifgs);
%H5     = repmat(Qecheck,1,Nifgs);
%H5(H5<0)=0;

Qecheck = repmat(Qecheck,1,Nifgs);
Qecheck(Qecheck<0)=0;

w = abs(echeck./sqrt(Qecheck)); % w-test
if ~isempty(index)
  for v = 1:length(index)
    w = [w(1:index(v)-1,:);NaN(1,Nifgs);w(index(v):end,:)];
  end
end

