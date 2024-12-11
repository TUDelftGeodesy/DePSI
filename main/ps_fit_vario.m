function [a,b,c1,c2,emp_vario] = ps_fit_vario(x,y,z,model,a0,b0,c10,c20)

% Estimation of variogram model
%
% Input:    - x                 x-coordinates
%           - y                 y-coordinates
%           - z                 z-values
%           - model             variogram model:
%                               2 = exponential
%                               3 = gaussian
%                               4 = spherical
%           - a0                initial value range
%           - b0                initial value hole effect
%           - c10               initial value sill
%           - c20               initial value nugget
%
% Output:   - a                 estimated range
%           - b                 estimated hole effect
%           - c1                estimated sill
%           - c2                estimated nugget
%           - emp_vario         empirical variogram
%
% ----------------------------------------------------------------------
% File............: ps_fit_vario.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Authors.........: Freek van Leijen
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
% v1.7.2.15, Freek van Leijen
% - change in distance computation and maximum distance to enable
% larger number of PS1
% v1.7.8.0, Freek van Leijen
% - use of hypot function to increase speed
%


% -----------------------------------------------------
% Initialize
% -----------------------------------------------------

global fig detail_plots visible_plots


% -----------------------------------------------------
% Calculate raw variogram
% -----------------------------------------------------

Nx = size(x,1);
h1 = NaN(Nx^2,1,'single');
h2 = NaN(Nx^2,1,'single');
count = 0;
for v = 1:Nx
  temp = hypot(x-x(v),y-y(v));
  temp_index = find(temp>0&temp<10000); %new, max 10 km
  Ntemp = length(temp_index);
  h1(count+1:count+Ntemp) = temp(temp_index);
  h2(count+1:count+Ntemp) = 0.5*((z(temp_index)-z(v)).^2);
  count = count+Ntemp;
end
h1 = h1(1:count);
h2 = h2(1:count);

if strcmp(detail_plots,'y')
  fig = fig+1;
  figure(fig);hold on
  if strcmp(visible_plots,'n')
    set(gcf,'visible','off');
  end
  plot(h1,h2,'*');
  xlabel('Distance [m]')
  ylabel('Variogram [rad^2]')
end


% -----------------------------------------------------
% Calculate experimental variogram
% -----------------------------------------------------

max_dist = max(h1(:));
lags = [0:max_dist/50:max_dist];
Nlags = length(lags)-1;
emp_vario = NaN(Nlags,3);
for v = 1:Nlags
    index = find((h1>=lags(v))&(h1<lags(v+1)));
    emp_vario(v,1) = length(index);
    if emp_vario(v,1)~=0
      emp_vario(v,2) = mean(h2(index));
    end
    emp_vario(v,3) = (lags(v)+lags(v+1))/2;
end
index = find(isnan(emp_vario(:,2)));
emp_vario(index,:) = [];

if strcmp(detail_plots,'y')
  figure(fig);
  plot(emp_vario(:,3),emp_vario(:,2),'r*','markersize',10)
end

% added by Prabu
if size(emp_vario,1) < round(Nlags/2)
    a = NaN;
    b = NaN;
    c1 = NaN;
    c2 = NaN;
    emp_vario = NaN;
    return;
end
% added by Prabu

emp_vario = emp_vario(1:round(Nlags/2),:);
Nlags = size(emp_vario,1);


% -----------------------------------------------------
% Estimate experimental variogram
% -----------------------------------------------------

lam = [a0 c10 c20];

switch model
  case 1
    error(['You can not estimate the nugget independently, choose ' ...
           'another model (nugget included)']);
  case 2 %exponential
    [lam,resnorm,residual,exitflag] = lsqcurvefit(@ps_variogram_exp,lam,emp_vario(:,3),emp_vario(:,2),[0 0 0],[],optimset('Display','off'));
  case 3 %gaussian
    [lam,resnorm,residual,exitflag] = lsqcurvefit(@ps_variogram_gaus,lam,emp_vario(:,3),emp_vario(:,2),[0 0 0],[],optimset('Display','off'));
  case 4 %spherical
    [lam,resnorm,residual,exitflag] = lsqcurvefit(@ps_variogram_spher,lam,emp_vario(:,3),emp_vario(:,2),[0 0 0],[],optimset('Display','off'));
  otherwise
    error('You specified a wrong variogram model');
end
    
%if exitflag==0
%  warning('Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.');
%elseif exitflag==-1
%  warning('Algorithm was terminated by the output function.');
%elseif exitflag==-2
%  warning('Problem is infeasible: the bounds lb and ub are inconsistent.');
%elseif exitflag==-3
%  warning('Optimization could not make further progress.');
%end

a = lam(1);
b = 0;
c1 = lam(2);
c2 = lam(3);


if strcmp(detail_plots,'y')
  [y_final,dyda] = ps_variogram(model,emp_vario(:,3),a,b,c1,c2);
  figure(fig);
  plot(emp_vario(:,3),y_final,'r','linewidth',2);
end



% subfunctions

function y = ps_variogram_exp(x,xdata);

y = x(2)*(1-exp(-xdata/x(1)))+x(3);


function y = ps_variogram_gaus(x,xdata);

y = x(2)*(1-exp(-(xdata/x(1)).^2))+x(3);


function y = ps_variogram_spher(x,xdata);

y = x(2)*(1.5*min(xdata/x(1),1)-0.5*min(xdata/x(1),1).^3)+x(3);
