function ps_data = ps_quadtree_freek(ps_data);

% Input:    - ps_data
%
% Output:   - ps_data
%
% ----------------------------------------------------------------------
% File............: ps_quadtree.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Authors.........: Astrid Humme
%                   Freek van Leijen
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


global project_id Npar_max az_spacing r_spacing fig

Nps_temp = size(ps_data,1);


% -----------------------------------------------------------------------
% Quadtree decomposition
% -----------------------------------------------------------------------

x = ps_data(:,2)*az_spacing;
y = ps_data(:,3)*r_spacing;
z = ps_data(:,7);

%Nps_block = max(round(Nps_temp/1000),10);
%[ind,bx,by,Nb,lx,ly] = quadtree_freek(x,y,[],Nps_block);
threshold = (max(z)-min(z))/10;
[ind,bx,by,Nb,lx,ly] = quadtree_freek(x,y,z,[],threshold);
lx_orig = lx;
ly_orig = ly;

fig = fig+1;
figure(fig);plot(lx_orig(:,[1 2 2 1 1])',ly_orig(:,[1 1 2 2 1])','k');hold on; scatter(ps_data(:,2)*az_spacing,ps_data(:,3)*r_spacing,5,ps_data(:,7),'filled');hold on

% -----------------------------------------------------------------------
% Data reduction
% -----------------------------------------------------------------------

index = 1; 
%while Nps_temp>5000& ~isempty(index)
while Nps_temp>1000 & ~isempty(index)
  block_area = (lx(:,2)-lx(:,1)).*(ly(:,2)-ly(:,1));
  ind2 = find(block_area==min(block_area));
  for v = 1:length(ind2)
    ind3 = find(ind == ind2(v));
    if ~isempty(ind3)
      median_data = median(ps_data(ind3,:),1);
      median_data(1:3) = round(median_data(1:3));
      ps_data(ind3,:) = [];
      ps_data(end+1,:) = median_data;
      ind(ind3) = [];
      ind(end+1) = ind2(v);
    end
  end
  lx(ind2,:) = NaN;
  ly(ind2,:) = NaN;
  index = find(~isnan(lx(:,1)));
  Nps_temp = size(ps_data,1);
end

fig = fig+1;
figure(fig);plot(lx_orig(:,[1 2 2 1 1])',ly_orig(:,[1 1 2 2 1])','k');hold on; scatter(ps_data(:,2)*az_spacing,ps_data(:,3)*r_spacing,5,ps_data(:,7),'filled');hold on
