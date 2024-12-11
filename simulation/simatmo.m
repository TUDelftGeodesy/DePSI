function aps=simatmo(size, beta, regime,scale,info);

% Generate an APS for SAR acquisition (1pix = 20x4 m) from fracsurfatmo.m
%
% aps           atmospheric phase screen [cyc]
%
% size          number of pixel
% beta          fractal dimension
% regime        cumulative percentage of each fractal 
%               dimension in the image
% info          'y' : graphic output
%               'n' : no graphic output
%
% try:
% aps=atmo(128, [5/3 8/3 2/3], [95 99 100],0.01,'y');

%if size~=128 fprintf('WARNING: the script works only wiiht size=128\n');return;end;

res = fracsurfatmo(size, beta, regime,'n',scale);

upbound=ceil(size/4.5);

restemp=res(:, 1:upbound);
upsample=5;

x  = [1:upsample:upsample*upbound];
y  = [1:1:size]';
xi = y';
yi = y;

% upsampling
aps =interp2(x,y,restemp,xi,yi,'linear');
aps=aps(:,1:size);


if strcmp(info,'y'), 
    figure;imagesc(aps);colormap(gray);colorbar;set(gca,'YDir','normal');title('Atmospheric Phase Screen [m]')
  end



