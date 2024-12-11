function [topo]= simtopo(N);

% [topo] = simtopo(N)
%
% INPUT
%
% N       size (NxN) of the image
%
% OUTPUT
% 
% topo    topography


% generate a 3D fractal surface of size NxN and fractal dim 2.2 
topo = fractopobis(1, N);








