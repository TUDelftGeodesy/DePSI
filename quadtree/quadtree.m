function  [ind,bx,by,Nb,lx,ly] = quadtree(x,y,s,n0)
% QUADTREE  Recursive division of a 2-dimensional set.
%	[IND,BX,BY,NB,LX,LY] = QUADTREE(X,Y,S,N0)
%	Performs recursive tree-like division of a set
%	of points with coordinates  X,Y.
%	S is binary mask showing which points of a set
%	are to be counted. N0 is maximum permissible
%	number of "counted" points in the elementary
%	block.
%	Returns vector IND of the same size as X, Y
%	showing which region each point of a set belongs
%	to; binary matrices BX, BY where each row shows
%	"binary address" of each region.
%	Also returns "Adjacency matrix" NB which is 1 if 
%	i and j regions are neighbours and 0 otherwise;
%	and matrices of limits for each region, LX and LY
%	so that  PLOT(LX(:,[1 2 2 1 1])',LY(:,[1 1 2 2 1])')
%	plots the boundaries of all regions.

%  Copyright (c) 1995  by Kirill K. Pankratov
%	kirill@plume.mit.edu
%	01/30/95

 % Default for maximal block size 
n0_dflt = 10;

 % Handle input ..............................
if nargin < 2
  error(' Not enough input arguments.')
end
if nargin<4, n0 = n0_dflt; end
if nargin<3, s = []; end
if isempty(s), s = ones(size(x)); end                   % s ==[] wordt isempty(s)

 % If no need to divide ......................
if length(x(find(s))) <= n0
  bx = []; by = []; Nb = [];
  ind = ones(size(x));
  return
end

 % Calculate limits ................
lim = [min(x) max(x) min(y) max(y)];

 % Recursively construct quadtree ............
[ind,bx,by] = qtree0(x,y,s,lim,n0);
bx = bx(:,1:size(bx,2)-1);
by = by(:,1:size(by,2)-1);

 % Compose "adjacency matrix" ................
szb = size(bx);
ob = ones(1,szb(1));
pow = 2.^(0:szb(2)-1);
pow = flipud(pow');

 % "Lower" and "upper" bounds for trees
bxmax = ceil(bx)*pow;
bxmin = floor(bx)*pow;
bymax = ceil(by)*pow;
bymin = floor(by)*pow;

 % "Physical" limits of each regions
lx = [bxmin bxmax+1];
ly = [bymin bymax+1];
lx = lim(1)+lx*(lim(2)-lim(1))/pow(1)/2;
ly = lim(3)+ly*(lim(4)-lim(3))/pow(1)/2;

B0 = bxmin(:,ob);
B1 = bxmax(:,ob);
Nb = (B0'-B1<=1) & (B1'-B0>=-1);

B0 = bymin(:,ob);
B1 = bymax(:,ob);
Nb = Nb & (B0'-B1<=1) & (B1'-B0>=-1);

