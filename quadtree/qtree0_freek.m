function  [ind,bx,by] = qtree0_freek(x,y,z,s,lim,n0)
%	QTREE0  Primitive for QUADTREE.
%	[IND,BX,BY] = QTREE0(X,Y,S,LIM,N0)	

%  Copyright (c) 1995  by Kirill K. Pankratov
%	kirill@plume.mit.edu
%	01/26/95


 % Divide and make further recursion if necessary
N = length(find(s));

 % If not to divide further .....................
ind = ones(size(x));
zmin = min(z);
zmax = max(z);
dz = zmax-zmin;

if isempty(dz)||dz<n0
  %if N <= n0
  bx = .5; by = .5;
  return
end

 % Further division is needed ...................

 % Calculate middle of the current range
x_mid = (lim(2)+lim(1))/2;
y_mid = (lim(4)+lim(3))/2;

 % Branch for current level
bx0 = [0 0 1 1]';
by0 = [0 1 0 1]';
 
n1 = find( x<=x_mid & y<=y_mid );   % Lower-left
lm = [lim(1)  x_mid  lim(3)  y_mid];
[ind1,bx1,by1] = qtree0_freek(x(n1),y(n1),z(n1),s(n1),lm,n0);

n2 = find( x<x_mid & y>y_mid );     % Upper-left
lm = [lim(1)  x_mid  y_mid  lim(4)];
[ind2,bx2,by2] = qtree0_freek(x(n2),y(n2),z(n2),s(n2),lm,n0);

n3 = find( x>x_mid & y<=y_mid );    % Lower-right
lm = [x_mid  lim(2)  lim(3)  y_mid];
[ind3,bx3,by3] = qtree0_freek(x(n3),y(n3),z(n3),s(n3),lm,n0);

n4 = find( x>x_mid & y>y_mid );     % Upper-right
lm = [x_mid  lim(2)  y_mid  lim(4)];
[ind4,bx4,by4] = qtree0_freek(x(n4),y(n4),z(n4),s(n4),lm,n0);

 % Make composite binary "tree matrices"
szB = [size(bx1); size(bx2); size(bx3); size(bx4)];
szv = cumsum(szB(:,1));
szh = max(szB(:,2))+1;

 % Initialize output matrices
bx = .5*ones(szv(4),szh);
by = bx;

 % Fill binary matrices ..................
nnv = (1:szv(1))';        nnh = 2:szB(1,2)+1;
bx(nnv,nnh) = bx1; by(nnv,nnh) = by1;
bx(nnv,1) = bx0(1)*ones(size(nnv));
by(nnv,1) = by0(1)*ones(size(nnv));

nnv = (szv(1)+1:szv(2))'; nnh = 2:szB(2,2)+1;
bx(nnv,nnh) = bx2; by(nnv,nnh) = by2;
bx(nnv,1) = bx0(2)*ones(size(nnv));
by(nnv,1) = by0(2)*ones(size(nnv));

nnv = (szv(2)+1:szv(3))'; nnh = 2:szB(3,2)+1;
bx(nnv,nnh) = bx3; by(nnv,nnh) = by3;
bx(nnv,1) = bx0(3)*ones(size(nnv));
by(nnv,1) = by0(3)*ones(size(nnv));

nnv = (szv(3)+1:szv(4))'; nnh = 2:szB(4,2)+1;
bx(nnv,nnh) = bx4; by(nnv,nnh) = by4;
bx(nnv,1) = bx0(4)*ones(size(nnv));
by(nnv,1) = by0(4)*ones(size(nnv));
  
 % Calculate indices .....................
ind2 = ind2+szv(1);
ind3 = ind3+szv(2);
ind4 = ind4+szv(3);
ind(n1) = ind1;
ind(n2) = ind2;
ind(n3) = ind3;
ind(n4) = ind4;
