function c = cone(nrows,ncolumns)
% CONE --  Generate unit cone.
%   c = CONE;        cone [0,1] of dimensions (512,512);
%   c = CONE(L);     cone [0,1] of dimensions (L,L);
%   c = CONE([L P]); cone [0,1] of dimensions (L,P);
%   c = CONE(L,P);   cone [0,1] of dimensions (L,P);
% 
%   CONE is sharp, not rounded.
%
%   See also PYRAMID, RAMP, SIMINTERF 
% 

%// $Revision: 1.3 $  $Date: 2001/09/28 14:24:30 $
%// Bert Kampes, 11-Dec-2000

%%% Handle input.
if     (nargin==2) ;
elseif (nargin==1) if (prod(size(nrows))==1) ncolumns=nrows;
		   else ncolumns=nrows(2); nrows=nrows(1); end;
elseif (nargin==0) nrows=512; ncolumns=nrows;
else   error('wrong number of input');
end;

%%% Cone generation (can smarter?)
[x,y] = meshgrid(1:ncolumns,1:nrows);
x0    = round(size(x,2)/2);%          top of cone
y0    = round(size(x,1)/2);%          top of cone
x     = (x-x0).^2;
y     = (y-y0).^2;
c     = -sqrt(x+y);

%%% Normalize [-sqrt(x0):0] --> [0,1]
%minc  = min(c(:));
minc  = -sqrt((size(x,2)-x0).^2 + (size(x,1)-y0).^2);
c     = (c-minc) ./ -minc;

%%% EOF

