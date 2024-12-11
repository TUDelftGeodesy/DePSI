%   VCEM3: Computes the (Co)Variance Factors Based on the Modified
%   Teunissen's Formula (1988) in Terms of A-Model. 
%   The Inpute of the Function Is as follows:
%       Function: vcem3(A,y,Qy1,sig0,Qy0)
%       A: Design (Observation) Matrix (m by n)
%       y: Vector of Observations (m)
%       Qy1: p-Number of m by m Cofactor Matrices in a Three Dimentional
%            Array as Qy1(:,:,k), k=1,2,...,p. (m by m by p)
%       sig0: A-Priori (Co)Variance Factors (Default Ones) (p)
%       Qy0: Constant Covariance Matrix (Default Zeros) (m by m)
%
% -------------------------------------------------------
% Function: vcem3
% Date    : 8 Jan. 2003
% Version : 1.1
% Author  : AliReza Amiri-Simkooei
%           Mathematical Geodesy and Positioning
%           Delft University of Technology
% --------------------------------------------------------

function sig2=vcem3(A,y,Qy1,sig0,Qy0)

%-----------------------------------------------
%--- Check Existance and Dimensions of Input ---
%-----------------------------------------------

if exist('A')~=1 error('Missing design matrix A.'); end
if exist('y')~=1 error('Missing vector of observations y.'); end
if exist('Qy1')~=1 error('Missing cofactor matrices Qy1.'); end
[m1,n]=size(A);
[m2,m3,p1]=size(Qy1);
if m2~=m3 error('Matrix Qy1 is not square.'); end
if m1~=m2 error('Incorrect size of matrices Qy1 and A.'); end
m=m1; clear m1 m2 m3
if exist('sig0')==1
    [p2,h1]=size(sig0);
    if p1~=p2 | h1~=1 error('Incorrect size of vector sig0.');
    else p=p1; clear p1; end
else p=p1; clear p1; sig0=ones(p,1); end
if exist('Qy0')==1
    [m1,m2]=size(Qy0);
    if m1~=m | m2~=m error('Incorrect size of matrix Qy0.');
    else clear m1 & m2; end
else Qy0=zeros(m,m); end

%-----------------------------------
%--- Check the validity of input ---
%-----------------------------------

Qy=Qy0;
for k=1:p
    Qy=Qy+sig0(k)*Qy1(:,:,k);
end
%---------------------------------
%--- Compute the normal matrix ---
%---------------------------------

Sigtest=inf;
Sighat=sig0;
sig1=sig0;
count1=0;
%for iiii=1:6
while abs(Sighat-Sigtest)>1e-6
    count1=count1+1;
    if count1>200 error('Too many number of iterations (larger than 200).');  end
Sigtest=Sighat;
Qy=Qy0;
for k=1:p
    Qy=Qy+sig0(k)*Qy1(:,:,k);
end
Qyinv=inv(Qy);
Pao=eye(m)-A*inv(A'*Qyinv*A)*A'*Qyinv;
for k=1:p
    l(k,1)=0.5*y'*Qyinv*Pao*Qy1(:,:,k)*Qyinv*Pao*y;
    for j=1:p
        N(k,j)=0.5*trace(Qyinv*Pao*Qy1(:,:,k)*Qyinv*Pao*Qy1(:,:,j));
    end
end
Ninv=inv(N);
Sighat=Ninv*l;
sig0=Sighat;

end

sig2=Sighat;