function [calfactors] = empcal_nonlin_Q(PSCmat_in,varda_in,mean_amp_in)

% Task: Matlab implementation of appendix B of Bianca Cassee's thesis
% By: GKE
% Date: 03/01/2005
%
% Input:
%
% Output:
%         
% ----------------------------------------------------------------------
% File............: empcal_nonlin_Q.m
% Version & Date..: 1.6, 31-JAN-2006
% Author..........: Gini Ketelaar
%                   Delft Institute of Earth Observation and Space Systems
%                   Delft University of Technology
% ----------------------------------------------------------------------
%
% This software is developed by Delft University of Technology and is
% intended for scientific use only. Applications for commercial use are
% prohibited.
%  
% Copyright (c) 2004-2005 Delft University of Technology, The Netherlands
%


% assumption that Qy is blkdiag

% Input from: Ampitude observation matrix
% Input data: A(p,k), K SAR images, P PSC's.

calfactors_array = NaN(10,size(PSCmat_in,2));

for v = 1:20
  
  rand('state',sum(100*clock))

  %cal_index = ceil(rand(100,1)*size(PSCmat_in,1));
  cal_index = ceil(rand(50,1)*size(PSCmat_in,1));
  PSCmat = PSCmat_in(cal_index,:);
  varda = varda_in(cal_index);
  mean_amp = mean_amp_in(cal_index);
  
%%%%%%%%%%%%%%%%%%%
% Data simulation %
%%%%%%%%%%%%%%%%%%%
varda_orig=varda;

it=1;
maxit=9;
Pall=size(PSCmat,1);
selected_P=(1:Pall)';
Tqka=9999;
sig2hat=0;
[P K]=size(PSCmat);
mintt=0.95; maxtt=ncfinv(1-0.001,K-1,1e5,0); % minimum and maximum
                                             % values for
                                             % testratio's
                                             % changed from 1e8 to
                                             % 1e5 to avoid warning [FvL]
maxtt;
maxomtt=1.05;

while and(  (it<=maxit)  ,or ( ((and(sig2hat<maxomtt,sig2hat>mintt))==0)     , (max(Tqka)>maxtt)    )   )
%%%disp(['Maxtqka ',num2str(max(Tqka))]);
    
close all;
[P K]=size(PSCmat);
A=zeros(K*P,K+P);
dy=zeros(K*P,1);
Qy=repmat( reshape(repmat(varda',K,1),1,K*P), size(A,2), 1);

x=zeros(K+P,1);
x(1:K,1)=1; x( (K+1):end, 1 )=mean_amp; %PSCmat(:,1); % initial unknown values
for p=1:P
    A( ((p-1)*K+1) : (p*K), 1:K ) = mean_amp(p)*eye(K); %PSCmat(:,1); % Aref    
    A( ((p-1)*K+1) : (p*K), K+p ) = ones(K,1); % take one image for initial values
    dy( ((p-1)*K+1) : (p*K), 1) = PSCmat(p,:)' - mean_amp(p); %PSCmat(p,1); % y-y0
end

% This matrix has a rankdefect of 1, take the scaling factor of the first
% image fixed
A = A(:,2:end);
Qy=Qy(2:end,:);
x = x(2:end,:);
%
% Initial values for unknowns
A=sparse(A);
Qyinv=1./Qy; % QY blockdiag
Qyinv=sparse(Qyinv);
clear Qy;
Qxinv=(A'.*Qyinv*A);
Qx = inv(Qxinv); %FvL
%Qx=inv(A'.*Qyinv*A);

xhat=Qx*A'.*Qyinv*dy;
criterion = xhat'*Qxinv*xhat; %FvL

while criterion > exp(-6)
    x=x+xhat;
    x_updated=[1;x]; % including basis (first scalefactor=1)
    A=zeros(K*P,P+K);
     for p=1:P
      A( ((p-1)*K+1) : (p*K), 1:K ) = x_updated(K+p)*eye(K); %PSCmat(:,1); % Aref    
      A( ((p-1)*K+1) : (p*K), K+p ) = x_updated(1:K); % take one image for initial values
      dy( ((p-1)*K+1) : (p*K), 1) = PSCmat(p,:)' - x_updated(1:K)*x_updated(K+p); % y-y0
     end
    A = A(:,2:end);
    Qxinv=(A'.*Qyinv*A);
    Qx = inv(Qxinv); % FvL
%    Qx=inv(A'.*Qyinv*A);
    xhat=Qx*A'.*Qyinv*dy;
    %adjusta; % lam0=17.075;testa;
    criterion = xhat'*Qxinv*xhat; % FvL
end

x=[1;x+xhat];

% ehat
ehat=zeros(K*P,1);
for p=1:P
    ehat( ((p-1)*K+1) : (p*K), 1) =PSCmat(p,:)' - x(1:K)*x(K+p); % y-y0
    yhat( ((p-1)*K+1) : (p*K), 1) =x(1:K)*x(K+p);
end

%Qehat
temp=zeros(P*K);
temp(sub2ind([P*K,P*K],1:P*K,1:P*K))=1./Qyinv(1,:);
Qehat=temp-A*Qx*A';

clear temp;

% w-tests
diagQehat=diag(Qehat);
w=ehat./sqrt(diagQehat); % assumption that Qy is blockdiag
%OMT
sig2hat=sum( ((ehat.^2).*(Qyinv(1,:))') ./(size(A,1)-size(A,2)) ); 
OMT=sum ((ehat.^2).*(Qyinv(1,:))');

% point-tests
temp=0;
clear Tq;
for p=1:P
 Qehat_sel= Qehat ( ((p-1)*K+2) : (p*K),  ((p-1)*K+2) : (p*K) );
 ehat_sel=   ehat ( ((p-1)*K+2) : (p*K), 1);
 Tq(p)=ehat_sel'*inv(Qehat_sel)*ehat_sel;
 clear temp Qehat_sel ehat_sel;
end

% testing procedure, pointtest, OMT
% determine the critical values
% for w^2
% for Tq(q=K-1)
% for sigp(q=round((K*P-(K+P-1))/P)))
%
% wtest, alfa=0.001
alfa=0.001;
pw=0.5;
ka_w2=ncx2inv(1-0.001,1,0);
la0=ka_w2;
ka_Tq=ncx2inv(pw,K-1,la0);
alfa_Tq=1-ncx2cdf(ka_Tq,K-1,0);
ka_OMT=ncx2inv(pw,K*P-(K+P-1),la0);

% Calculate testratio's
OMTka=OMT/ka_OMT;
Tqka=Tq/ka_Tq;

% plot w teststatistics
xp=-5:0.2:5;
[n,xp]=hist(w,xp);
n=n/(sum(n)*(xp(2)-xp(1)));
%figure(11);bar(xp,n,1);hold on;
%y=normpdf(-5:0.1:5,0,1); plot(-5:0.1:5,y,'r');
%xlabel('w-teststatistics');
% ADAPTION FROM TESTING
% removal of point directly influences stochastic model, therefore, point
% removal is followed by vce update
if rem(it,2)~=0
    varda=sig2hat*varda;
else
 if OMTka > max(Tqka)
    varda=sig2hat*varda;
 else
    [maxTq,ind]=max(Tqka); 
    mind=numel(Tq);
    if maxTq>maxtt
        if ind==1
            PSCmat=PSCmat(2:end,:);
            varda=varda(2:end,1);
            x=[x(1:K,1);x((K+2):end,1)];
            mean_amp=mean_amp(2:end,1);
            selected_P=selected_P(2:end);
        end
        if and(ind>1,ind<mind)
            PSCmat=[PSCmat(1:ind-1,:);PSCmat(ind+1:end,:)];
            varda=[varda(1:ind-1,1);varda(ind+1:end,1)];
            x=[x(1:(K+ind-1),1);x((K+ind+1):end,1)];
            mean_amp=[mean_amp(1:ind-1,1);mean_amp(ind+1:end,1)];
            selected_P=[selected_P(1:ind-1,1);selected_P(ind+1:end,1)];
        end
        if ind==mind
            PSCmat=PSCmat(1:(end-1),:);
            varda=varda(1:(end-1),1);
            x=x(1:(end-1),1);
            mean_amp=mean_amp(1:(end-1),1);
            selected_P=selected_P(1:(end-1),1);
        end
        varda=sig2hat*varda;
    else
        varda=sig2hat*varda;  % adapt stochastic model
    end
 end
end

calfactors=x(1:K,1)';
clear A Qyinv Qxinv y Pal Qyinv ehat Qehat vce_n vce_l vce_varhat

% check
[P K]=size(PSCmat);
clear p_da_cal p_da_nocal
for p=1:P
    temp=PSCmat(p,:)'./x(1:K);
    %p_da(p)=std(temp)/mean(temp); of course this goes wrong, especially
    %with the dopplers!
    p_da_cal(p)=sqrt(varda(p))/x(K+p);
    p_da_nocal(p)=std(PSCmat(p,:)')/mean(PSCmat(p,:)');
    clear temp;
end
%%%figure(99);hold on;
%%%plot(p_da_nocal,'y-*');
%%%plot(p_da_cal,'b-*');
%%%hold off;
%endcheck 

it=it+1;
end

% check if when calibrated, under threshold
[P K]=size(PSCmat);
clear p_da_cal p_da_nocal
for p=1:P
    temp=PSCmat(p,:)'./x(1:K);
    %p_da(p)=std(temp)/mean(temp); of course this goes wrong, especially
    %with the dopplers!
    p_da_cal(p)=sqrt(varda(p))/x(K+p);
    p_da_nocal(p)=std(PSCmat(p,:)')/mean(PSCmat(p,:)');
    clear temp;
end
 sel_index=find(p_da_cal<0.25);
 selected_P=selected_P(sel_index);
 p_da_cal=p_da_cal(sel_index);
 %%%disp(['Selected PS ',num2str(length(selected_P)),' of ',num2str(Pall)]);
    
clear Qx x;

%p_da
%p_da_nocal
%%%figure(100);hold on;
%%%plot(p_da_nocal,'r-o');
%%%plot(p_da_cal,'b-s');
%%%plot([0,40],[0.25,0.25],'k--','Linewidth',2);
%%%legend('No calibration','After calibration validation');
%%%xlabel('Point number');ylabel('Amplitude dispersion');
%%%hold off;

calfactors_array(v,:) = calfactors;

end

calfactors = mean(calfactors_array,1);



