function [K0,C0,coh0,phase_residual]=ps_topofit(cpxphase,bperp,n_trial_wraps,plot_label)
%PS_TOPOFIT find best-fitting range error 
% mcc if you want to use it to detect psc then:
% cpxphase: cpx phase - weigthed averaged
% bperp : perpendicualr baseline or h2ph
% n_trial_wraps : maximum number of cycles you expect the phase will be wrapped with the Bperp
% if zero it will consider linear relationship between wrapped phase and Bperp.
% if one it will consider linear relationship but one jump can happenp, and so on.
% [plot_label] used for outputing the plots 
% 
%   PS_TOPOFIT(cpxphase,bperp,n_trial_wraps,plotflag)
%
% OUTPUT= 
% [K0,C0,coh0,phase_residual]
% K0  The model considers angle(cpxphase) = K0 *Bperp + C0
% k0 is in rad. If you use Btemp instead of Bperp then it returns the deformation rate in rad/year
% if bperp = h2ph then K0 is the dem error, other wise K0 should be multiplied times Rp*sin(alpha)
% It estimates different K and selects the one that maximisez the cohenrence. Since infinite values
% cannot be used, it tries from  -8*n_trial_wraps to 8*n_trial_wraps (in meters). Then it finds the maximum of the cohe
% and by linear reggresion (using a the coherence as a function) gets a more precise value of where the coherece is maximum 
% , which is K0 
% C0 is the average of the residuals, therefore the noise of the Master
% coh0 is the coherence defined as abs( sum( phase_residual ) )/sum(abs(phase_residual));
%
%   Andy Hooper, June 2006
%
%   ==========================================================
%   04/2007 AH: Added 64-bit machine compatibility
%   04/2007 AH: Tightened up max topo error processing
%   ==========================================================


if nargin < 4 || isempty(plot_label)
 
  output_label=[];
  plotflag = 'n';

else 
plotflag = 'y';

end



if size(cpxphase,2)>1
   cpxphase=cpxphase.';
end

%ix gives the indeces where the signal is different than zero
%ix is zero when cpxphase else ix =1
ix=cpxphase~=0;  % if signal of one image is 0, dph set to 0

%removing values where ix=0, (i.e. cpxphase is zero)
cpxphase=cpxphase(ix);%removes values where ix=0
bperp=bperp(ix);
n_ix=length(ix);

bperp_range=max(bperp)-min(bperp);

%wrapped phases
wphase=angle(cpxphase);

%First uses periodogram
%then applies BLUE to the residuals and estimate the dem error
%total DEM error is teh sum of both

%%
%periodogram starts

%search space
trial_mult=[-ceil(8*n_trial_wraps):ceil(8*n_trial_wraps)];
n_trials=length(trial_mult);
%trial_phase=bperp/maxbperp*pi/4;
%this reads from left to right one by one
trial_phase=bperp/bperp_range*pi/4;
trial_phase_mat=exp(-j*trial_phase*trial_mult);
cpxphase_mat=repmat(cpxphase,1,n_trials);
phaser=trial_phase_mat.*cpxphase_mat;
phaser_sum=sum(phaser);
C_trial=angle(phaser_sum);
%
coh_trial=abs(phaser_sum)/sum(abs(cpxphase));

coh_diff=diff(coh_trial);
%coh_max_ix=1;                       % include 1st value in case minimum is actually off to left
coh_max_ix=[];
for i=2:length(coh_diff)
    if coh_diff(i)<0 & coh_diff(i-1)>0
        coh_max_ix=[coh_max_ix,i];
    end
end 
%coh_max_ix=[coh_max_ix,length(coh_trial)]; % include last value in case minimum is actually off to right     

coh_max=coh_trial(coh_max_ix); % maximum value of coherence


[dummy,coh_high_max_ix]=max(coh_trial); % only select highest

%DEM error values corresponding to the periodogram
%K0=pi/4/maxbperp*trial_mult(coh_high_max_ix);
K0=pi/4/bperp_range*trial_mult(coh_high_max_ix);
C0=C_trial(coh_high_max_ix);
coh0=coh_trial(coh_high_max_ix);


%BLUE starts
% linearise and solve
resphase=cpxphase.*exp(-j*(K0*bperp)); % subtract approximate fit
offset_phase=sum(resphase);
resphase=angle(resphase*conj(offset_phase)); % subtract offset, take angle (unweighted)

weighting=abs(cpxphase); %weihgts on the amplitude and therefore with the product of sin(Bperp)*sinc(Bdop)
%keyboard
mopt=double(weighting.*bperp)\double(weighting.*resphase);%remaining DEM error

K0=K0+mopt;%Total DEM_error

phase_residual=cpxphase.*exp(-j*(K0*bperp)); 

mean_phase_residual=sum(phase_residual); 
C0=angle(mean_phase_residual);    % static offset (due to noise of master + average noise of rest)
coh0=abs(mean_phase_residual)/sum(abs(phase_residual)); 

if plotflag=='y'
%		set(gcf,'visible','off')
    figure(999);
		set(gcf,'visible','off')
    subplot(2,1,2)
    bvec=linspace(min(bperp),max(bperp),200);
	wphase_hat=angle(exp(j*(K0(1)*bvec+C0(1))));
	p=plot(bvec,(wphase_hat),'r');
	hold on
    set(p,'linewidth',2)
    p=plot(bperp,wphase,'bo');
	grid on
    set(p,'linewidth',2)
	hold off
    set(gca,'ylim',[-pi,pi])
    set(gca,'fontsize',12,'fontweight','bold')
    ylabel('Wrapped Phase')
    %ylabel('\psi - \tilde{\psi} (rads)')
    xlabel('B_{\perp} (m)')
    %xlabel('Perpendicular baseline (m)')
    %title('DEM error')
    %legend('Data',['best fit, slope=',num2str(K0(1))],'location','best')
    subplot(2,1,1)
   % plot(pi/4/maxbperp*trial_mult,coh_trial)
		plot(pi/4/bperp_range/4/pi*0.05656*trial_mult,coh_trial,'g')
    grid on
    ylabel('\gamma_x')
    xlabel('\Delta \theta^{nc}_x (dem error) ')
    set(gca,'fontsize',12,'fontweight','bold')
    axis tight
    set(gca,'ylim',[0,1])
    print('-depsc', [ plot_label  '.eps']);
    close(999)
end

