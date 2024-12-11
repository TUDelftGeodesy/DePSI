function [fsurf] = fracsurfatmo(N,beta,regime,doplot,P0);
% [fsurf] = fracsurfatmo(N,beta,regime, [doplot,P0]);
%
%  Script to create an isotropic 2 dimensional fractal surface 
%  with a power law behavior which corresponds with the 
%  [-2/3, -8/3, -5/3] power law
%
%  beta : power law exponents for a 1D profile of the data
%         1D:  S(k) = C0 k^(-beta)
%         2D:  S(k) = C1 k^(-(beta+1))
%  regimes: cumulative percentage of spectrum covered by a specific beta
%
%  N    : size of the square area
%  [doplot] : 'y'/'n' show figure (default 'y')
%  [P0]     : multiply amplitude spectrum by P0 (default P0 = 1)
%
%  fsurf: real-valued 2D output surface
%        
% e.q.,  beta =  [5/3, 8/3, 2/3]; regime = [95, 99,100]; N=128; doplot='y'
% e.q.,  beta =  [5/3, 8/3, 2/3]; regime = [60, 90,100]; N=128; doplot='y'
%      [fsurf] = fracsurfatmo(N,beta,regime,doplot);
%      [fsurf] = fracsurfatmo(128, [5/3, 8/3, 2/3], [95, 99,100],'y');
%                                                        95% 4% 1%
%
% Ramon Hanssen, May 2000
% RH 24-May-2000 14:24 : Improved scaling with P0


if nargin==0,help fracsurfatmo;return;end
if nargin==3, doplot = 'n'          ;end
if nargin~=5,  
  P0 = 1;  % NOTE: for power we use P(k) = P0^2 x k^(-beta)
           %       Here we use amplitude only, so P0 in stead of P0^2
end

% Check if N is even
if rem(N,2), fprintf(1,'Error: use an even value for N\n');return;end

% Simulate a uniform random signal

  h      = rand(N);
  H      = fftshift(fft2(h));

% scale the spectrum with the power law
  x      = [-N/2: (N/2 -1)];
  [X,Y]  = meshgrid(x,x);
  k      = sqrt(X.^2 + Y.^2);
  beta   = beta+1;             % beta+1 is used as beta, since, the power exponent
                               % is defined for a 1D slice of the 2D spectrum:
                               % austin94: "Adler, 1981, shows that the surface profile 
                               %   created by the intersection of a plane and a
                               %   2-D fractal surface is itself fractal with 
                               %   a fractal dimension  equal to that of the 2D 
                               %   surface decreased by one."
  beta = beta./2;             % The power beta/2 is used because the power spectral
                               % density is proportional to the amplitude squared 
                               % Here we work with the amplitude, instead of the power
                               % so we should take sqrt( k.^beta) = k.^(beta/2)  RH
  noemer = zeros(size(k));
  mk      = max(X(:));
  k0      = 0;
  k1      = (regime(1)/100)*mk;
  k2      = (regime(2)/100)*mk;
  k3      = max(k(:));
  regime1 = find( k >  k0 & k <= k1 );
  regime2 = find( k >=  k1 & k <= k2 );
  regime3 = find( k >=  k2 & k <= k3 );

  noemer(regime1) = (k(regime1)   ).^beta(1);
  noemer(regime2) = (k(regime2)).^beta(2) ./ min((k(regime2)).^beta(2)) .* max(noemer(regime1));
  noemer(regime3) = (k(regime3)).^beta(3) ./ min((k(regime3)).^beta(3)) .* max(noemer(regime2));

aa=9; %TEST
if aa==2,figure(100);loglog(k(regime1),(k(regime1)   ).^beta(1),'x');
         figure(100);hold on;loglog(...
          k(regime2),...
         (k(regime2)).^beta(2) ./ min((k(regime2)).^beta(2)) .* max(noemer(regime1))...
       ,'cx');hold off;
         figure(100);hold on;loglog(...
           k(regime3),...
           (k(regime3)).^beta(3) ./ min((k(regime3)).^beta(3)) .* max(noemer(regime2))...
        ,'yx');hold off;
end;

% prevent dividing by zero;
  noemer(find(noemer==0))=1;

 
  Hnew   = P0 .* H ./noemer;

% Create spectral surface by ifft
  fsurf = abs(ifft2(Hnew));
% Remove mean to get zero-mean data
  fsurf = fsurf - mean(fsurf(:));

  if strcmp(doplot,'y'), 
    figure(1);imagesc(fsurf);colorbar
  end

