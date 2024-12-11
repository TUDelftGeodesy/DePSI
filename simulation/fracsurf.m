function [fsurf] = fracsurf(N,beta,doplot,P0);
% [fsurf] = fracsurf(N,beta,[doplot,P0]);
%
%  Script to create an isotropic 2 dimensional fractal surface 
%  with a power law behavior
%
%  beta : power law exponent for a 1D profile of the data
%         1D:  S(k) = C0 k^(-beta)
%         2D:  S(k) = C1 k^(-(beta+1))
%
%  N    : size of the square area
%  [doplot] : 'y'/'n' show figure (Default 'n')
%  [P0]     : multiply AMPLITUDE (!) spectrum by sqrt(P0) (default P0 = 1)
%
%  fsurf: real-valued 2D output surface
% 
% e.q., [fsurf] = fracsurf(128, 1.3 ,'y',1000);
%
% Ramon Hanssen, April 2000.
% RH 24-May-2000 14:25 : Added scaling with P0

if nargin==0,help fracsurf;return;end
if nargin<4, P0 = sqrt(1)  ;end  % NOTE: for power we use P(k) = P0 x k^(-beta)
                                 %    Here we use amplitude only, so sqrt(P0) in stead of P0
if nargin<3, doplot = 'n'  ;end

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
  noemer = k .^(beta/2);       % The power beta/2 is used because the power spectral
                               % density is proportional to the amplitude squared 
                               % Here we work with the amplitude, instead of the power
                               % so we should take sqrt( k.^beta) = k.^(beta/2)  RH

% prevent dividing by zero;
  noemer(find(noemer==0))=1;

% weight with "noemer" and with sqrt(P0) (sqrt because the operation is on the amplitude
%                                         spectrum instead of the power spectrum)
  Hnew   = sqrt(P0) .* H ./ noemer;

% Create spectral surface by ifft
  fsurf = abs(ifft2(Hnew));
% Remove mean to get zero-mean data
  fsurf = fsurf - mean(fsurf(:));


  if strcmp(doplot,'y'), 
    figure(1);imagesc(fsurf);colormap(gray);colorbar
  end

