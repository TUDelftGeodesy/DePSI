function [Bn, err, errBn] = baseline(N)

% Generate random orthogonal baseline values exspressed in methers 
%
% [Bn, err, errBn] = baseline(N, mean);
%
% INPUT
%
% N        number of requested orthogonal B
%
% OUTPUT
%
% Bn       vector [Nx1] cointaining baseline values
% err      error in orbit determination
% errBn    baseline values affected by orbit determination error

% baseline generation 3sigma= 500 m
Bn= randn(N,1)*500/1.5;


err = (rand(N,1) -0.5)*0.4; % accuracy in absolute orbit determination
                            % is 5-10 cm for each orbit. The master orbit is considered
                            % deterministic and the error is all related to the slave one.
                            % This brings the error in the interval [-0.1, 0.1] methers
errBn = Bn + err;



