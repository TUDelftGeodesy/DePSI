function [acqdates,Btemp] = ERS2_time_series(N)

%if N>60
%  error('Maxmimum 60 ERS2 images available')
%end

%ERS2_start = datenum('21-Mar-95');
%ERS2_stop = datenum('15-Jan-2001');
ERS2_start = datenum('21-Mar-90');
ERS2_stop = datenum('15-Jan-2015');

ERS2_range = ERS2_start:ERS2_stop;

RI = 35;  % repeat interval

refdate = ERS2_range(ceil(rand*length(ERS2_range)));

%ERS2_images = sort([refdate+RI:RI:ERS2_stop,refdate:-RI:ERS2_start])';
ERS2_images = sort([refdate+RI:RI:ERS2_stop,refdate-RI:-RI:ERS2_start])';
N_ERS2 = length(ERS2_images);

ERS2_idx = randperm(N_ERS2);
ERS2_idx = sort(ERS2_idx(1:N));

acqdates = ERS2_images(ERS2_idx);
Btemp = acqdates - refdate;


