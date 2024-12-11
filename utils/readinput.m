function value = readinput(parameter,param_file_content,param_file)

% readinput
%
% ----------------------------------------------------------------------
% File............: readinput.m
% Version & Date..: 1.7.2, 21-JAN-2007
% Authors.........: Freek van Leijen
%                   Delft Institute of Earth Observation and Space Systems
%                   Delft University of Technology
% ----------------------------------------------------------------------
%
% This software is developed by Delft University of Technology and is
% intended for scientific use only. Applications for commercial use are
% prohibited.
%  
% Copyright (c) 2004-2006 Delft University of Technology, The Netherlands
%

if ~isempty(param_file)
  param_file_content = textread(param_file,'%s');
end

index = strmatch(parameter,param_file_content,'exact');
Nindex = length(index);

if Nindex == 0
  error(['Parameter ' parameter ' not found in parameter file.']);
elseif Nindex == 1
  value_temp = char(param_file_content(index+2));
elseif Nindex > 1
  error(['Multiple occurrences of parameter ' parameter ' found.']);
end

[value,valid] = str2num(value_temp);
if valid==0
  value = value_temp;
  if strcmp(value(1),'''');
    value = value(2:end);
  end
  if strcmp(value(end),'''')
    value = value(1:end-1);
  end
end

