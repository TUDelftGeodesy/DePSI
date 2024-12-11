function slc_selection = ps_determine_slc_selection(orbitnr,slc_selection_file)

% Function to determine the slc images used for the amplitude
% dispersion based on input file.
%
% Input:    - orbitnr       slc filenames
%           - slc_selection_file  filename selected slc's
%
% Output:   - slc_selection       selected slc's (in order of slc filenames)
%
% ----------------------------------------------------------------------
% File............: ps_determine_slc_selection.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Author..........: Freek van Leijen
%                   Delft Institute of Earth Observation and Space Systems
%                   Delft University of Technology
% ----------------------------------------------------------------------
%
% This software is developed by Delft University of Technology and is
% intended for scientific use only. Applications for commercial use are
% prohibited.
%  
% Copyright (c) 2004-2009 Delft University of Technology, The Netherlands
%
% Change log
%


slaves = cellstr(orbitnr);

good_amplitudes = textread(slc_selection_file,'%s');
slc_index = NaN(1,size(good_amplitudes,1));

for v = 1:size(good_amplitudes,1)
  slc = find(strcmp(slaves,good_amplitudes(v)));
  if ~isempty(slc)
    slc_index(v) = slc;
  end
end

slc_selection = sort(slc_index(~isnan(slc_index)));

