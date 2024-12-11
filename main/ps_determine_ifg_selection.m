function ifg_selection = ps_determine_ifg_selection(orbitnr,ifg_selection_file)

% Function to determine the slc images used for the amplitude
% dispersion based on input file.
%
% Input:    - orbitnr       slc filenames
%           - ifg_selection_file  filename selected ifg's
%
% Output:   - ifg_selection       selected ifg's (in correct order)
%
% ----------------------------------------------------------------------
% File............: ps_determine_ifg_selection.m
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

selected_ifg = textread(ifg_selection_file,'%s');
ifg_index = NaN(1,size(selected_ifg,1));

for v = 1:size(selected_ifg,1)
  ifg = find(strcmp(slaves,selected_ifg(v)));
  if ~isempty(ifg)
    ifg_index(v) = ifg;
  end
end

ifg_selection = sort(ifg_index(~isnan(ifg_index)));

