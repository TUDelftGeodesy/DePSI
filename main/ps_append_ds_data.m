function Npsp = ps_append_ds_data(dsFileName,Npsp,Npsc_selections)

% Function to add DS to the PSP selection
%
% Input:  - dsFileNa           filename DS data file (.mat)
%         - Npsp               number of PSP
%         - Npsc_selections    number of PSC selections
%
% Output: - Npsp               updated number of PSP
%
% ----------------------------------------------------------------------
% File............: ps_append_ds_data.m
% Version & Date..: 23-SEP-2021
% Author..........: Philip Conroy
%                   Delft University of Technology
% ----------------------------------------------------------------------
%
% This software is developed by Delft University of Technology and is
% intended for scientific use only. Applications for commercial use are
% prohibited.
%  
% Copyright (c) 2004-2021 Delft University of Technology, The Netherlands
% 
% Change log
%


% ----------------------------------------------------------------------
% Initialize
% ----------------------------------------------------------------------

global project_id

% Load multilooked DS Data (virtual PS)
exportStruct = load(dsFileName); 
dsData = exportStruct.exportData;
nDs = size(dsData,1);

for z = 1:Npsc_selections
    psp_fid = fopen([project_id '_psp_2orig_sel' num2str(z) '.raw'],'a');
    fwrite(psp_fid,dsData','double');
    fclose(psp_fid);
    
    % Update number of secondary PSC with the additional DS "points"
    Npsp(z) = Npsp(z) + nDs;
end

end

