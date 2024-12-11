function [Ngroup, Np_group, arc_group, N_arc_avg] = ps_spatial_unwrapping_grouping(A)

% The function is used to group arcs, i.e. subdivide all arcs to a number
% of groups. The number in each group should not larger than the number of
% number of pscs.
%
% Input:    - A          design matrix which is sparse
%
% Output:   - Ngroup     number of groups in total
%           - Np_group   number of points per group
%           - arc_group  number of arcs per group
%           - N_arc_avg  average number of arcs per group
%
% ----------------------------------------------------------------------
% File............: ps_spatial_unwrapping_grouping.m
% Version & Date..: 1.7.2.16, 12-DEC-2009
% Author..........: Shizhuo Liu
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

global Nworkers unwrapping_mode
unwrapping_mode = 'sequential';
Nworkers=[];

if (numel(A) <= 3*10^4) ;   % A*Rinv is small, memory can hold
                            % the number is user dependent, based on
                            % experience, set it to smaller number might
                            % result in slowing down the process, set it
                            % to a too large number might run out of memory

  Ngroup = 1 ; 
  Np_group = size(A,2) ; % all pscs
  arc_group{1} = 1 : size(A,1) ; % all arcs
  N_arc_avg = [] ;
  
  fprintf(1,'The network has %4.0f arcs and %4.0f points \n',size(A,1), size(A,2));
  fprintf(1, 'No need for grouping, solving as normal \n') ;
  
  return
  
  
elseif(size(A,2) > 10^4)
    
  fprintf(1,['There are in total ', num2str(size(A,2)), ' pscs in the network',...
             'process might crash due to memory limit, Good luck!' ]) ;
  
else
  
  fprintf(1,'The network has %4.0f arcs and %4.0f points \n',size(A,1), size(A,2));
  fprintf(1, 'Needs grouping \n') ;
  fprintf(1, 'Grouping starts ...\n') ;
  
end

continue_flag = 1 ;

Np_total = size(A,2) ; % total number of points = Npscs
                       % A is transposed
                       
if(strcmp(unwrapping_mode, 'sequential'))
  
  Ngroup =  ceil(size(A,1)/Np_total) ; % initial
  
else
  
  Ngroup =  Nworkers ;
  
end

while(continue_flag)
  
  Np_group = floor(Np_total/Ngroup) ;  % number of pscs in the group
  
  for n = 1 : Ngroup
    
    if (n < Ngroup) 
      
      [I,J] =  find(A(:,(n-1)*Np_group+1:n*Np_group)~=0) ;
      arc_list = I' ;
      
    else
      
      [I,J] =  find(A(:,(n-1)*Np_group+1:Np_total)~=0) ; 
      arc_list = I' ;             
    end
    
    
    % arc_group_all{n} = arc_list ;     
    arc_group{n} = unique(arc_list) ; % commom arc_index in the group
				      % avoid multiple times
				      % computation with the normal
				      % matrix later                                   
    Narc_group(n) = length(arc_group{n}) ;
    
  end
    
  if(max(Narc_group) <= Np_total)
    
    continue_flag = 0 ;
    
  else
    
    Ngroup = Ngroup + 1 ;  
    
  end
  
end
N_arc_avg = mean(Narc_group) ;  


