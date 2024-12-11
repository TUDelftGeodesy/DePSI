function [cropIn,cropFinal,crop,nLines,nPixels,Bdop] = get_stack_parameters(processDir,slcs,crop,altimg,sensor,processor,master,swath_burst)

nSlc = size(slcs,1);
cropIn = NaN(nSlc,4);
Bdop = NaN(nSlc,1);
%Bperp = NaN(nSlc,1);
%Btemp = NaN(nSlc,1);

if strcmp(sensor,'s1')
  if strcmp(processor,'doris_rippl')
    if isempty(swath_burst) % merged image
      fileName = fullfile(processDir,master,'info.json');
    else % single burst
      fileName = fullfile(processDir,master,['slice_' num2str(swath_burst(2)) '_swath_' num2str(swath_burst(1))],'info.json');
    end
    masterData = jsondecode(fileread(fileName));
    masterfDC = masterData.readfiles.original.Xtrack_f_DC_constant_Hz_EarlyEdge_;
  elseif strcmp(processor,'doris_flinsar')
    if isempty(swath_burst) % merged image
      fileName = fullfile(processDir,master,'master.res');
    else % single burst
      fileName = fullfile(processDir,master,['swath_' num2str(swath_burst(1))],['burst_' num2str(swath_burst(2))],'slave.res');
    end
    masterData = textread(fileName,'%s','delimiter',':');
    masterfDC = get_parameter('Xtrack_f_DC_constant (Hz, early edge)',masterData,1);
  else
    if isempty(swath_burst) % merged image
      fileName = fullfile(processDir,master,'slave.res');
    else % single burst
      fileName = fullfile(processDir,master,['swath_' num2str(swath_burst(1))],['burst_' num2str(swath_burst(2))],'slave.res');
    end
    masterData = textread(fileName,'%s','delimiter',':');
    masterfDC = get_parameter('Xtrack_f_DC_constant (Hz, early edge)',masterData,1);
  end
else
  masterData = textread(fullfile(processDir,'master.res'),'%s','delimiter',':');
  masterfDC = get_parameter('Xtrack_f_DC_constant (Hz, early edge)',masterData,1);
end

for v = 1:nSlc

  if strcmp(processor,'doris_rippl')
    if isempty(swath_burst) % merged image
      fileName = fullfile(processDir,char(slcs(v)),'info.json');
    else % single burst
      fileName = fullfile(processDir,char(slcs(v)),['slice_' num2str(swath_burst(2)) '_swath_' num2str(swath_burst(1))],'info.json');
    end

    slaveData = jsondecode(fileread(fileName));

    % fDC
    slavefDC = slaveData.readfiles.original.Xtrack_f_DC_constant_Hz_EarlyEdge_;

    %if strcmp(slcs(v),master)
    %  % First_line:
    %  l0 = slaveData.processes.correct_phases.correct_phases__coor__radar__pol__VV.coordinates.first_line + 1; % +1 for python to matlab

    %  % Last_line:
    %  lN = l0 + slaveData.processes.correct_phases.correct_phases__coor__radar__pol__VV.coordinates.shape(1) - 1;

    %  % First_pixel:
    %  p0 = slaveData.processes.correct_phases.correct_phases__coor__radar__pol__VV.coordinates.first_pixel + 1; % +1 for python to matlab

    %  % Last_pixel:
    %  pN = p0 + slaveData.processes.correct_phases.correct_phases__coor__radar__pol__VV.coordinates.shape(2) - 1;
    %else
      % First_line (w.r.t. original_image):  % +1 for python to matlab
      l0 = eval(['slaveData.processes.correct_phases.correct_phases__coor__radar__pol__' altimg '.coordinates.first_line + 1']);

      % Last_line (w.r.t. original_image):
      lN = eval(['l0 + slaveData.processes.correct_phases.correct_phases__coor__radar__pol__' altimg '.coordinates.shape(1) - 1']);

      % First_pixel (w.r.t. original_image):  % +1 for python to matlab
      p0 = eval(['slaveData.processes.correct_phases.correct_phases__coor__radar__pol__' altimg '.coordinates.first_pixel + 1']);

      % Last_pixel (w.r.t. original_image):
      pN = eval(['p0 + slaveData.processes.correct_phases.correct_phases__coor__radar__pol__' altimg '.coordinates.shape(2) - 1']);
    %end

  else

    if strcmp(sensor,'s1')
      if isempty(swath_burst) % merged image
	if strcmp(processor,'doris_flinsar')
	  if strcmp(master,slcs(v))
	    fileName = fullfile(processDir,char(slcs(v)),'master.res');
	  else
	    fileName = fullfile(processDir,char(slcs(v)),'slave.res');
	  end
          line_fileName = fullfile(processDir,'nlines_crp.txt');
          pixel_fileName = fullfile(processDir,'npixels_crp.txt');
	else
          fileName = fullfile(processDir,char(slcs(v)),'slave.res');
        end
      else % single burst
        fileName = fullfile(processDir,char(slcs(v)),['swath_' num2str(swath_burst(1))],['burst_' num2str(swath_burst(2))],'slave.res');
      end
      slaveData = textread(fileName,'%s','delimiter',':');
    else
      slaveData = textread(fullfile(processDir,char(slcs(v)),['slave' altimg '.res']),'%s','delimiter',':');
    end

    % fDC
    slavefDC = get_parameter('Xtrack_f_DC_constant (Hz, early edge)',slaveData,1);

    if strcmp(processor,'doris_flinsar')
      lineID = fopen(line_fileName,'r');
      pixID = fopen(pixel_fileName,'r');
      line_info = fscanf(lineID,'%d');
      pix_info = fscanf(pixID,'%d');      
      fclose(lineID);
      fclose(pixID);
      l0 = line_info(2);
      lN = line_info(3);
      p0 = pix_info(2);
      pN = pix_info(3);
      %% First_line (w.r.t. original_image):
      %l0 = get_parameter('First_line (w.r.t. original_image)',slaveData,1);
    
      %% Last_line (w.r.t. original_image):
      %lN = get_parameter('Last_line (w.r.t. original_image)',slaveData,1);
    
      %% First_pixel (w.r.t. original_image):
      %p0 = get_parameter('First_pixel (w.r.t. original_image)',slaveData,1);
    
      %% Last_pixel (w.r.t. original_image):
      %pN = get_parameter('Last_pixel (w.r.t. original_image)',slaveData,1);
    else  
      % First_line (w.r.t. original_image):
      l0 = get_parameter('First_line (w.r.t. original_master)',slaveData,1);

      % Last_line (w.r.t. original_image):
      lN = get_parameter('Last_line (w.r.t. original_master)',slaveData,1);

      % First_pixel (w.r.t. original_image):
      p0 = get_parameter('First_pixel (w.r.t. original_master)',slaveData,1);

      % Last_pixel (w.r.t. original_image):
      pN = get_parameter('Last_pixel (w.r.t. original_master)',slaveData,1);
    end
  end

  cropIn(v,:) = [l0 lN p0 pN];
  
  Bdop(v) = masterfDC - slavefDC;

  %% Bperp
  %Bperp(v) = get_parameter('Bperp',ifgsData,1);

  %% Btemp
  %Btemp(v) = get_parameter('Btemp',ifgsData,2);

end

cropFinal = [max(cropIn(:,1)) min(cropIn(:,2)) max(cropIn(:,3)) min(cropIn(:,4))];

if ~isempty(crop)
  crop = [max(crop(1),cropFinal(1)) min(crop(2),cropFinal(2)) ...
          max(crop(3),cropFinal(3)) min(crop(4),cropFinal(4))];
else
  crop = cropFinal;
end
nLines = crop(2)-crop(1)+1;
nPixels = crop(4)-crop(3)+1;


