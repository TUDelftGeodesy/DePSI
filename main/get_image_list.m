function [slcs,startDate,stopDate] = get_image_list(processDir,startDate,stopDate,excludeDate,sensor)

% function to get list of slc images to include in the analysis

dirInfo = dir( processDir );
assert( ~isempty( dirInfo ), 'not a directory: %s', processDir );
slcs = regexp( { dirInfo(:).name }, '^(?<slc>[12][0-9]{7})$', 'names' );
slcs = [ slcs{:} ];
slcs = squeeze(struct2cell(slcs));

if ~isempty(startDate)
  startIdx = find(strcmp(startDate,slcs));
  if isempty(startIdx)
    error('You specified a wrong date for the first image.');
  end
else
  startIdx = 1;
  startDate = char(slcs(1));
end

if ~isempty(stopDate)
  stopIdx = find(strcmp(stopDate,slcs));
  if isempty(stopIdx)
    error('You specified a wrong date for the last image.');
  end
else
  stopIdx = size(slcs,1);
  stopDate = char(slcs(end));
end

slcs = slcs(startIdx:stopIdx);

if ~isempty(excludeDate)
  excludeDate = cellstr(excludeDate);
  excludeIdx = NaN(length(excludeDate),1);
  for v = 1:length(excludeDate)
    excludeIdxTemp = find(strcmp(excludeDate(v),slcs));
    if isempty(excludeIdxTemp)
      error(['You want to exclude an image that is not in the current' ...
	     ' selection.']);
    else
      excludeIdx(v) = excludeIdxTemp;
    end
  end
  slcs(excludeIdx) = [];
end


