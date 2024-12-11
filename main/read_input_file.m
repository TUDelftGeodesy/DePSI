  [orbitnr,dates,filenames_slc,filenames_ifgs,filenames_h2ph, ...
   filenames_output,Bperp,Btemp,Bdop,l0,lN,p0,pN] = ...
    textread(input_file,'%s %s %s %s %s %s %f %f %f %f %f %f %f');
  crop_in = [l0 lN p0 pN];
  orbitnr = strvcat(orbitnr);
  dates = strvcat(dates);
  
  if strcmp(filenames_ifgs(end),'dummy')
    filenames_ifgs(end) = [];
    filenames_h2ph(end) = [];
    filenames_output(end) = [];
    Bperp(end) = [];
    Btemp(end) = [];
    Bdop(end) = [];
  end
  
  Nslc = size(filenames_slc,1);
  Nifgs = size(filenames_ifgs,1);
  
  [Btemp,Btemp_index] = sort(Btemp);
  orbitnr(1:Nifgs,:) = orbitnr(Btemp_index,:);
  dates(1:Nifgs,:) = dates(Btemp_index,:);
  filenames_slc(1:Nifgs) = filenames_slc(Btemp_index);
  filenames_ifgs = filenames_ifgs(Btemp_index);
  filenames_h2ph = filenames_h2ph(Btemp_index);
  filenames_output = filenames_output(Btemp_index);
  Bperp = Bperp(Btemp_index);
  Bdop = Bdop(Btemp_index);
  crop_in(1:Nifgs,:) = crop_in(Btemp_index,:);
  
  %find master index
  master_index = find(Btemp<0);
  master_index = master_index(end)+1;

  crop_final = [max(crop_in(:,1)) min(crop_in(:,2)) max(crop_in(:,3)) min(crop_in(:,4))];

  if ~isempty(crop)
    crop = [max(crop(1),crop_final(1)) min(crop(2),crop_final(2)) ...
            max(crop(3),crop_final(3)) min(crop(4),crop_final(4))];
  else
    crop = crop_final;
  end
  Nlines = crop(2)-crop(1)+1;
  Npixels = crop(4)-crop(3)+1;


