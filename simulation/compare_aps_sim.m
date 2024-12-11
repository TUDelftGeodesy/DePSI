function compare_aps_sim(atmo,filenames_aps,Npsc_selections)

global m2ph

Nifgs = size(filenames_aps,1);
M = size(atmo,1);
for w = 1%:Npsc_selections
    for v = 1:Nifgs
        atmo1 = atmo(:,:,1)-atmo(:,:,v+1);
        fid2 = fopen(['/home/freek/projects/ps_software/' ...
                      'test_v1.7.0atmo/arc_network/' filenames_aps(v,:) num2str(w)], ...
                     'r');
        atmo2 = fread(fid2,[M M],'double');
        atmo2 = atmo2/m2ph;
        fid3 = fopen(['/home/freek/projects/ps_software/' ...
                      'test_v1.7.0atmo/star_network/' filenames_aps(v,:) num2str(w)], ...
                     'r');
        atmo3 = fread(fid3,[M M],'double');
        atmo3 = atmo3/m2ph;
        figure;
        subplot(2,3,1);imagesc(atmo1);colorbar
        title('atmo1');
        subplot(2,3,2);imagesc(atmo2);colorbar
        title('atmo2');
        subplot(2,3,3);imagesc(atmo3);colorbar
        title('atmo3');
        subplot(2,3,4);imagesc(atmo1-atmo2);colorbar
        title('atmo1-atmo2');
        subplot(2,3,5);imagesc(atmo1-atmo3);colorbar
        title('atmo1-atmo3');
        subplot(2,3,6);imagesc(atmo2-atmo3);colorbar
        title('atmo2-atmo3');
    end
end
