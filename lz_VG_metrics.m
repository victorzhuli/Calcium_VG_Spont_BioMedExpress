function [VG_degree, VG_clustering, VG_pathlength] = lz_VG_metrics(mouse_id, window_length, window_step)

%%%% extract metrics from visibility graph constructed from raw time series
%%%% output format: num_ROI X num_window X num_trial

nLn = 2000; % length of time-series

w_ln = window_length; % window length
w_st = window_step; % window step
nWn = (nLn - w_ln) / w_st+1; % # total window
nCh = 30;

% am data
if mouse_id == 02
    nTr = 15;
else
    nTr = 16;
end

VG_degree_am = nan(nCh, nWn, nTr);
VG_clustering_am = nan(nCh, nWn, nTr);
VG_pathlength_am = nan(nCh, nWn, nTr);

for iTr = 1:nTr
    loadName = sprintf('GC6f_emx_%02d_am_spont_VG_windowLen%03d_winStep_%03d_trial%02d_v2', mouse_id, w_ln, w_st, iTr);
    if ispc
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\VG_result\',loadName])
    else
        load(['/Users/lizhu/Dropbox/projects/calcium/VG_result/',loadName])
    end
    for iWn = 1:nWn
        parfor iCh = 1:nCh
            VG2comp = triu(VG(:,:,iCh,iWn),1);
            VG_degree_am(iCh,iWn,iTr) = mean(VG2comp(:));
            VG_clustering_am(iCh,iWn,iTr) = mean(clustering_coef_bu(VG(:,:,iCh,iWn)));
            VG_D = distance_bin(VG(:,:,iCh,iWn));
            [VG_pathlength_am(iCh,iWn,iTr)] = charpath(VG_D);
        end
    end
    fprintf('Converted mouse_%02d_am for Trial %02d\n',mouse_id, iTr)
end

% pm data
if mouse_id == 03
    nTr = 15;
else
    nTr = 16;
end
VG_degree_pm = nan(nCh, nWn, nTr);
VG_clustering_pm = nan(nCh, nWn, nTr);
VG_pathlength_pm = nan(nCh, nWn, nTr);

for iTr = 1:nTr
    loadName = sprintf('GC6f_emx_%02d_pm_spont_VG_windowLen%03d_winStep_%03d_trial%02d_v2', mouse_id, w_ln, w_st, iTr);
    if ispc
        load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\VG_result\',loadName])
    else
        load(['/Users/lizhu/Dropbox/projects/calcium/VG_result/',loadName])
    end
    for iWn = 1:nWn
        parfor iCh = 1:nCh
            VG2comp = triu(VG(:,:,iCh,iWn),1);
            VG_degree_pm(iCh,iWn,iTr) = mean(VG2comp(:));
            VG_clustering_pm(iCh,iWn,iTr) = mean(clustering_coef_bu(VG(:,:,iCh,iWn)));
            VG_D = distance_bin(VG(:,:,iCh,iWn));
            [VG_pathlength_pm(iCh,iWn,iTr)] = charpath(VG_D);
        end
    end
    fprintf('Converted mouse_%02d_pm for Trial %02d\n',mouse_id,iTr)
end

VG_degree = cat(3,VG_degree_am, VG_degree_pm);
VG_clustering = cat(3,VG_clustering_am, VG_clustering_pm);
VG_pathlength = cat(3,VG_pathlength_am, VG_pathlength_pm);

saveName = sprintf('GC6f_emx_%02d_spont_VG_metric_windowLen%03d_winStep_%03d_v2', mouse_id, w_ln, w_st);
if ispc
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\VG_metrics\',saveName], 'VG_degree', 'VG_clustering', 'VG_pathlength');
else
    save(['/Users/lizhu/Dropbox/projects/calcium/VG_metrics/',saveName], 'VG_degree', 'VG_clustering', 'VG_pathlength');
end