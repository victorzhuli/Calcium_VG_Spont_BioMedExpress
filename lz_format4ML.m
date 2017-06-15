%%%% organize the data into single matrix and save for .csv file for classification
function lz_format4ML(mouse_id, window_length, window_step)

nLn = 2000;
w_ln = window_length; % window length
w_st = window_step; % window step
nWn = (nLn - w_ln) / w_st+1; % # total window

%%%%% load labels
% am data -----------------------------------------------------------------
exp_session = 'am';
% whisker label  % added threshold on loadName 03/07/2017
loadName = sprintf('GC6f_emx_%02d_%s_spont_VG_windowLen%03d_winStep_%03d_whisker_labels_threshold_10', mouse_id, exp_session, window_length, window_step);
if ispc
    load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\whisker_pupil_auto_label\', loadName]);
else
    load(['/Users/lizhu/Dropbox/projects/calcium/whisker_pupil_auto_label/', loadName]);
end
VG_whk_cont_am = whisker_conv_win';
VG_whk_lab_am = whisker_label'; % NOTE: here, #window X #trial
% pupil label
loadName = sprintf('GC6f_emx_%02d_%s_spont_VG_windowLen%03d_winStep_%03d_pupil_labels', mouse_id, exp_session, window_length, window_step);
if ispc
    load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\whisker_pupil_auto_label\', loadName]);
else
    load(['/Users/lizhu/Dropbox/projects/calcium/whisker_pupil_auto_label/', loadName]);
end
VG_pup_cont_am = pupil_slope';
VG_pup_lab_am = pupil_label';
% resampling whisker continuous data
for iTr = 1: size(VG_whk_cont_am,2)
    VG_whk_cont_am = resample(VG_whk_cont_am, size(VG_pup_lab_am,1), size(VG_whk_cont_am,1));
end

% pm data -----------------------------------------------------------------
exp_session = 'pm';
% whisker label % added threshold on loadName 03/07/2017
loadName = sprintf('GC6f_emx_%02d_%s_spont_VG_windowLen%03d_winStep_%03d_whisker_labels_threshold_10', mouse_id, exp_session, window_length, window_step);
if ispc
    load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\whisker_pupil_auto_label\', loadName]);
else
    load(['/Users/lizhu/Dropbox/projects/calcium/whisker_pupil_auto_label/', loadName]);
end
VG_whk_cont_pm = whisker_conv_win';
VG_whk_lab_pm = whisker_label'; % NOTE: here, #window X #trial
% pupil label
loadName = sprintf('GC6f_emx_%02d_%s_spont_VG_windowLen%03d_winStep_%03d_pupil_labels', mouse_id, exp_session, window_length, window_step);
if ispc
    load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\whisker_pupil_auto_label\', loadName]);
else
    load(['/Users/lizhu/Dropbox/projects/calcium/whisker_pupil_auto_label/', loadName]);
end
VG_pup_cont_pm = pupil_slope';
VG_pup_lab_pm = pupil_label';
% resampling whisker continuous data
for iTr = 1: size(VG_whk_cont_pm,2)
    VG_whk_cont_pm = resample(VG_whk_cont_pm, size(VG_pup_lab_pm,1), size(VG_whk_cont_pm,1));
end

% combine 'am' and 'pm'
VG_whk_cont = cat(2, VG_whk_cont_am, VG_whk_cont_pm);                        % whisker continuous values: #windows X #trials
whk_cont = reshape(VG_whk_cont, 1, size(VG_whk_cont,1)*size(VG_whk_cont,2)); % concatenated whisker continuous values
VG_whk_lab = cat(2, VG_whk_lab_am, VG_whk_lab_pm);
whk_y = reshape(VG_whk_lab, 1, size(VG_whk_lab,1)*size(VG_whk_lab,2));

VG_pup_cont = cat(2, VG_pup_cont_am, VG_pup_cont_pm);                        % these four lines have the same format as whisker.
pup_cont = reshape(VG_pup_cont, 1, size(VG_pup_cont,1)*size(VG_pup_cont,2));
VG_pup_lab = cat(2, VG_pup_lab_am, VG_pup_lab_pm);
pup_y = reshape(VG_pup_lab, 1, size(VG_pup_lab,1)*size(VG_pup_lab,2));

loadName = sprintf('GC6f_emx_%02d_spont_VG_metric_windowLen%03d_winStep_%03d_v2', mouse_id, w_ln, w_st);
if ispc
    load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\VG_metrics\',loadName])
else
    load(['/Users/lizhu/Dropbox/projects/calcium/VG_metrics/',loadName])
end
VG_degree = reshape(VG_degree, 30, size(VG_degree,2)*size(VG_degree,3));
VG_clustering = reshape(VG_clustering, 30, size(VG_clustering,2)*size(VG_clustering,3));
VG_pathlength = reshape(VG_pathlength, 30, size(VG_pathlength,2)*size(VG_pathlength,3));

% combine them into one matrix
combined_features = [VG_degree', VG_clustering', VG_pathlength', whk_y', whk_cont', pup_y', pup_cont'];
    % degree from 30 columns (channel), clustering from 30  columns (channel), pathlength from
    % 30  columns (channel), label for 1 column

saveName = sprintf('format4ML_GC6f_emx_%02d_windowLen%03d_winStep_%03d_v2_threshold_10.csv', mouse_id, w_ln, w_st);
if ispc
    csvwrite(['C:\Users\Li_Lab537\Dropbox\projects\calcium\format4ML\', saveName], combined_features)
else
    csvwrite(['/Users/lizhu/Dropbox/projects/calcium/format4ML/', saveName], combined_features)
end