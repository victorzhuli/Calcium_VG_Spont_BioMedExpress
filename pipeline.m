%%%% pipeline for vg analysis
% global parameter:
window_length = 100; window_step = 50;     
threshold_whk = 10; threshold_pup = 10;

for mouse_id = 1:6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % convert calcium data to VG
    exp_session = 'am';
    lz_windowed_VG_conv(mouse_id, exp_session, window_length, window_step);
    
    exp_session = 'pm'; 
    lz_windowed_VG_conv(mouse_id, exp_session, window_length, window_step);
    
    % extract VG metrics
    lz_VG_metrics(mouse_id, window_length, window_step);
    
    % whisker labels
    exp_session = 'am';
    lz_windowed_whisker_conv(mouse_id, exp_session, window_length, window_step, threshold_whk);
%     lz_windowed_pupil_conv(mouse_id, exp_session, window_length, window_step, threshold_pup);A 
    
    exp_session = 'pm';
    lz_windowed_whisker_conv(mouse_id, exp_session, window_length, window_step, threshold_whk);
%     lz_windowed_pupil_conv(mouse_id, exp_session, window_length, window_step, threshold_pup);
    
    % organize data to the .csv format for classification
    %NOTE: before running the following function, all 6 mice should be
    %processed with the designated win_len and win_step
    lz_format4ML(mouse_id, window_length, window_step)
end

%% compute the correlation between behavior and cortex VG metrics 
% produce statisitcal hist for showing consistency across trials and subjects
w_ln = 200; w_st = 50; nSub = 6; nCh = 30;

metric1 = 'whisker'; metric2 = 'pathlength'; 

[pk_coef, pk_idx] = lz_exam_corr_VGmetrics_whk_pup(metric1, metric2, w_ln, w_st, 0); % last parameter: 1 - positive, 0 - negative
iCh = 6;
clc
figure; clf;
fprintf(['===========================\n=== ', metric1, ' ~ ', metric2, ' ===\n']);
for mouse_id = 1: nSub
    subplot(2, nSub, mouse_id); hist(pk_idx(mouse_id,:, iCh), 21);
    xlabel('Lag'); ylabel('Freq.'); title([metric1, ' ~ ', metric2, ', GC6f\_emx\_0',num2str(mouse_id)]);
    subplot(2, nSub, mouse_id+nSub); histfit(pk_coef(mouse_id,:, iCh),10); hold on; vline(0)
    xlabel('Coef.'); ylabel('Freq.');
    fprintf(['----------------------------\nMean peaked coefficient for mouse 0',num2str(mouse_id),': ', num2str(nanmean(pk_coef(mouse_id,:,iCh))),'\n']);
    fprintf(['STD of peaked coefficient for mouse 0',num2str(mouse_id),': ', num2str(nanstd(pk_coef(mouse_id,:,iCh))),'\n']);
end

%% compute the correlation between behavior and cortex VG metrics
metric1 = 'whisker'; metric2 = 'pathlength';
w_ln = 200; w_st = 50; nSub = 6; nCh = 30;
[pk_coef, pk_idx] = lz_exam_corr_VGmetrics_whk_pup(metric1, metric2, w_ln, w_st, 1);
    %%%% pk_coef: #mouse X #trial X #ROI

% project the correlation values on cortex image
pk_coef = abs(pk_coef);
lz_cortex_image_projecting(pk_coef, metric1, metric2)

%% load classification results
%%
lz_load_classification_results('LR');
lz_load_classification_results('RF');
lz_load_classification_results('KNN');
lz_load_classification_results('SVM');
lz_load_classification_results_across_subject('LR');
lz_load_classification_results_across_subject('KNN');
%% present individual classification results using Table and Figure
%%
classifier = 'KNN'; ifplot = 1;
lz_review_classify_result(classifier,ifplot)
%% Find AC and SP based on the value of SE
classifier = 'KNN';
mouse_id = 1;
feature_type = 'D+C';
wLn = 200;
% param = 1;
fprintf('===================\n');
lz_find_classify_results_based_on_condition(classifier, feature_type, wLn, mouse_id);
%% present across-subject classification results using Table and Figure
%%
classifier = 'KNN'; ifplot = 1;
lz_review_classify_result_across_subject(classifier,ifplot)
%% Find AC and SP based on the value of SE for across subject classification
classifier = 'KNN';
mouse_id = 6;
feature_type = 'D+C';
wLn = 200;
% param = 1;
fprintf('===================\n');
lz_find_classify_results_based_on_condition_across_subject(classifier, feature_type, wLn, mouse_id);

%% compare VG based methods with other framework
%% compare with FluoroSNNAP
%%%% organize raw data to .csv files as the FluoroSNNAP requires
for mouse_id =  1:6
    ROI_4_FluoroSNNAP = lz_ROI_FluoroSNNAP(mouse_id);
end
%% detect spikes using FluoroSNNAP (adapted)
thr = .85;
event_amplitude = 5;
parallel = 1; ifplot = 0;

for mouse_id = 2:6
    spks = lz_SpikesLib(mouse_id,'spikes_lz.mat', thr, event_amplitude, parallel, ifplot);
end
%% organize into .csv file for classification by Python
for mouse_id = 1:6
    lz_format4ML_spike(mouse_id)
end
%% examine classification results based on spikes and compare it with that using VG metrics
%%
exam_classification_results



























