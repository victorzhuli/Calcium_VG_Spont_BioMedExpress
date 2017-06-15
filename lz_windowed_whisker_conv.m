function whisker_label = lz_windowed_whisker_conv(mouse_id, exp_session, w_ln_VG, w_st_VG, threshold)

%%% label the whishker movement as the standard deviation of the windowed segments.
%%% INPUT: w_ln_VG and w_st_VG is the number of window used for the calcium data.
%%% INPUT: threshold is used for labeling the whisker data using the standard deviation value.
%%% OUTPUT: whisker_conv_win: continuous output of whisker movement (quantified by windowed standard deviation)
%%% OUTPUT: w_ln and w_st are window parameters that used for calculate std of whisker data, NOT THE SAME AS w_ln_VG and w_st_VG!!!
%%% OUTPUT: whisker_label: threshold labels based on whisker_conv_win and threshold
%%% update the VG_results folder

% mouse_id = 5; exp_session= 'am'; nWn_VG = 35; threshold =7;

nLn = 2000;
nWn_VG = (nLn - w_ln_VG) / w_st_VG +1; % # total window

mouse_name = sprintf('GC6f_emx_%02d',mouse_id);

if (mouse_id == 01 || mouse_id == 02 || mouse_id == 03) &&  strcmp(exp_session, 'am')
    folder_name = '150421am GC6-emx 1-3 spont';
elseif (mouse_id == 01 || mouse_id == 02 || mouse_id == 03) &&  strcmp(exp_session, 'pm')
    folder_name = '150421pm GC6-emx 1-3 spont';
elseif (mouse_id == 04 || mouse_id == 05 || mouse_id == 06) &&  strcmp(exp_session, 'am')
    folder_name = '150609am GC6-emx 4-6 spont';
else
    folder_name = '150609pm GC6-emx 4-6 spont';
end

w_ln = 150 ; % window length
w_st = 25 ; % window step

% windowed convertion
if ispc
    load(['C:\Users\Li_Lab537\Dropbox\GCaMP6f spont and tone reward\',folder_name,'\',mouse_name,'\anglekeeper.mat']);
else
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/',folder_name,'/',mouse_name,'/anglekeeper.mat']);
end

if (mouse_id == 2 && strcmp(exp_session, 'am')) || (mouse_id == 3 && strcmp(exp_session, 'pm')) 
    nTr = 15;
else
    nTr = 16;
end

nLn = 10000; % length of time-series
nWn = (nLn - w_ln) / w_st+1; % # total window

whisker_conv_win = nan(nTr, nWn);
whisker_label = nan(nTr, nWn_VG);
for iTr = 1: nTr
    for iWn = 1:nWn
        widx = [(iWn-1)*w_st+1, (iWn-1)*w_st+w_ln];
        data2compute = anglekeeper(iTr,widx(1):widx(end));
        data2compute = inpaint_nans(data2compute, 5); % interpolate NaN points using 8 neighbor average.
        whisker_conv_win(iTr, iWn) = std(data2compute - mean(data2compute)); % variance
    end
end
% organize to the same window length as calcium-VG data used.
border_ind = linspace(1, nWn, nWn_VG+1);
for iTr = 1: nTr
    for iWn = 1:nWn_VG
        data2compute4label = whisker_conv_win(iTr, ceil(border_ind(iWn)):ceil(border_ind(iWn+1)));
        whisker_label(iTr, iWn) = ( mean(data2compute4label)>threshold );
    end
end
% save results
saveName = sprintf('GC6f_emx_%02d_%s_spont_VG_windowLen%03d_winStep_%03d_whisker_labels_threshold_%02d', mouse_id, exp_session, w_ln_VG, w_st_VG,threshold);
if ispc
    save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\whisker_pupil_auto_label\', saveName], 'whisker_conv_win', 'whisker_label');
else
    save(['/Users/lizhu/Dropbox/projects/calcium/whisker_pupil_auto_label/', saveName], 'whisker_conv_win', 'whisker_label');
end
% % update the new whisker label to the saved dataset
% for iTr = 1:nTr
%     VG_lab_auto = whisker_label(iTr,:);
%     trial = sprintf('%02d',iTr);
%     load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\VG_result\GC6f_emx_0',num2str(mouse_id),'_',exp_session,'_spont_VG_windowLen300_winStep_050_trial',trial]);
%     save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\VG_result\GC6f_emx_0',num2str(mouse_id),'_',exp_session,'_spont_VG_windowLen300_winStep_050_trial',trial],'tidx_ca_di','tidx_whik','VG','VG_lab','VG_lab_auto');
% end
%%%%--------------------------------------------------------------------------
%%%% plot
%%
% iTr = 11;
% figure(1);clf;
% subplot(511); plot(anglekeeper(iTr,1:10000),'r')
%     border_ind = linspace(1, 10000, nWn_VG+1);
%     for iWn = 1: 35
%         vline(border_ind(iWn), 'g')
%     end
% subplot(512); stem(whisker_conv_win(iTr,:)); xlim([0 size(whisker_conv_win,2)]); hold on;
%     border_ind = linspace(1, size(whisker_conv_win,2), nWn_VG+1);
%     for iWn = 1: 35
%         vline(border_ind(iWn), 'g')
%     end
%     hline(threshold)
% subplot(513); stem(whisker_label(iTr,:),'linew',3); ylim([-.5 1.5]); title('auto-labeled')
% trial = sprintf('%02d',iTr);
% load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\VG_result\GC6f_emx_0',num2str(mouse_id),'_',exp_session,'_spont_VG_windowLen300_winStep_050_trial',trial]);
% subplot(514); stem(VG_lab, 'linew',3); ylim([-.5 1.5]); title('hand-labeled')

        
    
    


