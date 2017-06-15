function VG = lz_windowed_VG_conv(mouse_id, exp_session, window_length, window_step)

%%%% load data from dropbox --> load label --> segment the data and
%%%% construct VG for each segment
%%%% developed from my script "windowed_VG_conversion.m"
%%%% input: mouse_id: e.g., 02
%%%% input: exp_session: experimental session: am, pm
%%%% input: window_length: in time point, e.g., 300
%%%% input: window_length: in time point, e.g., 300

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

w_ln = window_length; % window length
w_st = window_step; % window step

% load raw data
if ispc
    load(['C:\Users\Li_Lab537\Dropbox\GCaMP6f spont and tone reward\',folder_name,'\',mouse_name,'\Ca.mat']);
    load(['C:\Users\Li_Lab537\Dropbox\GCaMP6f spont and tone reward\',folder_name,'\',mouse_name,'\diameter.mat']);
    load(['C:\Users\Li_Lab537\Dropbox\GCaMP6f spont and tone reward\',folder_name,'\',mouse_name,'\anglekeeper.mat']);
    load(['C:\Users\Li_Lab537\Dropbox\projects\calcium\wSeg\wSeg_',mouse_name,'_',exp_session,'.mat']); % whisker segments
else
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/',folder_name,'/',mouse_name,'/Ca.mat']);
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/',folder_name,'/',mouse_name,'/diameter.mat']);
    load(['/Users/lizhu/Dropbox/GCaMP6f spont and tone reward/',folder_name,'/',mouse_name,'/anglekeeper.mat']);
    load(['/Users/lizhu/Dropbox/projects/calcium/wSeg/wSeg_',mouse_name,'_',exp_session,'.mat']); % whisker segments
end
    
% Calcium dataset
% change cell to matrix: channel X time X trial
Cal = reshape(cell2mat(Ca.Ch0), 30, 2047, size(Ca.Ch0,2));
if (mouse_id == 03 && strcmp(exp_session, 'pm'))
    Cal = Cal(:, 1:2000, [1:11,13:16]);
else
    Cal = Cal(:,1:2000,:);
end
clear diamKeeper anglekeeper

tidx_ca_di = linspace(0,2000/Ca.sample_rate,2000);
tidx_whik  = linspace(0,2000*5/500,2000*5);

nLn = size(Cal,2); % length of time-series
nCh = size(Cal,1); % # of channel
nTr = size(Cal,3); % # trial
nWn = (nLn - w_ln) / w_st+1; % # total window
srate = Ca.sample_rate;

for iTr = 1: nTr
    VG = nan(w_ln, w_ln, nCh, nWn);
    % prepare for labeling for this window
    VG_lab = zeros(1,nWn);
    seg = wSeg{1,iTr}; nSeg = length(seg)/2;
    segseg = nan(nSeg, 2); % #WMsegment X 2 {(start, end) in second.}
    for iSeg = 1:nSeg
        segseg(iSeg,1) = seg(iSeg*2-1); segseg(iSeg,2) = seg(iSeg*2);
    end
    parfor iWn = 1: nWn
        widx = [(iWn-1)*w_st+1, (iWn-1)*w_st+w_ln];
        widx_s = widx/srate; % for labeling
        VG(:,:,:,iWn) = lz_VG_build_2(Cal(:,widx(1):widx(end),iTr)');  %#ok<PFBNS>
        for iSeg = 1:nSeg
            if widx_s(1) >= segseg(iSeg,1) && widx_s(1) <= segseg(iSeg,2)  %#ok<PFBNS>
                VG_lab(iWn) = 1;
            end
        end
        fprintf('Converted mouse %02d %2s window %02d for Trial %02d\n',mouse_id, exp_session, iWn,iTr)
    end
    saveName = sprintf('GC6f_emx_%02d_%2s_spont_VG_windowLen%03d_winStep_%03d_trial%02d_v2', mouse_id, exp_session, w_ln, w_st, iTr);
    
    if ispc
        save(['C:\Users\Li_Lab537\Dropbox\projects\calcium\VG_result\',saveName], 'VG', 'VG_lab','tidx_ca_di','tidx_whik');
    else
        save(['/Users/lizhu/Dropbox/projects/calcium/VG_result/',saveName], 'VG', 'VG_lab','tidx_ca_di','tidx_whik');
    end
end
