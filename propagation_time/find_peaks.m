%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%% find_peaks: find the peak of signals
%%
%% - Input:
%%
%%
%% - Output:
%%
%%
%% example:
%%
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [peak_idx] = find_peaks(ts, frame_size, opt)
    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results


    %% --------------------
    %% Constant
    %% --------------------


    %% --------------------
    %% Variable
    %% --------------------
    peak_idx = [];


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, error('not enough input'); end
    if nargin < 2, opt = ''; end

    [method, thresh] = get_opt(opt);


    %% --------------------
    %% Main starts
    %% --------------------
    if strcmp(method, 'max')
        peak_idx = find_peaks_max(ts, frame_size, thresh);
    elseif strcmp(method, 'first')
        peak_idx = find_peaks_first(ts, frame_size, thresh);
    elseif strcmp(method, 'max2')
        peak_idx = find_peaks_max2(ts, frame_size, thresh);
    elseif strcmp(method, 'first2')
        peak_idx = find_peaks_first2(ts, frame_size, thresh);
    else

    end

    
end


function [method, thresh] = get_opt(opt)
    method = 'max';
    thresh = 0.004;

    opts = regexp(opt, ',', 'split');
    for this_opt = opts
        eval([char(this_opt) ';']);
    end
end



%% =====================================================
%% find_peaks_max: find the max in each frame as peaks
function [peak_idx] = find_peaks_max(ts, frame_size, thresh)
    si = 1;
    peak_idx = [];
    
    while(si < length(ts)) 
        ts_seg = ts(si:min(end,si+frame_size-1));

        %% make sure there is a peak in this segment
        % max(ts_seg) - min(ts_seg)
        if max(ts_seg) - min(ts_seg) < thresh
            si = si + frame_size/2;
            continue;
        end
        
        [val, idx] = max(ts_seg);
        peak_idx = [peak_idx, si + idx - 1];

        si = si + idx - 1 + int32(floor(frame_size/2));
    end
end


%% =====================================================
%% find_peaks_fisrt: find the first one > thresh in the frame
function [peak_idx] = find_peaks_first(ts, frame_size, thresh)
    si = 1;
    peak_idx = [];
    

    %% find the first peak
    tmp_idx = find(ts > thresh);
    peak_idx = [tmp_idx(1)];


    %% find peaks every ~frame_size
    si = peak_idx + frame_size - 300;
    while(si < length(ts)) 
        ts_seg = ts(si:min(end,si+frame_size-1));

        tmp_idx = find(ts_seg > thresh);
        if length(tmp_idx) < 1
            si = si + frame_size;
        else
            idx = tmp_idx(1);
            peak_idx = [peak_idx, si + idx - 1];
            si = si + idx - 1 + frame_size - 300;
        end
    end
end




%% =====================================================
%% find_peaks_max: find the max in each frame as peaks
function [peak_idx] = find_peaks_max2(ts, frame_size, thresh)
    si = 1;
    peak_idx = [];
    
    while(si < length(ts)) 
        ts_seg = ts(si:min(end,si+frame_size-1));

        %% make sure there is a peak in this segment
        % max(ts_seg) - min(ts_seg)
        if max(ts_seg) - min(ts_seg) < thresh
            peak_idx = [peak_idx, 1];
            si = si + frame_size/2;
            continue;
        end
        
        [val, idx] = max(ts_seg);
        peak_idx = [peak_idx, si + idx - 1];

        si = si + idx - 1 + int32(floor(frame_size/2));
    end
end


%% =====================================================
%% find_peaks_fisrt: find the first one > thresh in the frame
function [peak_idx] = find_peaks_first2(ts, frame_size, thresh)
    si = 1;
    peak_idx = [];
    

    %% find the first 5 peaks
    num_calibration = 5;
    tmp_idx = find(ts > thresh);
    frame_stds = [tmp_idx(1)-300:frame_size:length(ts)];
    
    for fi = 1:length(frame_stds)
        ts_seg =ts(frame_stds(fi):min(end,frame_stds(fi)+frame_size-1));

        tmp_idx = find(ts_seg > thresh);
        if length(tmp_idx) > 0
            peak_idx = [peak_idx, frame_stds(fi) + tmp_idx(1) - 1];
        else
            peak_idx = [peak_idx, 1];
        end
    end

end
