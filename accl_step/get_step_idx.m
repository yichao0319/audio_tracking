%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% - Input:
%%
%% - Output:
%%
%%
%% example:
%%
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [step_idx, periodic_idx, autocorr, fig_idx] = get_step_idx(ts, win_ranges, fig_idx)
    % addpath('../utils');
    
    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results

    STFT = 0;


    %% --------------------
    %% Constant
    %% --------------------


    %% --------------------
    %% Variable
    %% --------------------
    

    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 2, win_ranges = [15:30]; end
    if nargin < 3, fig_idx = 0; end
    if length(win_ranges) < 2, win_ranges = [15:30]; end


    %% --------------------
    %% Main starts
    %% --------------------
    
    %% --------------------
    %% Find Periodic Pattern
    %% --------------------
    if DEBUG2, fprintf('Find Periodic Pattern\n'); end

    ts_len = length(ts);
    ts_r = [];
    ts_win = [];
    ti = 1;
    while(ti <= ts_len)
        r = [];
        for wi = 1:length(win_ranges)
            win = win_ranges(wi);
            if ti+2*win-1 > ts_len
                break;
            end

            corr = corrcoef(ts(ti:ti+win-1), ts(ti+win:ti+2*win-1));
            if numel(corr) > 1
                r(wi) = corr(1,2);
            end
        end

        %% find the best period length
        [max_r,idx] = max(r);
        max_win = win_ranges(idx);
        ts_r(ti:ti+max_win-1) = max_r;
        ts_win(ti:ti+max_win-1) = max_win;
    
        %% update per iteration    
        win_ranges = [max(1,max_win-5):max(15,max_win+5)];
        ti = ti + max_win;
    end

    % thresh = (mean(ts_r(1:10)) + max(ts_r)) * 1 / 2;
    thresh = 0.8;
    fprintf('  thresh=%f, min=%d, max=%d\n', thresh, min(ts_r), max(ts_r));

    %% thresholding to find periodic portions
    [v,std_idx] = find(ts_r > thresh);
    std_idx = std_idx(1);
    tmp_ts_r = ts_r(std_idx:end);
    [v,end_idx] = find(tmp_ts_r <= thresh);
    end_idx = end_idx(1) - 1;
    end_idx = end_idx + std_idx - 1;

    % std_idx = 100;
    % end_idx = 700;
    periodic_idx = [std_idx, end_idx];
    autocorr = ts_r;


    %% plot ts and corrcoef
    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % subplot(2,1,1);
    % plot(ts);
    % hold on;
    % plot([std_idx, end_idx], ts([std_idx, end_idx]), 'ro');
    % set(gca, 'XLim', [1 ts_len]);
    % subplot(2,1,2);
    % plot(ts_r);
    % hold on;
    % plot([std_idx, end_idx], ts_r([std_idx, end_idx]), 'ro');
    % set(gca, 'XLim', [1 ts_len]);


    %% --------------------
    %% Find Steps
    %% --------------------
    if DEBUG2, fprintf('Find Steps\n'); end

    tmp = ts(std_idx:end_idx);
    [xmax,imax,xmin,imin] = extrema(tmp);
    imax = sort(imax + std_idx - 1);
    imin = sort(imin + std_idx - 1);

    step_idx = [];
    ti = std_idx;
    while(ti <= end_idx - ts_win(end));
        win = ts_win(ti);
        % imin
        % fprintf('  range = %d-%d\n', ti, ti+win);
        ii = find(imin > ti & imin < ti + 1.5*win);
        % ii
        % ts(imin(ii))
        [v,idx] = min(ts(imin(ii)));

        step_idx = [step_idx, imin(ii(idx))];

        % fprintf('  > ti=%d: step=%d\n', ti, step_idx(end));
        ti = step_idx(end) + 1;
        % input('')
    end

end






