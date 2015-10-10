%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%% MAIN
%% - Input:
%%
%%
%% - Output:
%%
%%
%% example:
%%   cal_t1_itvl('0903.exp4.fix.pn2047.1m', 8000, 2047, 1, 1)
%%   cal_t1_itvl('0903.exp5.fix.pn2047.0.5m', 8000, 2047, 0.5, 1)
%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [base, t1, itvl] = cal_t1_itvl(peaks_idx, init_dist, Fs, config_itvl)
    % addpath('../utils');
    
    %% --------------------
    %% DEBUG
    %% --------------------
    DEBUG0 = 0;
    DEBUG1 = 1;
    DEBUG2 = 1;  %% progress
    DEBUG3 = 1;  %% verbose
    DEBUG4 = 1;  %% results

    warning('off','comm:commsrc:pn:GenPolyNotPrimitive');


    %% --------------------
    %% Constant
    %% --------------------
    fig_dir = './fig/';

    sound_speed = 331;

    thresh = 0.8;

    fig_idx = 0;
    font_size = 16;




    %% --------------------
    %% Variable
    %% --------------------
    


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, error('not enough input param'); end
    

    %% --------------------
    %% Main starts
    %% --------------------
    
    %% ====================================
    %% Find stable intervals
    %% ====================================
    scale = 10;
    base_confident = [];
    all_itvls = [];
    for base = 1:length(peaks_idx)
        selections = setxor([1:length(peaks_idx)], base);
        itvls = int32((peaks_idx(selections) - peaks_idx(base)) ./ (selections-base) * scale);

        all_itvls = [all_itvls itvls];

        ranges = [min(itvls):max(itvls)];
        counts = histc(itvls, ranges);
        idx = find(counts > 0);

        fprintf('----------------\nbase = %d\n', base);
        for i = 1:length(idx)
            fprintf('  %d: %d\n', ranges(idx(i)), counts(idx(i)));
        end

        [v,idx] = max(counts);
        best_itvl = ranges(idx);
        idx = find(itvls == best_itvl);
        best_itvl = double(best_itvl) / scale;
        % stable_idx = sort([base selections(idx)]);
        stable_idx = selections(idx);
        base_confident(base) = v;
        base_selections{base} = selections(idx);
        base_itvl(base) = best_itvl;

        % tb = peaks_idx(stable_idx)/Fs - (stable_idx-base)*best_itvl/Fs
        % dists = [0:1:10];
        % for di = 1:length(dists)
        %     this_tb = tb - dists(di) / sound_speed;
        %     errs(di) = sum(this_tb - this_tb(1));
        %     fprintf('  dist=%f: err=%.10f\n', dists(di), errs(di));
        % end
        % [v,idx] = sort(errs);

        % input('')
    end

    %% select the most common interval as the interval
    ranges = [min(all_itvls):max(all_itvls)];
    counts = histc(all_itvls, ranges);

    fprintf('----------------\nall\n', base);
    idx = find(counts > 0);
    for i = 1:length(idx)
        fprintf('  %d: %d\n', ranges(idx(i)), counts(idx(i)));
    end

    %% select the base with the same interval and the highest confident
    [v,idx] = max(counts);
    % ranges(idx)
    best_itvl = double(ranges(idx)) / scale;
    idx = find(base_itvl == best_itvl);
    if length(idx) > 0
        [v,iidx] = max(base_confident(idx));
        base = idx(iidx);
    else
        [v,base] = max(base_confident);
    end


    if abs(base_itvl(base)-config_itvl*Fs) > config_itvl * 50;
        itvl = config_itvl;
    else
        itvl = base_itvl(base) / Fs;
    end
    fprintf('  itvl = %d (use %d)\n', base_itvl(base), itvl*Fs);
    
    
    tbs = [];
    %% t_b = t_i - (i-b)*itvl - d_i / v
    for si = 1:length(base_selections{base})
        idx = base_selections{base}(si);
        tbs(si) = peaks_idx(idx)/Fs - (idx-base)*itvl - init_dist/sound_speed;
        
        fprintf('  base%d - peak%d: t1 = %.10f\n', base, idx, tbs(si));
    end

    t1 = mean(tbs);
    
end

