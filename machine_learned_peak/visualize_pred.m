%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% example:
%%   visualize_pred('rx.3.1')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize_pred(filename)
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
    spike_dir  = './spikes/';
    pred_dir   = './svm_models/';
    output_dir = '';

    sample        = 5;
    num_neighbors = 5;


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;
    fs = 44100;


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, filename = 'rx.3.1'; end


    %% --------------------
    %% Main starts
    %% --------------------

    %% --------------------
    %% Load Data
    %% --------------------
    if DEBUG2, fprintf('Read data\n'); end

    pred = load(sprintf('%s%s.samp%d.neighbor%d.txt.predict', pred_dir, filename, sample, num_neighbors));

    corr   = load(sprintf('%s%s.corr.txt', spike_dir, filename));
    status = load(sprintf('%s%s.status.txt', spike_dir, filename));
    peaks  = load(sprintf('%s%s.gt.txt', spike_dir, filename));

    peaks = peaks(1:length(status));
    itvl = mean([peaks(2:end) - peaks(1)] ./ [1:(length(peaks)-1)]');
    row  = 1;


    if DEBUG3,
        fprintf('  corr size: %dx%d\n', size(corr));
        fprintf('  status size: %dx%d\n', size(status));
        fprintf('  peaks size: %dx%d\n', size(peaks));
        fprintf('  pred size: %dx%d\n', size(pred));

        fprintf('  itvl: mean=%.2f\n', itvl);
    end


    %% --------------------
    %% Sync data
    %% --------------------
    if DEBUG2, fprintf('Sync Data\n'); end

    range_idx = [];
    for si = 1:length(status)
        peak_idx = peaks(si);
        std_idx = peak_idx - 200;
        end_idx = peak_idx + 200;
        range_idx = [range_idx unique(sort([std_idx:sample:end_idx, peak_idx]))];
    end
    range_idx = range_idx';


    pred_idx = range_idx(pred == 1);
    pred_corr = corr(pred_idx);

    if DEBUG3,
        fprintf('  range size: %dx%d\n', size(range_idx));
        fprintf('  pred peak size: %dx%d\n', size(pred_idx));
    end



    %% --------------------
    %% Plot data
    %% --------------------
    if DEBUG2, fprintf('Plot data\n'); end

    fig_idx = fig_idx + 1;
    fh = figure(fig_idx); clf;

    idx = 1:length(corr);
    plot(idx/fs, corr, '-k');
    hold on;
    lh = plot(peaks/fs, corr(peaks), 'ro');
    set(lh, 'MarkerSize', 10);
    lh = plot(pred_idx/fs, pred_corr, 'bx');
    set(lh, 'MarkerSize', 10);

    waitforbuttonpress


    for si = 1:length(status)
        fh = figure(11); clf;

        peak_idx = peaks(si);
        idx = 1:length(corr);
        plot(idx/fs, corr, '-k.');
        hold on;
        lh = plot(peak_idx/fs, corr(peak_idx), 'ro');
        set(lh, 'MarkerSize', 10);
        lh = plot(pred_idx/fs, pred_corr, 'bx');
        set(lh, 'MarkerSize', 10);

        xlim([peak_idx-200 peak_idx+200]/fs);
        title(sprintf('signal%d: status=%d', si, status(si)));

        waitforbuttonpress
    end
end