%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% example:
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize_spikes(filename)
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
    input_dir  = './spikes/';
    output_dir = '';


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

    corr = load(sprintf('%s%s.corr.txt', input_dir, filename));
    status = load(sprintf('%s%s.status.txt', input_dir, filename));
    peaks = load(sprintf('%s%s.gt.txt', input_dir, filename));


    if DEBUG3,
        fprintf('  corr size: %dx%d\n', size(corr));
        fprintf('  status size: %dx%d\n', size(status));
        fprintf('  peaks size: %dx%d\n', size(peaks));

        % itvls = peaks(2:end) - peaks(1:end-1);
        itvls = [peaks(2:end) - peaks(1)] ./ [1:(length(peaks)-1)]';
        fprintf('  itvl: mean=%.2f, median=%.2f, min=%2.f, max=%.2f\n', mean(itvls), median(itvls), min(itvls), max(itvls));
    end


    %% --------------------
    %% Plot data
    %% --------------------
    if DEBUG2, fprintf('Plot data\n'); end

    for si = 1:length(status)
        fh = figure(11); clf;

        peak_idx = peaks(si);
        idx = 1:length(corr);
        plot(idx/fs, corr, '-b.');
        hold on;
        plot(peak_idx/fs, corr(peak_idx), 'ro');
        xlim([peak_idx-200 peak_idx+200]/fs);

        title(sprintf('signal%d: status=%d', si, status(si)));

        waitforbuttonpress
    end
end