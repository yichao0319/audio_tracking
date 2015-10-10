%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
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

function visualize_combine()
    % addpath('../utils');
    
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
    input_dir  = '';
    fig_dir = './fig/';

    thresh = 0.1;  %% torrlerate error < thresh. i.e. when error < thresh, error = 0


    %% --------------------
    %% Variable
    %% --------------------
    filename = '0904.acc';
    fig_idx = 0;

    dists = [sqrt(2.8^2+1^2), 1, sqrt(3.1^2+1^2), sqrt(3.1^2+2^2), 2, sqrt(2.8^2+2^2), sqrt(2.8^2+1^2)] * 1.02;
    % displacements = [2.8, 3.1, 3, 3.1, 2.8, 3] * 1.02;
    % directions = [0,1,0; 0,1,0; 1,0,0; 0,-1,0; 0,-1,0; -1,0,0];

    %% pn, mean
    % dists = [3.0775, 1.0023, 3.2573, 3.6892, 2.0023, 2.3526, 3.0775];
    %% pn, median
    % dists = [1.5291, 1.0016, 3.0369, 3.5316, 1.9073, 2.5240, 1.5291];
    %% sinc, mean
    % dists = [2.9736, 1.7843, 3.4112, 10.8496, 2.0038, 3.4409, 2.9736];
    %% sinc, median
    % dists = [2.9735, 1.7343, -5.4669, 9.4542, 2.0038, 3.4428, 2.9735];
    displacements = [5, 6, 5, 5, 5, 6] * 0.7;
    directions = [0.7739, -0.6298, 0.0662; 0.7701, -0.6304, 0.0974; 0.6124, 0.7828, 0.1103; -0.7809, 0.6166, -0.0995; -0.7889, 0.6060, -0.1022; -0.5324, -0.8447, -0.0544];

    loc{1} = [0,0];

    gran = 0.01;
    
    
    %% --------------------
    %% Check input
    %% --------------------
    % if nargin < 1, arg = 1; end


    %% --------------------
    %% Main starts
    %% --------------------
    fh = figure(1); clf;
    
    for si = 1:length(dists)
        if si > 1
            % loc{si} = loc{si-1} + displacements(si-1) * [sin(directions(si-1)), cos(directions(si-1))];
            loc{si} = loc{si-1} + displacements(si-1)*directions(si-1,1:2);
        end

        fh = figure(1)
        plot(loc{si}(1), loc{si}(2), 'ro'); 
        hold on;
        ang=0:0.01:2*pi; 
        xp=dists(si)*cos(ang);
        yp=dists(si)*sin(ang);
        plot(loc{si}(1)+xp,loc{si}(2)+yp);
        if si > 1, arrow(loc{si-1}, loc{si}); end
        print(fh, '-dpsc', [fig_dir filename '.' num2str(si) '.line.eps']);


        %% recalculate probability
        xlim = get(gca,'xlim');
        ylim = get(gca,'ylim');

        xrange = xlim
        yrange = ylim
        prob_map_x = [xrange(1):gran:xrange(2)];
        prob_map_y = [yrange(1):gran:yrange(2)];
        prob_map = zeros(length(prob_map_y), length(prob_map_x));
        % return;

        fh = figure(2); clf;
        for sj = 1:si
            distx = (prob_map_x - loc{sj}(1));
            disty = (prob_map_y - loc{sj}(2))';
            distxy = sqrt(repmat(distx .^ 2, length(disty), 1) + repmat(disty .^ 2, 1, length(distx)));
            dist_dif = abs(distxy - dists(sj));
            idx = find(dist_dif < thresh);
            dist_dif(idx) = thresh;
            prob_map = prob_map - log(dist_dif);

            norm_prob_map = prob_map - min(min(prob_map));
            norm_prob_map = norm_prob_map / sum(sum(norm_prob_map));
            x = repmat(prob_map_x, length(prob_map_y), 1);
            y = repmat(prob_map_y', 1, length(prob_map_x));
            est_loc{sj}(1) = sum(sum(norm_prob_map .* x));
            est_loc{sj}(2) = sum(sum(norm_prob_map .* y));
            % est_loc{sj}
        end

        imagesc(prob_map_x, prob_map_y, norm_prob_map); % frequency-time Plot of the signal
        colorbar;
        set(gca,'YDir','normal');
        hold on;

        for sj = 1:si
            plot(loc{sj}(1), loc{sj}(2), 'ro'); 
            hold on;
            if sj > 1, arrow(loc{sj-1}, loc{sj}); end
        end
        print(fh, '-dpsc', [fig_dir filename '.' num2str(si) '.map.eps']);

        % waitforbuttonpress
    end
    
    fig_idx = 2;
end