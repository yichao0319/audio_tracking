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
    
    fig_idx = 0;

    %% exp
    % filename = '0904.acc';
    % dists = [sqrt(2.8^2+1^2), 1, sqrt(3.1^2+1^2), sqrt(3.1^2+2^2), 2, sqrt(2.8^2+2^2), sqrt(2.8^2+1^2)] * 1.02;
    % displacements = [5, 6, 5, 5, 5, 6] * 0.7;
    % directions = [0.7739, -0.6298, 0.0662; 0.7701, -0.6304, 0.0974; 0.6124, 0.7828, 0.1103; -0.7809, 0.6166, -0.0995; -0.7889, 0.6060, -0.1022; -0.5324, -0.8447, -0.0544];

    %% 1023.exp8
    % filename = '1023.exp8';
    % dists = sqrt(sum(([0,0; 10,0; 20,0; 20,10; 20,20; 10,20; 0,20; 0,10; 0,0] - repmat([10,10],9,1)).^2,2)) .* (rand(9,1)/10+1);
    % displacements = ones(1,8)*13 * 0.7; 
    % displacements([2,5]) = 12; displacements([3]) = 13.2;
    % directions = [1.00,0.00;1.00,-0.01;0.08,-1.00;0.09,-1.00;-1.00,-0.09;-1.00,-0.09;-0.15,0.99;-0.16,0.99];
    %% 1023.exp9
    % filename = '1023.exp9';
    % dists = sqrt(sum(([0,0; 10,0; 20,0; 30,0; 30,10; 30,20; 20,20; 10,20; 0,20; 0,10; 0,0] - repmat([10,10],11,1)).^2,2)) .* (rand(11,1)/1+1);
    % displacements = ones(1,10)*13 * 0.7; 
    % displacements([2,5]) = 12; displacements([3]) = 13.2;
    % directions = [1.00,0.00;1.00,-0.05;1.00,-0.07;-0.10,-0.99;-0.19,-0.98;-0.90,-0.44;-0.87,-0.49;-0.86,-0.50;0.55,0.83;0.59,0.81;];
    filename = '1023.exp9';
    dists = sqrt(sum(([0,0; 10,0; 20,0; 30,0; 30,10; 30,20; 20,20; 10,20; 0,20; 0,10; 0,0] - repmat([10,10],11,1)).^2,2)) .* (rand(11,1)/0.8+1);
    displacements = ones(1,10)*13 * 0.7; 
    displacements([2,5]) = 12 * 0.7; displacements([3]) = 13.2 * 0.7;
    directions = [1.00,0.00;1.00,-0.05;1.00,-0.07;-0.10,-0.99;-0,-0.98;-0.90,-0;-0.87,-0;-0.86,-0.10;0,0.93;0.1,0.91;];



    loc{1} = [0,0];

    % gran = 0.01;
    gran = 0.5;
    
    
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
        % print(fh, '-dpsc', [fig_dir filename '.' num2str(si) '.line.eps']);
        print(fh, '-dpng', [fig_dir filename '.' num2str(si) '.line.png']);
        

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
        % print(fh, '-dpsc', [fig_dir filename '.' num2str(si) '.map.eps']);
        print(fh, '-dpng', [fig_dir filename '.' num2str(si) '.map.png']);

        % waitforbuttonpress
    end
    
    fig_idx = 2;
end