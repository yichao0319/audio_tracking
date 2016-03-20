%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Yi-Chao Chen @ UT Austin
%%
%% example:
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function preprocess_svm_input()

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
    output_dir = './processed/';


    %% --------------------
    %% Variable
    %% --------------------
    fig_idx = 0;
    sample = 5;
    num_neighbors = 5;


    %% --------------------
    %% Check input
    %% --------------------
    if nargin < 1, arg = 1; end


    %% --------------------
    %% Main starts
    %% --------------------



    for ei = 1:23
        filename = sprintf('rx.3.%d', ei);
        fprintf('  file=%s\n', filename);


        fileID = fopen(sprintf('%s%s.samp%d.neighbor%d.txt', output_dir, filename, sample, num_neighbors), 'w');
        corr   = load(sprintf('%s%s.corr.txt', input_dir, filename));
        status = load(sprintf('%s%s.status.txt', input_dir, filename));
        peaks  = load(sprintf('%s%s.gt.txt', input_dir, filename));

        itvl = mean([peaks(2:end) - peaks(1)] ./ [1:(length(peaks)-1)]');
        row  = 1;

        for si = 1:length(status)
            peak_idx = peaks(si);
            std_idx = peak_idx - 200;
            end_idx = peak_idx + 200;
            range_idx = unique(sort([std_idx:sample:end_idx, peak_idx]));

            for ri = 1:length(range_idx)
                idx = range_idx(ri);
                if idx == peak_idx
                    fprintf(fileID, '1 ');
                    % fprintf('  > peak\n');
                else
                    fprintf(fileID, '0 ');
                end

                std_idx2 = idx - 200;
                end_idx2 = idx + 200;
                range_idx2 = [std_idx2:sample:end_idx2];
                range_idx2 = setdiff(range_idx2, idx);
                invalid_range = find(range_idx2 < 1);
                valid_range   = find(range_idx2 >= 1);

                features(row, 1) = corr(idx);
                col = 1;
                if length(invalid_range) > 0
                    features(row, col+invalid_range) = 0;
                end
                features(row, col+valid_range) = corr(range_idx2(valid_range));
                col = col + length(range_idx2);
                % fprintf('  row=%d, col=%d\n', row, col);


                sel_neighbors = [-num_neighbors:-1 1:num_neighbors];
                sel_neighbors = round(sel_neighbors*itvl + idx);
                sel_neighbors = sel_neighbors(sel_neighbors>0);
                sel_neighbors = sel_neighbors(1:num_neighbors);

                features(row, col+[1:num_neighbors]) = corr(sel_neighbors);
                col = col + num_neighbors;
                % fprintf('     col=%d\n', col);


                for ci = 1:col
                    fprintf(fileID, '%d:%d ', ci, features(row,ci));
                end
                fprintf(fileID, '\n');

                row = row + 1;
            end
        end

        fclose(fileID);
    end

end