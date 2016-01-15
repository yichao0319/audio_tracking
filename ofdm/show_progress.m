%% show_progress
function show_progress(cur_idx, len, std_idx)
    persistent progress;

    if cur_idx == std_idx
        progress = ones(1,10);
    end

    state = (cur_idx-1)/len;
    idx = find(progress > 0);
    if length(idx) > 1
        if state*10 >= idx(1)
            fprintf('  %d%%\n', ceil(state*100));
            progress(idx(1)) = 0;
        end
    end
end
