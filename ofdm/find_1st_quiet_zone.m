%% find_1st_quiet_zone:
%%   to aviod the signal is received at the beginning, find a quiet zone to start a audio file
function [std_idx] = find_1st_quiet_zone(analogData)
    win = 1000;
    avg_mag(1) = mean(abs(analogData(1:win)));
    for ti = 2:length(analogData)/10
        avg_mag(ti) = (avg_mag(ti-1) * win - abs(analogData(ti-1)) + abs(analogData(ti+win-1))) / win;
    end
    idx_low = find(avg_mag < mean(avg_mag));
    idx_low_itvl = idx_low(2:end) - idx_low(1:end-1);
    idx = find(idx_low_itvl > 100);
    if length(idx) > 0
        idx = idx(1);
        std_idx = ceil((idx_low(idx) + idx_low(1)) *1/3);
    else
        std_idx = ceil((idx_low(end) + idx_low(1)) *1/3);
    end

    % fig_idx = 0;
    % fig_idx = fig_idx + 1;
    % fh = figure(fig_idx); clf;
    % subplot(2,1,1);
    % plot(analogData(1:length(avg_mag)));
    % subplot(2,1,2);
    % plot(avg_mag); hold on;
    % plot([1,length(avg_mag)], [1,1]*mean(avg_mag), '-r');
    % plot(std_idx, analogData(std_idx), 'go');
end