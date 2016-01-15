function [new_analogData] = correct_freq_offset(analogData, gt_idx, preamble)
    % len = length(preamble);
    % for si = 1:length(gt_idx)
    %     seg = analogData(gt_idx(si):gt_idx(si)+len-1);
    %     theta(si) = angle(seg * preamble') / (2*pi);
    % end
    % theta
    % % freq_offset = mean(theta) / 128;
    % freq_offset = mean(theta)

    preamble_len = length(preamble);
    
    %% ----------------------
    %% Find H from the first rx preamble
    % rx = analogData(gt_idx(1):gt_idx(1)+preamble_len-1);

    % H = rx*(preamble')*pinv(preamble*(preamble'));
    % new_analogData = H*analogData;

    %% ----------------------
    %% Find H from each rx preamble
    new_analogData = analogData;
    for si = 1:length(gt_idx)
        rx = analogData(gt_idx(si):gt_idx(si)+preamble_len-1);
        H = rx*(preamble')*pinv(preamble*(preamble'));
        new_analogData(gt_idx(si):gt_idx(si)+preamble_len-1) = H*rx;
    end
end