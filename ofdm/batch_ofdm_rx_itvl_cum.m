%% sum up different number of signals to see the improvement 
function batch_ofdm_rx_itvl_cum
    font_size = 18;
    % audio_lens = [3.4:1:40.4];
    audio_lens = [3.4:1:20.4];
    
    for li = 1:length(audio_lens)
        [orig(li), cum(li), err_orig_max(li), err_orig_median(li), err_orig_avg(li), err_cum(li)] = ofdm_rx_itvl_cum('1111.exp2', 128, 1, audio_lens(li));
    end

    fh = figure(1); clf;
    bh = bar([orig; cum]');
    set(bh(1), 'FaceColor', 'b');
    set(bh(2), 'FaceColor', 'r');
    set(gca, 'FontSize', font_size);
    legend('Orig', 'Sum', 'Location', 'NorthOutside');
    xlabel('# Accumulated Signals', 'FontSize', font_size);
    ylabel('Corr Magnitude', 'FontSize', font_size);
    set(gca, 'XLim', [0.5 length(audio_lens)+0.5]);

    fh = figure(2); clf;
    bh = bar([err_orig_max; err_orig_median; err_orig_avg; err_cum]'+1);
    set(bh(1), 'FaceColor', 'b');
    set(bh(2), 'FaceColor', 'g');
    set(bh(3), 'FaceColor', 'y');
    set(bh(4), 'FaceColor', 'r');
    set(gca, 'FontSize', font_size);
    legend('Orig-max', 'Orig-median', 'Orig-avg', 'Sum', 'Location', 'NorthOutside');
    xlabel('# Accumulated Signals', 'FontSize', font_size);
    ylabel('Error (samples)', 'FontSize', font_size);
    % set(gca, 'XLim', [0.5 length(audio_lens)+0.5]);

    orig
    cum
end