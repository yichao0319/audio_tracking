function [idx_offset] = sync_pn_by_fmcw(y, y_pn, y_fmcw, corr, idx, win, tx_chirp, new_chirp_len)
    fs = 44100;
    chirp_len = length(tx_chirp);


    %% PN
    [~,pn_idx] = max(corr);

    %% FMCW
    for ci = 1:chirp_len
        tmp = y_fmcw(ci:ci+chirp_len-1);
        r = corrcoef(tmp, tx_chirp);
        corr_fmcw(ci) = r(1,2);
    end
    [~,fmcw_idx] = max(corr_fmcw);
    num_chirps = floor((length(y)-fmcw_idx+1) / new_chirp_len);
    fmcw_idxs = ones(1,num_chirps) * fmcw_idx + [0:num_chirps-1]*new_chirp_len;

    dif = pn_idx-fmcw_idxs;
    dif = dif(dif>0);
    idx_offset = min(dif);
    fprintf('  offset = %.2f\n', idx_offset);


    fh = figure(31); clf;
    spectrogram(y(idx:idx+2*fs-1),256,64,256,fs,'yaxis');
    hold on;
    plot(ones(1,2)*(pn_idx-idx+1)/fs, get(gca, 'ylim'), '-r');
    for ci = 1:1:num_chirps
        plot(ones(1,2)*(fmcw_idxs(ci)-idx+1)/fs, get(gca, 'ylim'), '-g');
    end
    xlim([0,2]);

end