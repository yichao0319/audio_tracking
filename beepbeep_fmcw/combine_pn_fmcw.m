
function [y_comb] = combine_pn_fmcw(y_pn, y_fmcw, method)
    if strcmp(method, 'avg')
        y_comb = (y_pn+y_fmcw)/2;
    elseif strcmp(method, 'mrc_var')
        var_pn = std(y_pn) ^ 2;
        var_fmcw = std(y_fmcw) ^ 2;
        y_comb = (var_fmcw * y_pn + var_pn * y_fmcw) / (var_pn+var_fmcw);

    elseif strcmp(method, 'mrc_spk')
        thresh = 0.6;

        [pks, locs] = findpeaks(y_pn);
        var_pn = length(find(pks > (max(y_pn)*thresh)));

        [pks, locs] = findpeaks(y_fmcw);
        var_fmcw = length(find(pks > (max(y_fmcw)*thresh)));

        y_comb = (var_fmcw * y_pn + var_pn * y_fmcw) / (var_pn+var_fmcw);

    elseif strcmp(method, 'mrc_spk2')
        thresh = 0.6;

        [pks, locs] = findpeaks(y_pn);
        var_pn = length(find(pks > (max(y_pn)*thresh))) ^ 2;

        [pks, locs] = findpeaks(y_fmcw);
        var_fmcw = length(find(pks > (max(y_fmcw)*thresh))) ^ 2;

        y_comb = (var_fmcw * y_pn + var_pn * y_fmcw) / (var_pn+var_fmcw);

    end
end

