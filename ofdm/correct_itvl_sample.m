%% Because of the time skew is different among different machines, 
%% the received signal interval is not necessary the same as the interval in tx machine.
%% The function find the interval which minimize the error
function [itvl_samp, base_sig] = correct_itvl_sample(peak_idx)
    for base_sig = 1:length(peak_idx)
        est_itvls = (peak_idx - peak_idx(base_sig)) ./ ((1:length(peak_idx)) - base_sig);
        est_itvls = est_itvls(find(est_itvls<Inf));

        itvl_samps(base_sig) = median(est_itvls);
        est_peak_idx = round(peak_idx(base_sig) + ([1:length(peak_idx)]-base_sig)*itvl_samps(base_sig));
        errs(base_sig) = mean(abs(peak_idx - est_peak_idx));

        fprintf('base%d: ', base_sig);
        % fprintf('%4.2f,', est_itvls);
        fprintf('err=%.2f, itvl=%4.2f', errs(base_sig), itvl_samps(base_sig));
        fprintf('\n');

    end

    [v,base_sig] = min(errs);
    itvl_samp = itvl_samps(base_sig);
end
