function [r] = my_xcorr(ts1, ts2)
    ts2_len = length(ts2);
    r = [];
    for ti = 1:length(ts1)-ts2_len+1
        ts1_seg = ts1(ti:ti+ts2_len-1);

        rr = corrcoef(ts1_seg, ts2);
        r(ti) = rr(1,2);
    end
end
