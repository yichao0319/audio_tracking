function [ys, idx]=split(y,index,span,interval,N)

ys=zeros(N,span);
idx = [];
for i=1:N
    idx_std = index+(i-1)*(span+interval);
    idx_end = index+span-1+(i-1)*(span+interval);
    if idx_std > length(y) | idx_end > length(y)
        ys = ys(1:i-1,:);
        break;
    end

    idx = [idx idx_std];
    ys(i,:)=y(idx_std:idx_end);
end