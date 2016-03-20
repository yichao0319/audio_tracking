function [maxIndex, corr2] = findStartIndexByDoubleWin(y,windowSize,detectLength)

corr2 = zeros(1, detectLength-2*windowSize);
y=abs(y);
maxCorr=0;
maxIndex=0;
denom=sum(y(1:windowSize));
nom=sum(y(windowSize+2:2*windowSize+1));
for i=windowSize+1:detectLength-windowSize-1
    corr2(i-windowSize) = nom/denom;;
    corr=nom/denom;
    denom=denom+y(i)-y(i-windowSize);
    nom=nom+y(i+windowSize+1)-y(i+1);
    if corr>maxCorr
        maxCorr=corr;
        maxIndex=i-windowSize;
    end
end

