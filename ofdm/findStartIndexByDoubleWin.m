function maxIndex=findStartIndexByDoubleWin(y,windowSize,detectLength)

% corr=zeros(detectLength-2*windowSize);
y=abs(y);
maxCorr=0;
maxIndex=0;
denom=sum(y(1:windowSize));
nom=sum(y(windowSize+2:2*windowSize+1));
for i=windowSize+1:detectLength-windowSize
    % corr(i-windowSize)=sum(abs(y(i+1:i+windowSize)))/sum(abs(y(i-windowSize:i-1)));
    corr=nom/denom;
    denom=denom+y(i)-y(i-windowSize);
    nom=nom+y(i+windowSize+1)-y(i+1);
    if corr>maxCorr
        maxCorr=corr;
        maxIndex=i-windowSize;
    end
end
% plot(corr);


% function maxIndex=findStartIndexByDoubleWin(y,windowSize,detectLength)

% fprintf('  start find index\n');

% % corr=zeros(detectLength-2*windowSize, 1);
% y=abs(y);
% maxCorr=0;
% maxIndex=0;
% % denom=sum(y(1:windowSize));
% % nom=sum(y(windowSize+2:2*windowSize+1));

% for i=windowSize+1:windowSize:detectLength-windowSize
%     corr=sum(abs(y(i+1:i+windowSize)))/sum(abs(y(i-windowSize:i-1)));
%     if corr>maxCorr
%         maxCorr=corr;
%         maxIndex=i-windowSize;
%     end
% end

% for i=max(1+windowSize,maxIndex):maxIndex+1*windowSize
%     corr=sum(abs(y(i+1:i+windowSize)))/sum(abs(y(i-windowSize:i-1)));
%     if corr>maxCorr
%         maxCorr=corr;
%         maxIndex=i-windowSize;
%     end
% end
% % plot(corr);

