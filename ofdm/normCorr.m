function dist=normCorr(corrHeight,v)

% correct the norm by considering the difference of the height
dist=sqrt(corrHeight^2+sum(v.^2));