function sol=trainGD(anchor,x,y,c,N,n)

z=0.1*rand([1,3*n+1]);
alpha=0.05;
p=numel(x);
for i=1:N
    sumG=0;
    for j=1:p
        d=[x(j),y(j),c(j)];
        sumG=sumG+calcGradient(z,anchor,n,d);
    end
    meanG=sumG/p;
    z=z-alpha*meanG;
end
sol=z;