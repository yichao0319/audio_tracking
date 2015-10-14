function sol=trainSGD(anchor,x,y,c,N,n)

z=0.1*rand([1,3*n+1]);
z=[5 1 0 1];
alpha=0.01;
k=1;
p=numel(x);
for i=1:N
    d=[x(k),y(k),c(k)];
    z=z-alpha*calcGradient(z,anchor,n,d);
    k=k+1;
    if k>p
        k=1;
    end
end
sol=z;